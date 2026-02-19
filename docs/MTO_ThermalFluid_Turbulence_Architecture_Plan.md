# MTO_ThermalFluid: k-ω SST + Curvature Correction — Architecture Analysis & Implementation Plan

**Document Version:** 1.0  
**Date:** February 2025  
**Scope:** MTO_ThermalFluid module only  

---

## Executive Summary

This document provides the architecture analysis, mathematical formulation, and stepwise implementation plan for integrating a **fully consistent primal–adjoint k-ω SST turbulence model** (standard and curvature-corrected variants) into the MTO_ThermalFluid topology optimization framework. The implementation targets research-grade quality suitable for high-impact publication, HPC deployment, and optimization of turbulent microchannels with strong curvature and secondary flows.

---

## Part 1: Current MTO_ThermalFluid Implementation Analysis

### 1.1 Solver Structure Mapping

| Component | Location | Current State |
|-----------|----------|---------------|
| **Main loop** | `MTO_ThermalFluid.C` | SIMPLE iteration: `update.H` → Primal (U, T) → Adjoint (Tb, Ub, Ua) → `costfunction.H` → `sensitivity.H` |
| **Primal momentum** | `Primal_U.H` | Uses `turbulence->divDevReff(U)` + Brinkman `α(xh)U` |
| **Primal energy** | `Primal_T.H` | Advection–diffusion with effective diffusivity `DT(xh)`; **no turbulent thermal diffusivity (αt)** |
| **Adjoint momentum (power)** | `AdjointFlow_Ua.H` | Uses `fvm::laplacian(nu, Ua)` — **laminar viscosity only** |
| **Adjoint momentum (heat)** | `AdjointHeat_Ub.H` | Same: `fvm::laplacian(nu, Ub)` — **laminar viscosity only** |
| **Adjoint energy** | `AdjointHeat_Tb.H` | Adjoint of T equation (no turbulent diffusion) |
| **Turbulence** | `readTransportProperties.H` | `incompressible::turbulenceModel::New()` — currently **laminar** (`turbulenceProperties`: `laminar`) |
| **Sensitivity** | `sensitivity.H` | Volume-based: `gsensPowerDiss`, `fsenshMeanT` from U·Ua, U·Ub, grad(T)·grad(Tb); **no turbulence sensitivities** |
| **BCs** | Custom patch fields | `adjointOutletVelocityPower`, `adjointOutletPressurePower`, `adjointOutletVelocityHeat`, `adjointOutletPressureHeat` — use `nu` (laminar) |

### 1.2 Critical Inconsistencies (Must Fix)

1. **Adjoint viscous term mismatch**
   - Primal: `divDevReff(U)` → `div(ν_eff dev2(∇U))` with ν_eff = ν + ν_t when turbulence is on
   - Adjoint: `fvm::laplacian(nu, Ua)` — uses ν only, and Laplacian structure (not div-of-deviatoric)
   - **Consequence:** Adjoint is not the true transpose of the primal viscous operator; sensitivities will be wrong for any RAS model

2. **fvSchemes vs. source code**
   - `fvSchemes` defines `div((nuEff*dev2(T(grad(Ua)))))` for Ua
   - Source uses `fvm::laplacian(nu,Ua)` — the assembled matrix does not match fvSchemes
   - With laminar, ν_eff = ν and `div(ν dev2(∇Ua))` ≈ `ν laplacian(Ua)` for incompressible flow (∇·Ua=0)
   - With turbulence, ν_eff ≠ ν; adjoint must use ν_eff and include adjoint turbulence contributions

3. **Frozen turbulence assumption**
   - Adjoint BCs have commented-out `rasModel.nuEff()`; they use laminar `nu`
   - No adjoint turbulence variables (ka, ωa)
   - No differentiation of ν_t, production, or destruction terms

4. **Thermal turbulence**
   - Primal T equation uses `DT` only; no turbulent thermal diffusivity `αt = νt/Prt`
   - For turbulent flows, energy equation should include `α_eff = α + αt`
   - Adjoint Tb must be updated accordingly

5. **Brinkman–turbulence interaction**
   - Brinkman term α(xh)U is correctly in primal and adjoint
   - Turbulence equations must remain stable for large α (porous/solid regions)
   - k and ω must be bounded (e.g. k→0, ω→large in solid) to avoid stiffness

### 1.3 Design Variable and Sensitivity Chain

- **Design variable:** pseudo-density `xh` (after filtering, Heaviside)
- **Interpolation:** `α(xh) = α_max·q·(1-xh)/(q+xh)`, `DT(xh)` similarly
- **Sensitivities:** `dJ/dxh` obtained via continuous adjoint; chain rule through filter + Heaviside in `filter_chainrule.H`
- **Cost functions:** Power dissipation (pressure work), mean temperature (p-norm)

For consistent sensitivities with RANS:
- ν_t = ν_t(U, k, ω); k and ω satisfy transport equations involving U, xh (via α)
- Adjoint chain: dJ/dxh includes contributions from d(ν_t)/d(U,k,ω) · d(U,k,ω)/dxh
- This requires adjoint turbulence equations (ka, ωa) and their coupling to Ua, pa

---

## Part 2: Mathematical Formulation

### 2.1 Primal: Incompressible RANS + k-ω SST

**Momentum:**
$$\nabla\cdot(\mathbf{u}\mathbf{u}) = -\nabla p + \nabla\cdot[(\nu+\nu_t)(\nabla\mathbf{u}+\nabla\mathbf{u}^T)] - \alpha(\rho)\mathbf{u}$$

**Continuity:** ∇·u = 0  

**k-equation:**
$$\nabla\cdot(\mathbf{u} k) = \mathcal{P}_k - \beta^* \omega k + \nabla\cdot[(\nu+\sigma_k\nu_t)\nabla k]$$

**ω-equation:**
$$\nabla\cdot(\mathbf{u}\omega) = \frac{\gamma}{\nu_t}\mathcal{P}_k - \beta\omega^2 + \nabla\cdot[(\nu+\sigma_\omega\nu_t)\nabla\omega] + \mathcal{D}_\omega$$

where:
- $\mathcal{P}_k = \min(\nu_t S^2, 10\beta^*\omega k)$ (shear stress limiter)
- $S^2 = 2S_{ij}S_{ij}$, $S_{ij} = \frac{1}{2}(\partial_j u_i + \partial_i u_j)$
- $\mathcal{D}_\omega$ = cross-diffusion term (blending)
- $\nu_t = a_1 k / \max(a_1\omega, F_2 \|\nabla u\|)$

**Blending functions F1, F2**, cross-diffusion, and constants follow Menter (1994) / OpenFOAM kOmegaSST.

### 2.2 Curvature Correction (Spalart–Shur / SST-CC)

Modify production terms with factor $f_{r1}$:

$$f_{r1} = (1 + c_{r1})\frac{2r^*}{1+r^*} - c_{r1}, \quad r^* = \frac{\tilde{S}}{\tilde{\Omega}}$$

where $\tilde{S}$, $\tilde{\Omega}$ are strain and rotation invariants (Spalart–Shur formulation). Typical bounds: $f_{r1} \in [0, 1.25]$. Apply to $\mathcal{P}_k$ and $\mathcal{P}_\omega$.

**Full differentiation of $f_{r1}$** w.r.t. U (and hence strain/rotation) is required for adjoint consistency.

### 2.3 Primal Energy Equation (Turbulent)

$$(\rho c_p)\nabla\cdot(\mathbf{u} T) = \nabla\cdot[(\kappa/c_p + \alpha_t)\nabla T] + Q$$

with $\alpha_t = \nu_t / Pr_t$, $Pr_t \approx 0.85$.

### 2.4 Continuous Adjoint Formulation

**Adjoint momentum (Ua):**
$$-\nabla\cdot(\mathbf{u}\otimes\mathbf{u}_a) + (\nabla\mathbf{u})^T\mathbf{u}_a = -\nabla p_a + \nabla\cdot[(\nu+\nu_t)(\nabla\mathbf{u}_a + \nabla\mathbf{u}_a^T)] - \alpha\mathbf{u}_a + \mathbf{S}_{ka} + \mathbf{S}_{\omega a}$$

where $\mathbf{S}_{ka}$, $\mathbf{S}_{\omega a}$ are source terms from differentiation of ν_t and production terms w.r.t. U, appearing in the adjoint momentum via the chain rule.

**Adjoint pressure:** ∇·ua = 0 (enforced via adjoint continuity).

**Adjoint k (ka):**
Transport equation with:
- Convection: $-\nabla\cdot(\mathbf{u} k_a)$
- Diffusion: $\nabla\cdot[(\nu+\sigma_k\nu_t)\nabla k_a]$
- Source terms from $\partial\mathcal{L}/\partial k$ (Lagrangian derivatives of objective + constraints w.r.t. k)

**Adjoint ω (ωa):**
Similarly, with sources from $\partial\mathcal{L}/\partial\omega$.

**Key:** The Lagrangian $\mathcal{L} = J + \langle \lambda, \mathcal{N}(\mathbf{u},p,k,\omega) \rangle$ leads to cross-coupling: Ua receives sources from ka, ωa; ka and ωa receive sources from Ua and from the objective (e.g. via ν_t in momentum and in thermal diffusion).

**No frozen turbulence:** All derivatives $\partial\nu_t/\partial k$, $\partial\nu_t/\partial\omega$, $\partial\nu_t/\partial U$, $\partial\mathcal{P}_k/\partial U$, etc., must be included.

### 2.5 Sensitivity Contributions

For design variable xh:
$$\frac{dJ}{d\rho} = \frac{\partial J}{\partial \rho} + \mathbf{u}_a^T \frac{\partial \mathcal{R}_u}{\partial \rho} + p_a \frac{\partial \mathcal{R}_p}{\partial \rho} + k_a \frac{\partial \mathcal{R}_k}{\partial \rho} + \omega_a \frac{\partial \mathcal{R}_\omega}{\partial \rho}$$

where $\mathcal{R}_*$ are residuals. For Brinkman: $\partial\mathcal{R}_u/\partial\rho \propto \partial\alpha/\partial xh$. Turbulence residuals depend on ρ only through α (and possibly effective diffusivity in porous regions).

---

## Part 3: Architecture Design

### 3.1 Modular Structure (OpenFOAM Conventions)

```
MTOwithTurbulence/
├── source code/
│   └── MTO_ThermalFluid/
│       ├── src/
│       │   ├── MTO_ThermalFluid.C          # Main solver (minimal changes)
│       │   ├── createFields.H
│       │   ├── readTransportProperties.H   # Create turbulence + adjoint turbulence
│       │   ├── Primal_U.H                  # Uses turbulence->divDevReff, turbulence->correct()
│       │   ├── Primal_T.H                  # + alphat if turbulent
│       │   ├── Primal_kOmegaSST.H          # NEW: k, ω solve (or delegated to model)
│       │   ├── AdjointFlow_Ua.H            # Uses adjointTurbulence->divDevReff(Ua) + sources
│       │   ├── AdjointHeat_Ub.H            # Same
│       │   ├── AdjointHeat_Tb.H            # + adjoint turbulent diffusion
│       │   ├── Adjoint_kOmegaSST.H         # NEW: ka, ωa solve
│       │   ├── sensitivity.H               # + turbulence-related sensitivity terms
│       │   └── ...
```

**Option A — Custom turbulence library (recommended for full control):**
- Create `$FOAM_USER_LIBBIN/kOmegaSSTCC` (primal k-ω SST with curvature correction)
- Create `$FOAM_USER_LIBBIN/adjointKOmegaSST` (adjoint turbulence model)
- Link into MTO_ThermalFluid

**Option B — Extend OpenFOAM incompressible turbulence:**
- If OpenFOAM provides `kOmegaSST` and (in v2206+) `adjointKOmegaSST`, configure turbulenceProperties to use them
- Add curvature correction as a derived class or patch
- **Risk:** OpenFOAM adjoint may target shape optimization, not density-based topology; BCs and objective formulation may differ

**Recommendation:** Option A for full control, compatibility with older OpenFOAM (e.g. 6, 8), and topology-specific adjoint BCs. Option B can be evaluated if using OpenFOAM v2206+ and if its adjoint formulation matches.

### 3.2 Integration Points

| Primal | Adjoint | Sensitivity |
|--------|---------|-------------|
| `turbulence->divDevReff(U)` | `adjointTurbulence->divDevReff(Ua)` + source terms from ka, ωa | — |
| `turbulence->correct()` (k, ω) | `adjointTurbulence->correct()` (ka, ωa) | — |
| — | — | Add terms from $\partial\mathcal{R}_{k,\omega}/\partial xh$ |

### 3.3 Discretization Strategy

| Term | Scheme | Justification |
|------|--------|---------------|
| Convection (U, k, ω) | Bounded second-order (e.g. `limitedLinear`, `vanLeer`) | Stability, positivity for k, ω |
| Diffusion | `Gauss linear corrected` | Non-orthogonal correction |
| Pressure–velocity | SIMPLEC (current) | Adequate; PIMPLE if transience needed later |
| Turbulence convection | `bounded upwind` or `limitedLinear 1` | Strict positivity |
| Adjoint convection | `bounded upwind` (reverse flow) | Stability for adjoint |

**Positivity:** Use `max(k, k_min)`, `max(ω, ω_min)` in production/destruction to avoid singularities. Clip ν_t to reasonable bounds.

**Brinkman:** In solid-like cells (xh→0), α→∞; ensure k, ω solvers handle extreme source terms (implicit treatment of destruction terms).

### 3.4 Boundary Conditions

**Primal:**
- Walls: wall functions or low-Re (depending on y+)
- Inlet: k, ω from turbulence intensity and length scale
- Outlet: zeroGradient

**Adjoint:**
- Walls: homogeneous or adjoint wall functions (Kavvadias-style)
- Inlet/outlet: from objective (e.g. power dissipation, heat flux)

Existing adjoint outlet BCs must be extended to use `nuEff` when turbulence is on.

---

## Part 4: Stepwise Implementation Plan

### Phase 1: Foundation (No Turbulence Yet)
1. **Align adjoint viscous term with primal structure**
   - Replace `fvm::laplacian(nu,Ua)` with `divDevReff(Ua)` using `nuEff = nu` (laminar)
   - Ensure fvSchemes and assembled matrix match
2. **Verify sensitivity consistency** with finite differences in laminar mode
3. **Document current baseline** (laminar) for regression

### Phase 2: Primal k-ω SST
1. **Enable RAS turbulence** in turbulenceProperties (`kOmegaSST`)
2. **Primal:** Ensure k, ω solve and `turbulence->correct()` are called; add k, ω to `update.H` if needed
3. **Primal T:** Add `alphat` (turbulent thermal diffusivity) to energy equation
4. **Brinkman stability:** Enforce k, ω bounds in porous regions; test with intermediate densities
5. **Validation:** Channel flow, backward-facing step (compare with literature)

### Phase 3: Primal k-ω SST with Curvature Correction
1. **Implement curvature correction** (fr1) in production terms
2. **Validation:** Curved duct, rotating flow
3. **Ensure differentiability** of fr1 w.r.t. U (for later adjoint)

### Phase 4: Adjoint Turbulence (Consistent)
1. **Derive continuous adjoint** of k and ω equations (full linearization)
2. **Implement adjoint turbulence model** (ka, ωa transport + sources)
3. **Couple to adjoint momentum:** Add source terms from d(ν_t)/d(U,k,ω)
4. **Update adjoint BCs** to use `nuEff`
5. **Sensitivity:** Add turbulence residual contributions to `gsensPowerDiss`, `fsenshMeanT`

### Phase 5: Adjoint Curvature Correction
1. **Differentiate fr1** and production terms w.r.t. U
2. **Add adjoint sources** from curvature terms
3. **Validate** sensitivities (FD, complex step)

### Phase 6: Integration & Robustness
1. **End-to-end** topology optimization test (turbulent microchannel)
2. **Convergence studies** (mesh refinement, solver settings)
3. **Documentation** (equations, discretization, validation)

---

## Part 5: Validation Strategy

| Test | Purpose |
|------|---------|
| Manufactured solution | Verify primal + adjoint consistency |
| Channel flow (Re~1e4) | Primal turbulence, wall functions |
| Backward-facing step | Separation, recirculation |
| Curved duct | Curvature correction |
| Finite difference sensitivity | dJ/dxh vs adjoint |
| Topology optimization loop | Full cycle, convergence |

---

## Part 6: Architectural Risks and Mitigations

| Risk | Mitigation |
|------|------------|
| **OpenFOAM version drift** | Target OpenFOAM 6–8; avoid v2206-specific APIs unless adopted |
| **Stability in porous regions** | Bounded k, ω; implicit destruction; careful α_max scaling |
| **Curvature correction stiffness** | Limit fr1; smooth strain/rotation computation |
| **Adjoint solver convergence** | Relaxation; proper preconditioning; solve order (Ua, pa, ka, ωa) |
| **Thermal coupling** | Include αt in primal T; full adjoint differentiation |
| **MPI / scalability** | Use standard OpenFOAM decomposition; avoid global reductions in hot path |

---

## Part 7: Deliverables Checklist

- [ ] Primal k-ω SST (standard)
- [ ] Primal k-ω SST with curvature correction
- [ ] Adjoint k-ω SST (sensitivity-consistent)
- [ ] Adjoint curvature correction terms
- [ ] MTO_ThermalFluid solver integration
- [ ] Validation cases (channel, BFS, curved duct)
- [ ] Documentation (equations, discretization, stability)
- [ ] Sensitivity verification (FD comparison)

---

## Appendix A: Key References

1. Kavvadias et al. (2014), "The continuous adjoint approach to the k–ω SST turbulence model...", *Engineering Optimization*, 47(11), 1523–1542.
2. Spalart & Shur (1997), "On the sensitization of turbulence models to rotation and curvature", *Aerospace Science and Technology*.
3. Menter (1994), "Two-equation eddy-viscosity turbulence models...", *AIAA Journal*.
4. OpenFOAM v2206 adjoint k-ω SST (differentiated turbulence).
5. DAFoam: adjoint framework for OpenFOAM.

---

## Appendix B: Solver Execution Order (Proposed)

```
Per SIMPLE iteration:
  1. update.H (xh, alpha, DT, etc.)
  2. turbulence->correct()           # Primal k, ω
  3. Primal_U.H
  4. Primal_T.H
  5. AdjointHeat_Tb.H
  6. Adjoint_kOmegaSST.H             # ka, ωa (sources depend on Ua, Ub, Tb)
  7. AdjointHeat_Ub.H
  8. AdjointFlow_Ua.H
  9. costfunction.H
 10. sensitivity.H
```

Order of adjoint solves may require iteration (Ua ↔ ka, ωa) if strongly coupled; experience suggests sequential solve with under-relaxation can suffice.
