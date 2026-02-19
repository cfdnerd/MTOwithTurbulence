# MTO_ThermalFluid Turbulence Implementation Report

**Date:** February 2025  
**Scope:** Phases 1–6 implementation

---

## Summary of Implemented Changes

### Phase 1: Foundation — Adjoint Viscous Operator Consistency ✅

**Changes:**
- Replaced `fvm::laplacian(nu, Ua)` with `fvm::laplacian(turbulence->nuEff(), Ua)` in `AdjointFlow_Ua.H`
- Replaced `fvm::laplacian(nu, Ub)` with `fvm::laplacian(turbulence->nuEff(), Ub)` in `AdjointHeat_Ub.H`
- Added `volScalarField nuEff` in `readTransportProperties.H`, updated each iteration in `update.H`
- Updated all four adjoint outlet BCs to use `nuEff` from objectRegistry instead of laminar `nu`
- Added `div((nuEff*dev2(T(grad(Ub)))))` to fvSchemes

**Result:** Primal and adjoint viscous operators are structurally aligned. When laminar, ν_eff = ν; when turbulent, ν_eff = ν + ν_t.

---

### Phase 2: Primal k-ω SST Integration ✅

**Changes:**
- Set `simulationType RAS` and `RASModel kOmegaSST` in turbulenceProperties
- Created `0/k` and `0/omega` with appropriate BCs (kqRWallFunction, omegaWallFunction)
- Added `turbulence->correct()` before Primal_U in main loop
- Added `alphaEff = DT + turbulence->alphat()` in Primal_T.H for turbulent thermal diffusion
- Updated AdjointHeat_Tb.H to use `alphaEff`
- Added omega to fvSolution solvers and relaxation

**Result:** Primal flow uses k-ω SST; energy equation includes turbulent thermal diffusivity.

---

### Phase 3: Curvature Correction — Deferred ⚠️

**Status:** Not implemented. Curvature correction requires a custom `kOmegaSSTCC` library extending OpenFOAM's kOmegaSST. The turbulence model hierarchy differs between OpenFOAM versions (5/6 vs 8+), making a portable library non-trivial.

**Recommendation:** Use ancolli/kOmegaSSTCC or c-schubert/openfoam_kOmegaSST_curvate_corr as reference. Add `kOmegaSSTCC` to run-time selection when the library is built. Placeholder comments added in turbulenceProperties.

---

### Phase 4: Adjoint Turbulence (ka, ωa) ✅

**Changes:**
- Added `ka` and `omegaa` fields in `readTransportProperties.H`
- Created `Adjoint_kOmegaSST.H` with adjoint k and omega transport equations
- Adjoint equations solved only when `k` and `omega` exist (RAS active)
- Added `0/ka` and `0/omegaa` initial conditions
- Inserted `#include "Adjoint_kOmegaSST.H"` in main loop (before AdjointHeat_Ub, AdjointFlow_Ua)
- Added ka, omegaa to fvSolution

**Known limitations:**
1. **Source-term coupling:** Adjoint momentum (Ua) does not yet receive source terms from `d(ν_t)/dk * ka` and `d(ν_t)/dω * ωa`. Full Kavvadias formulation requires these for strict sensitivity consistency. Current implementation uses ν_eff in the viscous operator but omits the cross-coupling.
2. **Adjoint BCs:** Adjoint k and ω use zeroGradient at walls; adjoint wall functions (Kavvadias-style) are not implemented.
3. **Solve order:** ka and ωa are solved sequentially; iterative coupling may improve convergence for stiff cases.

---

### Phase 5: Adjoint Curvature Correction — Cancelled

**Status:** Not applicable until Phase 3 (primal curvature correction) is implemented.

---

### Phase 6: Integration, Validation, Documentation

**Files modified/created:**
- Architecture plan: `docs/MTO_ThermalFluid_Turbulence_Architecture_Plan.md`
- Equations reference: `docs/Adjoint_SST_Equations_Reference.md`
- Implementation report: `docs/Implementation_Report.md` (this file)

**Validation checklist (user action required):**
- [ ] Build MTO_ThermalFluid against OpenFOAM (wmake)
- [ ] Run with laminar (simulationType laminar) — baseline
- [ ] Run with RAS kOmegaSST — verify k, omega, nut, nuEff
- [ ] Compare sensitivities via finite differences
- [ ] Topology optimization loop with turbulent case

---

## Inconsistencies and Drawbacks for Further Correction

| Issue | Severity | Action |
|-------|----------|--------|
| **Adjoint momentum–turbulence coupling** | High | Add source terms from `∂ν_t/∂k` and `∂ν_t/∂ω` into Ua equation (Kavvadias Sec. 3) |
| **Curvature correction missing** | Medium | Build kOmegaSSTCC library; add to turbulenceProperties when available |
| **Adjoint wall BCs** | Medium | Implement adjoint wall functions for ka, ωa at walls |
| **Sensitivity to alphat** | Low | If objective depends on T, include ∂(αt)/∂k, ∂(αt)/∂ω in adjoint thermal chain |
| **fvSchemes for k/ω** | Low | Consider `bounded Gauss upwind` for k, ω if positivity issues arise |
| **Brinkman–turbulence stability** | Medium | Monitor k, ω in porous regions; enforce bounds in update if needed |

---

## Laminar Fallback

To run laminar (e.g., for validation): set `simulationType laminar;` and `RASModel laminar;` in turbulenceProperties. Remove or avoid `0/k` and `0/omega` in the case, or ensure the laminar model does not require them.

---

## References

1. Kavvadias et al. (2014), Engineering Optimization, 47(11), 1523–1542
2. OpenFOAM v2206 adjoint k-ω SST
3. MTO_ThermalFluid_Turbulence_Architecture_Plan.md
4. Adjoint_SST_Equations_Reference.md
