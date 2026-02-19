# k-ω SST Adjoint Equations — Implementation Reference

This document provides the detailed adjoint equation structure and source-term expressions needed for implementation. It complements `MTO_ThermalFluid_Turbulence_Architecture_Plan.md`.

---

## 1. Primal k-ω SST (Summary)

**Eddy viscosity:**
$$\nu_t = \frac{a_1 k}{\max(a_1\omega, \; F_2 \cdot \|\nabla\times\mathbf{u}\|)}$$

**Production (with curvature correction):**
$$\mathcal{P}_k = f_{r1} \cdot \min\left(\nu_t S^2,\; 10\beta^* \omega k\right)$$

**k-equation:**
$$\nabla\cdot(\mathbf{u} k) - \nabla\cdot[(\nu + \sigma_{k}\nu_t)\nabla k] = \mathcal{P}_k - \beta^* \omega k$$

**ω-equation:**
$$\nabla\cdot(\mathbf{u}\omega) - \nabla\cdot[(\nu + \sigma_{\omega}\nu_t)\nabla\omega] = \frac{\gamma}{\nu_t}\mathcal{P}_k - \beta\omega^2 + \mathcal{D}_\omega$$

**Momentum viscous term:** $\nabla\cdot[(\nu+\nu_t)\mathrm{dev}(\nabla\mathbf{u} + \nabla\mathbf{u}^T)]$

---

## 2. Adjoint Momentum Source Terms from Turbulence

The adjoint momentum equation receives extra source terms because $\nu_t = \nu_t(k,\omega,U)$ and appears in:
- Viscous flux: $\nabla\cdot[(\nu+\nu_t)\mathrm{dev}(\nabla\mathbf{u}_a + \nabla\mathbf{u}_a^T)]$
- Additional terms from $\partial\nu_t/\partial U$ and $\partial\nu_t/\partial k$, $\partial\nu_t/\partial\omega$ when varying U, k, ω in the Lagrangian.

Define $\mathcal{R}_U$ = momentum residual. The adjoint momentum right-hand side includes:
$$\mathbf{S}_{\text{turb}} = -\frac{\partial\mathcal{R}_U}{\partial k}\bigg|^T k_a - \frac{\partial\mathcal{R}_U}{\partial\omega}\bigg|^T \omega_a + \text{(terms from objective)}$$

For $\mathcal{R}_U \ni \nabla\cdot[(\nu+\nu_t)\mathrm{dev}(\nabla\mathbf{u})]$:
- Variation in k: $\nabla\cdot\left(\frac{\partial\nu_t}{\partial k}\mathrm{dev}(\nabla\mathbf{u})\right) k_a$ contributes a volume source in the Ua equation after integration by parts. The adjoint operator acts on Ua; the k-variation gives a source term proportional to $\frac{\partial\nu_t}{\partial k}\mathrm{dev}(\nabla\mathbf{u}) k_a$.
- Similarly for ω.

**Implementationally:** The adjoint momentum equation is:
$$-\mathbf{u}\cdot\nabla\mathbf{u}_a + (\nabla\mathbf{u})^T\mathbf{u}_a - \nabla p_a + \nabla\cdot[(\nu+\nu_t)(\nabla\mathbf{u}_a + \nabla\mathbf{u}_a^T)] - \alpha\mathbf{u}_a = \mathbf{S}_{\text{turb}}$$

where $\mathbf{S}_{\text{turb}}$ aggregates all sources from turbulence differentiation (see Kavvadias et al. for full expressions).

---

## 3. Adjoint k-Equation (ka)

Form: $-\nabla\cdot(\mathbf{u} k_a) - \nabla\cdot[(\nu+\sigma_k\nu_t)\nabla k_a] + \text{sources} = \text{RHS}$

**Source terms:**
- From $\partial\mathcal{P}_k/\partial k$: contributes to diagonal / explicit term in ka equation
- From $\partial\nu_t/\partial k$ in momentum: when we take $\langle \mathbf{u}_a, \delta_k\mathcal{R}_U \rangle$, we get terms like $\mathbf{u}_a \cdot \nabla\cdot\left(\frac{\partial\nu_t}{\partial k}\mathrm{dev}(\nabla\mathbf{u})\right)$ → after IBP, contributes to ka equation
- From $\partial(\beta^*\omega k)/\partial k = \beta^*\omega$
- From objective J if J depends on k (e.g. via ν_t in dissipation)

---

## 4. Adjoint ω-Equation (ωa)

Form: $-\nabla\cdot(\mathbf{u} \omega_a) - \nabla\cdot[(\nu+\sigma_\omega\nu_t)\nabla\omega_a] + \text{sources} = \text{RHS}$

**Source terms:**
- From $\partial\nu_t/\partial\omega$ in momentum
- From $\partial(\beta\omega^2)/\partial\omega = 2\beta\omega$
- From $\partial(\gamma\mathcal{P}_k/\nu_t)/\partial\omega$
- From cross-diffusion $\mathcal{D}_\omega$ differentiation
- From curvature correction $f_{r1}$ if it depends on ω (usually via strain/rotation, not ω directly)

---

## 5. Curvature Correction Adjoint

$f_{r1}$ depends on $\tilde{S}$, $\tilde{\Omega}$ (functions of velocity gradient). So:
$$\frac{\partial f_{r1}}{\partial U} = \frac{\partial f_{r1}}{\partial \tilde{S}}\frac{\partial\tilde{S}}{\partial(\nabla U)} + \frac{\partial f_{r1}}{\partial \tilde{\Omega}}\frac{\partial\tilde{\Omega}}{\partial(\nabla U)}$$

This contributes to:
1. Adjoint momentum (via $\partial\mathcal{P}_k/\partial U$ and $\partial\mathcal{P}_\omega/\partial U$)
2. Adjoint k and ω (via $\partial\mathcal{P}_k/\partial f_{r1}$, $\partial\mathcal{P}_\omega/\partial f_{r1}$)

**Implementation:** Compute $\partial f_{r1}/\partial r^*$ analytically, then $r^* = \tilde{S}/\tilde{\Omega}$, and use chain rule for $\partial r^*/\partial(\nabla U)$.

---

## 6. Objective-Dependent Terms

**Power dissipation:** $J = -\int_{\Gamma_o} \phi\left(p + \frac{1}{2}|\mathbf{u}|^2\right) dS$

- Direct dependence on U, p → standard boundary terms in Ua, pa BCs
- No direct k, ω dependence → objective does not add RHS to ka, ωa from J itself
- Indirect: ν_t affects U, p; so adjoint sources come through PDE residuals

**Mean temperature (p-norm):** $J \propto \left(\int T^p\right)^{1/p}$

- T depends on αt = νt/Prt → J depends on νt → ∂J/∂k, ∂J/∂ω appear as sources in ka, ωa

---

## 7. Discretization Notes

- **Convection:** Adjoint has reverse flow; use `bounded upwind` or equivalent for ka, ωa
- **Diffusion:** Same structure as primal (symmetric)
- **Sources:** Assemble explicitly; use `fvm::Sp` for linearized diagonal contributions where appropriate
- **Solve order:** Typically: Ua, pa → ka → ωa (or iterative coupling)

---

*For complete derivation, consult Kavvadias et al. (2014) and OpenFOAM adjoint source code (v2206+).*
