# Strategic Numerical Revamp — Adjoint Turbulence Stability

**Context**: MTO_ThermalFluid adjoint solver FPE at iteration ~148 in the Ua phase, specifically during the turbSourceUa cap loop (between `uaTrace_beforeCapLoop` and `turbSourceUa_magMin`). Parameter tuning and local guards have not resolved the issue. This document outlines strategic alternatives based on literature and best practices.

---

## 1. Root Cause Summary

- **Where**: Cap loop over `turbSourceUa` cells (forAll)
- **Symptom**: FPE (signal 8) on ranks 3, 6, etc., at the same code offset
- **Likely cause**: Operations on partition-specific data (e.g. `mag(sv)`, `mag(U[cellI])`) with extreme values that overflow or trigger FPE on some ranks
- **Current mitigations**: Diagonal guards, nutSafe, k/omega denominator floors, nuEff floor, outlet BC magSf guards, simplified global cap (no per-cell Um/hChar)

---

## 2. Literature-Based Strategies

### A. Implicit vs Explicit Source Treatment

**Observation** (OpenFOAM/SU2): Explicit adjoint turbulence sources add stiffness and can destabilize. Implicit contributions (e.g. `fvm::Sp`) improve conditioning.

**Current**: `turbSourceUa` is fully explicit on the RHS:  
`uaRHS = -fvc::grad(pa) + turbSourceUa`  
`solve(UaEqn == uaRHS)`

**Options**:
- **Diagonal stabilization**: Add `fvm::Sp(sigmaUa, Ua)` with `sigmaUa` chosen to dominate unstable contributions. This damps high-frequency modes.
- **Partial implicit linearization**: If turbSourceUa could be partially linearized in Ua (e.g. via ka/omegaa dependence), move that part into the matrix.

**Limitation**: turbSourceUa depends on ka, omegaa, U, grad(U), nut, k, omega—not on Ua. So a true implicit linearization in Ua is not straightforward. The diagonal `fvm::Sp` stabilization remains viable.

---

### B. Frozen Turbulence Fallback (Graceful Degradation)

**Observation** (OpenFOAM v2206, Kavvadias): Frozen turbulence (nut constant) is less accurate but more stable. Full differentiation can diverge; frozen gives suboptimal but bounded behavior.

**Strategy**: At each Ua pass, decide:
- If trust gates pass (Tb bounded, etc.) → use full turbSourceUa.
- Otherwise → use `turbSourceUa = 0` (or heavily damped).

**Implementation**: Extend the existing `adjUaAllPassTbMax` logic. Add a second gate (e.g. `adjUaTurbSourceMax` or residual-based) and, when tripped, set `turbSourceUa = 0` and optionally reduce `adjTurbBetaMin` for that iteration.

---

### C. Source Term Relaxation / Blending

**Observation** (homotopy, continuation): Sudden changes in source terms cause instability. Blending with previous values improves robustness.

**Strategy**:  
`turbSourceUa = (1 - omegaSource)*turbSourceUa_prev + omegaSource*turbSourceUa_new`  
with `omegaSource` in (0, 1).

**Implementation**: Store `turbSourceUa` from the previous UaLoop or opt iteration and blend. Start with `omegaSource = 0.3–0.5` and tune.

---

### D. Spatial Smoothing / Filtering of turbSourceUa

**Observation** (artificial viscosity, gradient smoothing): High-frequency components in the turbulence source drive instability; smoothing damps them.

**Strategy**: Apply a Laplacian or explicit filter to turbSourceUa before adding to the RHS:
```
turbSourceUa_smooth = fvc::laplacian(alphaSmooth*mesh.delta()^2, turbSourceUa)
turbSourceUa = (1 - betaSmooth)*turbSourceUa + betaSmooth*turbSourceUa_smooth
```
or a simpler relaxation: `turbSourceUa = turbSourceUa + alphaSmooth*fvc::laplacian(turbSourceUa)`.

**Implementation**: Add a diffusion term for turbSourceUa with a small diffusivity (e.g. `alphaSmooth * sqr(h)`).

---

### E. Two-Phase Ua Solve (Predictor–Corrector)

**Observation** (IMEX, operator splitting): Splitting stiff and non-stiff parts improves stability.

**Strategy**:
1. **Phase 1**: Solve Ua without turbSourceUa (frozen turbulence):  
   `solve(UaEqn == -fvc::grad(pa))`
2. **Phase 2**: Compute turbSourceUa from Phase 1 fields, then solve again with under-relaxation:  
   `Ua = (1 - relax)*Ua + relax*solve(UaEqn == uaRHS_full)`

**Effect**: First phase yields a stable base; second phase adds the turbulence source gradually.

---

### F. Safe Magnitude and Overflow Guards

**Observation**: FPE in the cap loop likely comes from `mag(sv)` or `mag(U[cellI])` when components are very large, causing overflow in `x² + y² + z²` or `sqrt`.

**Strategy**: Use a robust magnitude:
```cpp
inline scalar magSafe(const vector& v)
{
    const scalar x = v.x(), y = v.y(), z = v.z();
    const scalar m2 = x*x + y*y + z*z;
    if (!std::isfinite(m2) || m2 > 1e300) return 1e150;  // cap to avoid overflow
    return Foam::sqrt(Foam::max(m2, 0.0));
}
```

Replace `mag(sv)` and `mag(U[cellI])` in the cap loop with `magSafe(...)`.

---

### G. Bypass Cap Loop Entirely (Emergency Fallback)

**Strategy**: If the cap loop continues to FPE, avoid per-cell operations:
- Cap turbSourceUa only via a global operation (e.g. `minMod` or scaling of the whole field).
- Use `turbSourceUa = min(mag(turbSourceUa), srcMaxAbs) * sign(turbSourceUa)` conceptually, or scale the field so `max(mag(turbSourceUa)) <= srcMaxAbs` without a cell-wise loop.

**Implementation**: Replace the forAll cap loop with:
```cpp
scalar maxMag = gMax(mag(turbSourceUa.primitiveField()));
if (maxMag > capVal)
    turbSourceUa *= (capVal / Foam::max(maxMag, VSMALL));
```

This uses collective `gMax` and avoids partition-specific per-cell logic.

---

### H. Solve Order and Coupling

**Observation** (SU2, discrete adjoint): Block-coupled or more tightly coupled schemes are more stable than fully segregated ones.

**Strategy**:
- Use `simple.consistent()` (SIMPLE-C) for Ua–pa if not already.
- Consider solving Ua and pa in a block-coupled manner for the adjoint.
- Optionally solve Ua before Ub/Tb in some passes to reduce feedback from thermal adjoints.

---

### I. Homotopy on Turbulence Source

**Observation** (monolithic homotopy, continuation): Gradually turning on the turbulence source improves robustness.

**Strategy**: Ramp turbSourceUa from zero to full over optimization iterations:
- `turbSourceUa *= homotopy(1.0 - xi, xi)` where `xi` goes from 0 to 1 over `adjTurbRampStart`–`adjTurbRampEnd`.

**Current**: `adjTurbBetaMin` and `adjTurbRamp*` already implement a ramp. Consider an additional homotopy specifically on the magnitude of turbSourceUa (e.g. cap by `xi * srcMaxAbs`).

---

### J. Discrete Adjoint (Long-Term)

**Observation**: Discrete adjoint (AD) is more consistent and robust than continuous adjoint for complex turbulence models.

**Strategy**: Derive the adjoint via automatic differentiation of the discretized primal. This is a major refactor but would align with SU2-style practices.

---

## 3. Recommended Implementation Order

| Priority | Strategy | Effort | Impact |
|----------|----------|--------|--------|
| 1 | **G. Bypass cap loop** (use gMax + global scale) | Low | Directly avoids FPE in cap loop |
| 2 | **F. Safe mag()** in cap loop and related code | Low | Removes overflow/FPE in mag |
| 3 | **B. Frozen turbulence fallback** when gates fail | Medium | Graceful degradation |
| 4 | **C. Source blending** with previous turbSourceUa | Medium | Smoother source history |
| 5 | **D. Spatial smoothing** of turbSourceUa | Medium | Dampens high-frequency modes |
| 6 | **A. Diagonal fvm::Sp** for extra stabilization | Low | Better conditioning |
| 7 | **E. Two-phase Ua solve** | Higher | Structural change |
| 8 | **J. Discrete adjoint** | Very high | Long-term architecture |

---

## 4. References

- Kavvadias et al., “The continuous adjoint approach to the k–ω SST turbulence model,” *Engineering Optimization* (2015)
- OpenFOAM v2206 numerics: differentiated vs frozen turbulence
- SU2: duality-preserving discrete adjoint, AIAA-2016-3518
- Wang & Akbarzadeh, “Stabilisation of discrete adjoint solvers for incompressible flow”
- Artificial viscosity / IMEX for adjoint equations (MIT, Michigan)
- Brinkman penalization in topology optimization (arXiv:2302.14156)
