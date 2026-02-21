# Crash Investigation: PowerDiss0 Trial at Iteration 313

## Executive Summary

The solver crashed during **optimization iteration 313** in the **Ua (adjoint momentum) solve**, in the `UaLoop_0` stage immediately **after** `uaTrace_beforeCapLoop`. With `FOAM_SIGFPE` enabled, the crash is consistent with a floating-point exception. The failure chain is: **Tb non-convergence** → Tb clipped to ±1000 → very large `grad(Tb)` at clip boundaries → **ka/omegaa thermal RHS** (for objFunction=1) blows up → ka/omegaa become extreme → **turbSourceUa** overflows or yields NaN/Inf → crash in sanitization or in `gMax(mag(turbSourceUa))`.

**Critical correction**: This run uses `objFunction=1` (mean temperature). **Pnorm is not used**; it only applies when `objFunction=2` (p-norm max temperature). The Tb equation RHS for objFunction=1 is the constant Q from thermalProperties.

---

## 1. Where the Crash Occurred

### 1.1 Last Successful Operation

All 60 rank diagnostic logs show the same last entry:

```
313 UaLoop_0 uaTrace_beforeCapLoop 1.000000e+00
```

The crash happens **after** this log and **before** the next diagnostic (`turbSourceUa_magMin`, etc.).

### 1.2 Failing Code Location

In `AdjointFlow_Ua.H`, the execution order is:

1. `MTOdiag::logMetric(..., "uaTrace_beforeCapLoop", 1.0);` ✓ (completed)
2. Sanitization loop over `turbSourceUa`:
   ```cpp
   forAll(turbSourceUa, cellI) {
       vector& sv = turbSourceUa[cellI];
       if (!(std::isfinite(sv.x()) && std::isfinite(sv.y()) && std::isfinite(sv.z())))
           sv = vector::zero;
       else {
           sv.x() = Foam::max(-compCap, Foam::min(compCap, sv.x()));
           sv.y() = Foam::max(-compCap, Foam::min(compCap, sv.y()));
           sv.z() = Foam::max(-compCap, Foam::min(compCap, sv.z()));
       }
   }
   ```
3. Global scaling: `maxMag = gMax(mag(turbSourceUa.primitiveField()));`

Most likely failure point: use of `Foam::min`/`Foam::max` on **signaling NaNs** in `turbSourceUa`, which can trigger an FPE on architectures that treat NaN comparisons as invalid operations.

---

## 2. Root Cause: Adjoint Thermal Field Instability

### 2.1 Tb (Adjoint Temperature) Divergence

`solver_status.log` shows the Tb final residual at iteration 313 and prior iterations:

| optIter | Tb maxFinRes   | conv | Notes                          |
|---------|----------------|------|--------------------------------|
| 310     | 9.04e+08       | 0    | Not converged                  |
| 311     | 2.54e+01       | 0    | Still high                     |
| 312     | 1.69e+08       | 0    | Not converged                  |
| 313     | (crash before write) | - | Never written (crash in Ua)    |

- Tb consistently **does not converge** (max iterations 1500 hit with very large residuals).
- Tb is **clipped** to ±1000 (`tb_hardCap`, `tb_absMax`).
- At 313, the Tb RHS was capped (`tb_rhs_capped` logged), indicating extreme source terms.

### 2.2 Propagation Path (objFunction=1: Mean Temperature)

1. **Tb instability** → Tb does not converge; Tb is clipped to ±1000. `grad(Tb)` becomes very large at clip boundaries.
2. **Tb RHS capping** (`tb_rhs_capped`): the assembled Tb source exceeds `tbRhsCapMax` (likely from the div term using the clipped Tb).
3. **ka/omegaa thermal source** (Adjoint_kOmegaSST.H, `objFunction == 1` only):
   ```cpp
   kaRhsThermal   = (k+kSmall)²/Tref² * (dNutdk/Prt) * (grad(T) & grad(Tb));
   omegaaRhsThermal = ... * (grad(T) & grad(Tb));
   ```
   When Tb is clipped, `grad(Tb)` can be huge → `gradTgradTb` explodes → ka/omegaa RHS becomes extreme.
4. **turbSourceUa** in `AdjointFlow_Ua.H`:
   ```cpp
   turbSourceUa = -fvc::div(ka*2*nutSafe*twoSymm(fvc::grad(U)))/kDenom
                  - nutSafe*fvc::div(omegaa*2*gamma1*twoSymm(fvc::grad(U)))/omegaDenom;
   ```
   If ka or omegaa contain NaN or extreme values, this expression can overflow or produce NaN/Inf.
5. **Crash**: sanitization loop or `gMax(mag(turbSourceUa))` → FPE (e.g. overflow in mag, or signaling NaN in min/max).

### 2.3 Special Condition at 313: `doTurbAdjSolve = 0`

- `optProperties`: `adjTurbLagEvery = 2`, `adjTurbRampEnd = 80`.
- At opt 313: `doTurbAdjSolveThisOpt = false` (313 > 80 and 313 % 2 ≠ 0).
- However, pass 0 still runs the turbulence adjoint (`doTurbAdjSolveThisOpt || pass == 0`), so ka/omegaa **are** solved at 313.
- They are driven by Ua from opt 312 and Tb from opt 313 (with RHS capping).

---

## 3. Supporting Evidence

### 3.1 runLOG

- Ends mid-line during `"time step continuity errors :"` (primal U/p output).
- `sigFpe : Enabling floating point exception trapping (FOAM_SIGFPE)` is enabled.
- Typical of an FPE abort during a later stage (adjoint) while primal output is still in the buffer.

### 3.2 solver_status.log

- Last complete iteration: **312**.
- No entries for 313 (crash before `writeSolverStatus.H`).

### 3.3 Design State and Run Configuration

- `test_area = yes` → xh fixed at 1 in the test zone (all fluid).
- Rank0: `xh_min = xh_max = 1`, `alpha = 0`, `alpha_nBad = 0`.
- **objFunction = 1** (mean temperature). Pnorm is not used; Tb RHS is the constant Q from thermalProperties.
- Design is not the direct cause; the problem is in the adjoint chain.

---

## 4. Recommended Fixes

### 4.1 Immediate: Guard NaN/Inf Before Sanitization

Before the sanitization loop, ensure `turbSourceUa` is checked and zeroed if non-finite:

```cpp
// Before the forAll sanitization loop, add:
forAll(turbSourceUa, cellI) {
    vector& sv = turbSourceUa[cellI];
    if (!(std::isfinite(sv.x()) && std::isfinite(sv.y()) && std::isfinite(sv.z()))) {
        sv = vector::zero;
    }
}
// Then do the clamping in a separate loop, or combine but avoid min/max on potentially bad values
```

Alternatively, use a compare-only path for non-finite values to avoid FPE in min/max:

```cpp
forAll(turbSourceUa, cellI) {
    vector& sv = turbSourceUa[cellI];
    const bool bad = !(std::isfinite(sv.x()) && std::isfinite(sv.y()) && std::isfinite(sv.z()));
    if (bad) {
        sv = vector::zero;
    } else {
        sv.x() = Foam::max(-compCap, Foam::min(compCap, sv.x()));
        sv.y() = Foam::max(-compCap, Foam::min(compCap, sv.y()));
        sv.z() = Foam::max(-compCap, Foam::min(compCap, sv.z()));
    }
}
```

(Current code already branches on `isfinite`; the important part is that no min/max is applied to non-finite components.)

### 4.2 Short-term: Stabilize Tb (objFunction=1)

**Note**: Pnorm applies only when objFunction=2; it is irrelevant for objFunction=1.

- Increase `adjPseudoInvDtTb` (e.g. 1e5 → 2e5) to improve diagonal dominance in the adjoint convection.
- Use `tbRhsCapMax` more aggressively (e.g. lower it if RHS is still extreme after capping).
- Tune `tbPreSolveReset`: reset Tb to zero when |Tb| exceeds this limit before the solve.
- Consider relaxing Tb or resetting it when `tb_finalRes` exceeds a trust threshold (e.g. 1e4) instead of clipping.

### 4.3 Medium-term: Decouple When Tb Is Bad

- When `tb_finalRes > threshold` (e.g. 1e6), skip or heavily damp the ka/omegaa thermal source (objFunction=1 only) in Adjoint_kOmegaSST.H.
- Apply `tbPreSolveReset` earlier when Tb grows beyond a trust region.

### 4.4 Optional: Disable FPE During Critical Sections

- For forensics or robustness, consider temporarily disabling `FOAM_SIGFPE` around the turbSourceUa cap loop.
- Prefer fixing NaNs/Infs at the source rather than relying on this.

---

## 5. Conclusion

| Item | Finding |
|------|---------|
| **Crash point** | Ua solve, UaLoop_0, right after `uaTrace_beforeCapLoop` |
| **Likely trigger** | FPE from NaN/Inf or overflow in `turbSourceUa` during sanitization or `gMax(mag(...))` |
| **Upstream cause** | Unstable Tb → clipping → large grad(Tb) → ka/omegaa thermal RHS (objFunction=1) → ka/omegaa extreme → turbSourceUa overflow/NaN |
| **objFunction** | 1 (mean temperature); Pnorm is irrelevant |
| **test_area** | Not causal; xh=1, alpha=0 in test zone |

Addressing Tb stability (adjPseudoInvDtTb, tbRhsCapMax, tbPreSolveReset) and decoupling ka/omegaa thermal source when Tb is bad, plus robust NaN/Inf handling in the turbSourceUa path, should prevent this crash.

---

## 6. Patches Applied (Post-Investigation)

The following patches were implemented to avoid this class of crash:

| Location | Change |
|----------|--------|
| **AdjointFlow_Ua.H** | Two-pass sanitization with compCap=1e50; guard gMax(mag(...)) and prevTurbSourceUaMax against non-finite values |
| **Adjoint_kOmegaSST.H** | Cap gradTgradTb per cell (adjTbGradTgradTbCapMax) and zero non-finite before use in ka/omegaa thermal source |
| **AdjointHeat_Tb.H** | Zero non-finite Tb RHS entries before gMax; guard srcMax in RHS cap block |
| **createFields.H** | Add adjTbGradTgradTbCapMax parameter |
| **optProperties** | Add adjTbGradTgradTbCapMax 1e10 |
