# Adjoint Options for Convergence Tuning

Reference for adjoint tuning parameters (formerly optProperties sections 8–10). These options have been removed from optProperties for the default case but are **retained here for future reference** when convergence issues arise in other problem cases.

**Status**: Options pruned from optProperties (run to opt 172 stable). Use this doc to restore or tune parameters if adjoint solvers diverge, crash, or fail to converge on different geometries, Reynolds numbers, or mesh sizes.

---

## When to Use This Reference

If you encounter:
- **Tb divergence** (large residuals, non-convergence)
- **Ua/Ub/ka/omegaa instability** (singular solvers, NaN/Inf)
- **Crashes** (FPE, e.g. at turbSourceUa or Tb RHS)
- **Frequent** `tb_preSolve_reset`, `tb_rhs_capped`, `uaAllPass_skipped`, or `uaHigherPass_skipped` in solver_diagnostics

→ Add the relevant parameters back to optProperties and tune per the guidance below.

---

## 1. Parameters Reference (optProperties)

All parameters are supported via `lookupOrDefault` in `createFields.H`. Add any of these to `optProperties` to override code defaults.

| Section | Option | Code Default | Purpose |
|---------|--------|--------------|---------|
| 8 | adjPseudoInvDtUa | 50 | Pseudo-transient 1/dt for Ua (diagonal dominance) |
| 8 | adjPseudoInvDtUb | 25 | Pseudo-transient 1/dt for Ub |
| 8 | adjPseudoInvDtTb | 100000 | Pseudo-transient 1/dt for Tb (critical for Tb stability) |
| 8 | adjPseudoInvDtTurb | 100 | Pseudo-transient 1/dt for ka/omegaa |
| 8 | tbPreSolveReset | 1e6 | Reset Tb to zero if max\|Tb\| exceeds (crash safeguard) |
| 8 | tbHardCapMin | 1e4 | Lower bound on Tb hard cap (per-cell clipping) |
| 8 | tbHardCapMax | 1e6 | Upper bound on Tb hard cap |
| 8 | tbRhsCapMax | 1e12 | Cap Tb RHS source to avoid FPE |
| 8 | adjTbGradTgradTbCapMax | 1e10 | Cap grad(T)·grad(Tb) in ka/omegaa thermal source |
| 9 | adjTurbBetaMin | 0.2 | Initial turbulence-adjoint continuation factor |
| 9 | adjTurbRampStart | 1 | Iteration to start turbulence beta ramp |
| 9 | adjTurbRampEnd | 80 | Iteration to finish turbulence-adjoint beta ramp |
| 9 | adjTurbLagEvery | 2 | Solve turbulence adjoint every N opt iterations |
| 9 | adjHeatBetaMin | 0.2 | Initial heat-adjoint continuation factor |
| 9 | adjHeatRampStart | 1 | Iteration to start heat-adjoint beta ramp |
| 9 | adjHeatRampEnd | 120 | Iteration to finish heat-adjoint beta ramp |
| 9 | adjUbPhiHbyAbCapFloor | 1e-7 | phiHbyAb cap floor (Ub suppression) |
| 10 | adjUaPassDecay | 0.35 | Ua-source damping per pass |
| 10 | adjUaPassMinGain | 0.10 | Min gain floor for later passes |
| 10 | adjUaLocalCapCoeff | 2000 | Turb source cap coeff on U²/L |
| 10 | adjUaRhsCapCoeff | 500 | Ua RHS cap coeff |
| 10 | adjUaAllPassTbMax | 1e74 | Skip Ua if max\|Tb\| exceeds (decouple when Tb bad) |
| 10 | adjUaHigherPassTbMax | 1e60 | Skip higher-pass Ua if max\|Tb\| exceeds |

---

## 2. Diagnostics to Check

When investigating convergence, inspect `solver_diagnostics_rank*.log` and `solver_status.log`:

| Diagnostic | Meaning | Tune |
|------------|---------|------|
| `tb_preSolve_reset` | Tb reset to zero before solve | tbPreSolveReset, adjPseudoInvDtTb |
| `tb_rhs_capped` | Tb RHS scaled down | tbRhsCapMax |
| `uaAllPass_skipped` | Ua solve skipped (Tb too large) | adjUaAllPassTbMax |
| `uaHigherPass_skipped` | Higher-pass Ua skipped | adjUaHigherPassTbMax |
| Tb maxFinRes high / conv=0 | Tb not converging | adjPseudoInvDtTb (raise to 1e5–2e5) |
| Ua/Ub/ka/omegaa sing=1 | Singular solver | adjPseudoInvDt* |

---

## 3. Crash/Stability Safeguards (ESSENTIAL when unstable)

| Option | Suggested value | Role |
|--------|-----------------|------|
| **adjPseudoInvDtTb** | 100000–200000 | Strong diagonal for Tb; CRASH_INVESTIGATION_opt313 |
| **tbPreSolveReset** | 1e5–1e6 | Reset Tb when it blows up |
| **tbRhsCapMax** | 1e12 | Cap Tb RHS to avoid FPE |
| **adjTbGradTgradTbCapMax** | 1e10 | Cap ka/omegaa thermal source |
| **adjUaAllPassTbMax** | 1e6 | Skip Ua when Tb is bad |

---

## 4. Example: Restoring Tuning for a Difficult Case

```cpp
// Add to optProperties when adjoint convergence is poor:

adjPseudoInvDtTb    150000;   // Tb not converging
tbPreSolveReset     1e5;      // Tb blowing up
tbRhsCapMax         1e11;     // Tb RHS extreme
adjUaAllPassTbMax   1e6;      // Ua decouple when Tb bad
```

---

## 5. Log Analysis (baseline run, opt 172)

| Diagnostic | Count | Result |
|------------|-------|--------|
| tb_preSolve_reset | 0 | Never triggered |
| tb_rhs_capped | 0 | Never triggered |
| uaAllPass_skipped = 1 | 0 | Ua never skipped |
| uaHigherPass_skipped = 1 | 0 | Higher-pass Ua never skipped |

Baseline (test_area=yes, filterR=2.5, U=1 m/s, Dh=8 mm) completed 172 iterations with no safeguard triggers. Code defaults are sufficient for that case.

---

*Reference document for convergence tuning. Parameters removed from optProperties but available via lookupOrDefault when added to optProperties.*
