# Adjoint Options Evaluation — optProperties Sections 8–10

Evaluation of adjoint tuning options to identify essential vs removable parameters.  
**Status**: Preliminary (awaiting solver_diagnostics logs from runs).

---

## 1. Options Under Scrutiny

| Section | Option | Current Default (optProperties) | Code Default (lookupOrDefault) | Purpose |
|---------|--------|--------------------------------|--------------------------------|---------|
| 8 | adjPseudoInvDtUa | 50 | 50 | Pseudo-transient 1/dt for Ua (diagonal dominance) |
| 8 | adjPseudoInvDtUb | 40 | 25 | Pseudo-transient 1/dt for Ub |
| 8 | adjPseudoInvDtTb | 100000 | 10 | Pseudo-transient 1/dt for Tb |
| 8 | adjPseudoInvDtTurb | 100 | 100 | Pseudo-transient 1/dt for ka/omegaa |
| 8 | tbPreSolveReset | 1e5 | 1e6 | Reset Tb to zero if max\|Tb\| exceeds this (crash safeguard) |
| 8 | tbHardCapMin | 1e5 | 1e4 | Lower bound on Tb hard cap (per-cell clipping) |
| 8 | tbHardCapMax | 1e6 | 1e6 | Upper bound on Tb hard cap |
| 8 | tbRhsCapMax | 1e12 | 1e12 | Cap Tb RHS source to avoid FPE |
| 8 | adjTbGradTgradTbCapMax | 1e10 | 1e10 | Cap grad(T)·grad(Tb) in ka/omegaa thermal source |
| 9 | adjTurbBetaMin | 0.2 | 0.2 | Initial turbulence-adjoint continuation factor |
| 9 | adjTurbRampStart | 1 | 1 | Iteration to start turbulence beta ramp |
| 9 | adjTurbRampEnd | 80 | 80 | Iteration to finish turbulence beta ramp |
| 9 | adjTurbLagEvery | 2 | 2 | Solve turbulence adjoint every N opt iterations |
| 9 | adjHeatBetaMin | 0.1 | 0.2 | Initial heat-adjoint continuation factor |
| 9 | adjHeatRampStart | 1 | 1 | Iteration to start heat beta ramp |
| 9 | adjHeatRampEnd | 220 | 120 | Iteration to finish heat beta ramp |
| 9 | adjUbPhiHbyAbCapFloor | 1e-7 | 1e-7 | phiHbyAb cap floor (Ub suppression) |
| 10 | adjUaPassDecay | 0.35 | 0.35 | Ua-source damping per pass |
| 10 | adjUaPassMinGain | 0.10 | 0.10 | Min gain floor for later passes |
| 10 | adjUaLocalCapCoeff | 2000 | 2000 | Turb source cap coeff on U²/L |
| 10 | adjUaRhsCapCoeff | 500 | 500 | Ua RHS cap coeff |
| 10 | adjUaAllPassTbMax | 1e6 | 1e74 | Skip Ua if max\|Tb\| exceeds (decouple when Tb bad) |
| 10 | adjUaHigherPassTbMax | 1e5 | 1e60 | Skip higher-pass Ua if max\|Tb\| exceeds |

---

## 2. Diagnostics to Check (solver_diagnostics / solver_status)

When you provide `solver_diagnostics_rank*.log` and `solver_status.log`, look for:

| Diagnostic | Meaning | Informs |
|------------|---------|---------|
| `tb_preSolve_reset` | Tb reset to zero before solve | tbPreSolveReset necessity |
| `tb_rhs_capped` | Tb RHS scaled down | tbRhsCapMax necessity |
| `uaAllPass_skipped` | Ua solve skipped (Tb too large) | adjUaAllPassTbMax necessity |
| `uaHigherPass_skipped` | Higher-pass Ua skipped | adjUaHigherPassTbMax necessity |
| Tb maxFinRes | Tb convergence quality | adjPseudoInvDtTb, tbHardCap* |
| Ua/Ub/Tb/ka/omegaa conv/sing | Solver convergence/singular | adjPseudoInvDt* |
| `tb_maxTAbs_reset`, `tb_sumT_reset` | Tb normalization resets | Tb stability |

---

## 3. Preliminary Classification

### 3.1 Likely ESSENTIAL (crash/stability safeguards)

| Option | Reason |
|--------|--------|
| **tbPreSolveReset** | Resets Tb when it blows up; CRASH_INVESTIGATION_opt313 identifies Tb divergence as root cause. Without it, Tb can drive ka/omegaa → turbSourceUa → FPE. |
| **tbRhsCapMax** | Caps Tb RHS to avoid FPE in PBiCGStab. Explicitly tied to crash at opt 313. |
| **adjTbGradTgradTbCapMax** | Caps grad(T)·grad(Tb) in ka/omegaa thermal source. Prevents explosion when Tb clipped. |
| **adjUaAllPassTbMax** | Skips Ua when Tb is bad; decouples Ua from divergent Tb. |

### 3.2 Likely IMPORTANT (numerical conditioning)

| Option | Reason |
|--------|--------|
| **adjPseudoInvDtTb** | Code default (10) differs from optProperties (100000). Higher values strongly improve Tb diagonal dominance; crash doc recommends 1e5–2e5. |
| **adjPseudoInvDtUa** | Affects Ua solver stability. |
| **adjPseudoInvDtUb** | Affects Ub solver stability. |
| **adjPseudoInvDtTurb** | Affects ka/omegaa solver stability. |

### 3.3 MAY BE REMOVABLE (sensible hardcoded defaults)

| Option | Reason |
|--------|--------|
| **tbHardCapMin / tbHardCapMax** | Used for per-cell Tb clipping; defaults (1e4–1e6) are reasonable. Could hardcode. |
| **adjUaHigherPassTbMax** | Similar role to adjUaAllPassTbMax but for higher passes. If adjUaAllPassTbMax is kept, this could be derived or hardcoded. |
| **adjTurbBetaMin, adjTurbRampStart, adjTurbRampEnd** | Continuation ramp. Sensible defaults; rarely tuned. |
| **adjTurbLagEvery** | Performance vs accuracy trade-off; 2 is typical. |
| **adjHeatBetaMin, adjHeatRampStart, adjHeatRampEnd** | Heat continuation ramp. |
| **adjUbPhiHbyAbCapFloor** | Ub cap floor; 1e-7 is standard. |
| **adjUaPassDecay, adjUaPassMinGain** | Ua multi-pass damping; defaults appear stable. |
| **adjUaLocalCapCoeff, adjUaRhsCapCoeff** | Turb source caps; scale with U²/L. Could hardcode. |

---

## 4. Next Steps

1. **Share solver_diagnostics logs**: `solver_diagnostics_rank*.log` and `solver_status.log` from representative runs.
2. **Count trigger frequency**:
   - How often does `tb_preSolve_reset` fire?
   - How often does `tb_rhs_capped` fire?
   - How often do `uaAllPass_skipped` / `uaHigherPass_skipped` fire?
3. **Refine classification**: If triggers are rare, more options can move to “removable” with hardcoded defaults.
4. **Remove options**: After evaluation, remove non-essential options from optProperties and use code defaults.

---

## 5. Proposed Minimal Set (to be confirmed with logs)

**Keep in optProperties** (user-tunable):

- `adjPseudoInvDtTb` — critical for Tb stability; user may need to raise for difficult cases.
- `tbPreSolveReset` — explicit crash safeguard.
- `tbRhsCapMax` — explicit FPE safeguard.
- `adjTbGradTgradTbCapMax` — ka/omegaa thermal source safeguard.
- `adjUaAllPassTbMax` — decoupling when Tb bad.

**Remove** (use hardcoded defaults in createFields.H):

- All others in sections 8–10, unless logs show they are frequently triggered or necessary for convergence.

---

*Last updated: preliminary; awaiting solver diagnostics.*
