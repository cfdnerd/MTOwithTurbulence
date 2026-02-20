# MTO_ThermalFluid Current Status (Single Reference)

Last updated: 2026-02-19  
Scope: `src/MTO_ThermalFluid` (OpenFOAM 6 branch)  

## 1) What this solver currently is

`MTO_ThermalFluid` is a density-based topology optimization solver for forced convection/thermal-fluid design with:
- Primal SIMPLE flow + heat solve
- Adjoint solves for power and heat objectives
- Optional turbulent closure based on primal `k-omega SST`
- MMA design update loop with PDE/filter chain

Main loop order (implemented):
1. `update.H`
2. `turbulence->correct()`
3. `Primal_U.H`
4. `Primal_T.H`
5. `AdjointHeat_Tb.H`
6. Adjoint coupling passes (`nAdjTurbPasses`):
   - `Adjoint_kOmegaSST.H` (with lag policy)
   - `AdjointHeat_Ub.H`
   - `AdjointFlow_Ua.H`
7. `costfunction.H`
8. `sensitivity.H`

## 2) Turbulence model status (primal and adjoint)

### Primal turbulence (current)
- Active model: incompressible RAS `k-omega SST` via `turbulence->correct()`.
- Primal momentum uses turbulence stress through OpenFOAM turbulence model path.
- Primal thermal diffusion uses turbulent contribution:
  - `alphaEff = DT + nut/Prt` in `Primal_T.H`.
  - `alphaEff` is floored each iteration to avoid thermal matrix ill-conditioning.

### Adjoint turbulence (current)
- Adjoint turbulence fields `ka` and `omegaa` are implemented and solved in `Adjoint_kOmegaSST.H`.
- Adjoint turbulence includes:
  - pseudo-transient diagonal stabilization
  - continuation factor (`betaAdj`) ramped by optimization iteration
  - momentum-coupled RHS terms using `dNut/dk`, `dNut/dOmega`
  - thermal RHS coupling when heat objective path is active
- Turbulence-adjoint solve frequency is lag-controlled (`adjTurbLagEvery`) and can be skipped on selected outer iterations except pass 0.

### Not implemented / deferred in turbulence
- Curvature-corrected SST (`kOmegaSSTCC`) is not integrated.
- Full adjoint wall-function treatment for turbulence adjoints is not implemented.
- The implementation is robust/stabilized for optimization use; it is not a symbolic exact discrete-adjoint of every turbulence closure term.

## 3) Objectives and adjoint state

Implemented adjoint fields:
- `Ua`, `pa` (power dissipation path)
- `Ub`, `pb`, `Tb` (thermal path)
- `ka`, `omegaa` (turbulence adjoint path)

Implemented objective modes:
- Objective 1: mean-temperature style path
- Objective 2: p-norm temperature path (`Pnorm`)

## 4) Stability and convergence frameworks currently implemented

The codebase has a layered stabilization framework designed for turbulent topology optimization:

### A) Pseudo-transient continuation
- Added as `fvm::Sp(adjPseudo..., field)` in:
  - `Ua`, `Ub`, `Tb`, `ka`, `omegaa` equations.
- Controls:
  - `adjPseudoInvDtUa`
  - `adjPseudoInvDtUb`
  - `adjPseudoInvDtTb`
  - `adjPseudoInvDtTurb`

### B) Continuation/homotopy ramps
- Turbulence-adjoint continuation:
  - `adjTurbBetaMin`, `adjTurbRampStart`, `adjTurbRampEnd`
- Heat-adjoint continuation:
  - `adjHeatBetaMin`, `adjHeatRampStart`, `adjHeatRampEnd`

### C) Lagged adjoint turbulence coupling
- Adjoint turbulence solve can be lagged:
  - `adjTurbLagEvery`
- Pass scheduling in main loop via `doTurbAdjSolveThisOpt`.

### D) Trust-region caps and finite guards
- Broad finite checks and clamping on:
  - reciprocal diagonals (`rAU*`, `rAt*`)
  - pressure fluxes (`phiHbyA*`)
  - adjoint fields (`Ua`, `Ub`, `Tb`, `pa`, `pb`)
- Hard-cap rollback to last safe field states when limits are violated.
- Residual drift detectors for `pa`/`pb` to prevent runaway iterations.

### E) Ua-specific robust source/RHS control (latest)
- `turbSourceUa` pass-dependent damping:
  - `adjUaPassDecay`, `adjUaPassMinGain`
- Local physical cap on `turbSourceUa` magnitude:
  - `adjUaLocalCapCoeff` (convective-scale based)
- Pre-solve Ua RHS sanitization and cap:
  - `adjUaRhsCapCoeff`
- Tb-based Ua trust gates:
  - all-pass gate: `adjUaAllPassTbMax`
  - higher-pass gate: `adjUaHigherPassTbMax`
- Safe-state restore and skip behavior when trust gates trip.

### F) Diagnostic heartbeat framework
- Fine-grained phase heartbeats and metric logging for localization of failures by iteration/pass/rank.
- Diagnostic controls:
  - `diagLevel`, `diagEvery`, `diagRank`, `diagOptStart`

## 5) Current numerical behavior (latest observed state)

- Recent runs progressed beyond earlier stop points and reached iteration ~141.
- Higher-pass Ua guard is actively triggering when `|Tb|` becomes extreme.
- Remaining instability signature: pass-0 `Ua` can still become stiff when `Tb` and `Ub` are highly saturated; new all-pass gate and RHS capping are in place to mitigate this.
- This is currently a robustness-prioritized optimization workflow for difficult turbulent regimes.

## 6) Files most relevant to turbulence and stability

- Main orchestration:
  - `src/MTO_ThermalFluid/src/MTO_ThermalFluid.C`
- Primal:
  - `src/MTO_ThermalFluid/src/Primal_U.H`
  - `src/MTO_ThermalFluid/src/Primal_T.H`
- Adjoint:
  - `src/MTO_ThermalFluid/src/AdjointFlow_Ua.H`
  - `src/MTO_ThermalFluid/src/AdjointHeat_Ub.H`
  - `src/MTO_ThermalFluid/src/AdjointHeat_Tb.H`
  - `src/MTO_ThermalFluid/src/Adjoint_kOmegaSST.H`
- Controls:
  - `src/MTO_ThermalFluid/src/createFields.H`
  - `src/MTO_ThermalFluid/app/constant/optProperties`

## 7) Build and run (OpenFOAM 6)

- Source OpenFOAM 6 environment.
- Build:
  - `cd $MTO_SOURCE/MTO_ThermalFluid`
  - `wmake`
- Use the same MPI family/version across OpenFOAM, PETSc, and runtime (`mpirun`) to avoid ABI issues.

## 8) Practical limitations and next priorities

1. Keep robustness controls enabled for turbulent optimization runs with strong penalization.
2. If convergence stalls, tune in this order:
   - `adjUaAllPassTbMax`, `adjUaRhsCapCoeff`, `adjUaHigherPassTbMax`
   - then continuation ramps and pseudo-transient strengths.
3. For future physics fidelity upgrades:
   - curvature-corrected SST integration
   - improved adjoint turbulence wall treatment
   - tighter consistency checks against finite-difference sensitivities on reduced test cases.

