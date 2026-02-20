scalar prevPaInitRes = 1e300;
label paDivergeCount = 0;
const label maxAdjUaIter = 15;
// Keep turbulence forcing growth bounded across inner iterations/opt steps.
scalar prevTurbSourceUaMax = -1.0;
volVectorField UaSafe(Ua);
volScalarField paSafe(pa);
surfaceScalarField phiaSafe(phia);
for(i=0;i<maxAdjUaIter;i++)
{ 
  const word uaStage(word("UaLoop_") + Foam::name(i));
  const scalar rampDen = Foam::max(scalar(adjTurbRampEnd - adjTurbRampStart), scalar(1));
  const scalar rampXi = Foam::min
  (
      scalar(1),
      Foam::max(scalar(0), scalar(opt - adjTurbRampStart)/rampDen)
  );
  const scalar turbAdjBeta = Foam::min(scalar(1), Foam::max(adjTurbBetaMin + (1.0 - adjTurbBetaMin)*rampXi, scalar(0)));
  const dimensionedScalar betaAdj("betaAdj", dimless, turbAdjBeta);
  const dimensionedScalar adjPseudoUa("adjPseudoUa", dimless/dimTime, adjPseudoInvDtUa);
  const volScalarField magUStage("magUStage", mag(U));
  const scalar UavgStage = Foam::max(gAverage(magUStage.primitiveField()), 1e-8);
  const scalar UmaxStage = Foam::max(gMax(magUStage.primitiveField()), UavgStage);
  const volScalarField magUaSafe("magUaSafe", mag(UaSafe));
  const volScalarField magPaSafe("magPaSafe", mag(paSafe));
  const scalar UaSafeAvg = Foam::max(gAverage(magUaSafe.primitiveField()), 1e-8);
  const scalar paSafeAvg = Foam::max(gAverage(magPaSafe.primitiveField()), 1e-8);
  // Cap with adjoint-history scale plus primal-based ceiling.
  const scalar UaHardCap = Foam::min(Foam::max(50.0*UaSafeAvg, 1e3), Foam::max(1e6*UmaxStage, 1e3));
  const scalar paHardCap = Foam::min(Foam::max(50.0*paSafeAvg, 1e4), Foam::max(1e6*UmaxStage*UmaxStage, 1e5));
  volVectorField adjointTransposeConvectiona((fvc::grad(U) & Ua));

  volVectorField turbSourceUa(fvc::grad(pa) * scalar(0));  // zero by default
  if (mesh.thisDb().foundObject<volScalarField>("k") && mesh.thisDb().foundObject<volScalarField>("omega"))
  {
      const volScalarField& ka = mesh.thisDb().lookupObject<volScalarField>("ka");
      const volScalarField& omegaa = mesh.thisDb().lookupObject<volScalarField>("omegaa");
      const volScalarField& k = mesh.thisDb().lookupObject<volScalarField>("k");
      const volScalarField& omega = mesh.thisDb().lookupObject<volScalarField>("omega");
      const volScalarField& nut = turbulence->nut();
      const dimensionedScalar gamma1("gamma1", dimless, 5.0/9.0);
      // Adaptive denominator floors to avoid division-driven source blowup
      const scalar kAvg = Foam::max(gAverage(k.primitiveField()), 1e-12);
      const scalar omegaAvg = Foam::max(gAverage(omega.primitiveField()), 1e-12);
      const scalar kFloorVal = Foam::max(1e-3*kAvg, 1e-8);
      const scalar omegaFloorVal = Foam::max(1e-3*omegaAvg, 1e-8);
      const dimensionedScalar kSmall("kSmall", dimArea/dimTime/dimTime, kFloorVal);
      const dimensionedScalar omegaSmall("omegaSmall", dimless/dimTime, omegaFloorVal);
      const volScalarField kDenom(max(k, kSmall));
      const volScalarField omegaDenom(max(omega, omegaSmall));
      // Adjoint momentum source: scale to [0 1 -2]. ka-term: div/(k); omegaa-term: div/(omega)*nut
      turbSourceUa =
          -fvc::div(ka*2*nut*twoSymm(fvc::grad(U)))/kDenom
          -nut*fvc::div(omegaa*2*gamma1*twoSymm(fvc::grad(U)))/omegaDenom;
      turbSourceUa *= betaAdj;

      // Robust bound independent of pa to avoid source-cap feedback blow-up:
      // cap by physical convective scale U^2/L and by historical growth.
      const scalar nCell = Foam::max(scalar(mesh.nCells()), scalar(1));
      const scalar volTot = Foam::max(gSum(mesh.V()), scalar(1e-30));
      const scalar hChar = Foam::pow(volTot/nCell, 1.0/3.0);
      const scalar convScale = Foam::max(UmaxStage*UmaxStage/Foam::max(hChar, 1e-8), 1e-8);
      const scalar srcCapPhys = Foam::max(1e3*convScale, 1e-2);
      const scalar srcCapHist = (prevTurbSourceUaMax > 0.0)
          ? Foam::max(5.0*prevTurbSourceUaMax, 1e-2)
          : srcCapPhys;
      const scalar srcMaxAbs = Foam::min(srcCapPhys, srcCapHist);
      scalar srcUsedMax = 0.0;
      forAll(turbSourceUa, cellI)
      {
          vector& sv = turbSourceUa[cellI];
          if (!(std::isfinite(sv.x()) && std::isfinite(sv.y()) && std::isfinite(sv.z())))
          {
              sv = vector::zero;
          }
          else
          {
              const scalar sm = mag(sv);
              if (sm > srcMaxAbs)
              {
                  sv *= (srcMaxAbs/(sm + VSMALL));
              }
              srcUsedMax = Foam::max(srcUsedMax, mag(sv));
          }
      }
      prevTurbSourceUaMax = srcUsedMax;
      MTOdiag::logVectorFieldHealth(runTime, opt, uaStage, "turbSourceUa", turbSourceUa.primitiveField());
      MTOdiag::logMetric(runTime, opt, uaStage, "turbSourceUa_capPhys", srcCapPhys);
      MTOdiag::logMetric(runTime, opt, uaStage, "turbSourceUa_capHist", srcCapHist);
      MTOdiag::logMetric(runTime, opt, uaStage, "turbSourceUa_capUsed", srcMaxAbs);
      MTOdiag::logMetric(runTime, opt, uaStage, "turbSourceUa_usedMax", srcUsedMax);
      MTOdiag::logMetric(runTime, opt, uaStage, "k_floor", kFloorVal);
      MTOdiag::logMetric(runTime, opt, uaStage, "omega_floor", omegaFloorVal);
      MTOdiag::logMetric(runTime, opt, uaStage, "turbAdjBeta", turbAdjBeta);
  }

  tmp<fvVectorMatrix> tUaEqn
  (
      fvm::div(-phi, Ua)
    + adjointTransposeConvectiona
    - fvm::laplacian(turbulence->nuEff(), Ua)
    + fvm::Sp(adjPseudoUa, Ua)
    + fvm::Sp(alpha, Ua)
   ==
      fvOptions(Ua)
  );
  fvVectorMatrix& UaEqn = tUaEqn.ref();

  UaEqn.relax();

  fvOptions.constrain(UaEqn);

  SolverPerformance<vector> spUa = solve(UaEqn == -fvc::grad(pa) + turbSourceUa);
  logSolverStats(runTime, "AdjointFlow_Ua", i, "Ua", spUa, opt, runTime.timeName());

  fvOptions.correct(Ua);
//****************************************
  volScalarField rAUa(1.0/UaEqn.A());
  const scalar rAUaAvg = Foam::max(gAverage(rAUa.primitiveField()), 1e-12);
  const scalar rAUaMin = 1e-6*rAUaAvg;
  const scalar rAUaMax = 1e6*rAUaAvg;
  forAll(rAUa, cellI)
  {
      scalar& rv = rAUa[cellI];
      if (!std::isfinite(rv))
      {
          rv = rAUaMin;
      }
      else
      {
          rv = Foam::min(Foam::max(rv, rAUaMin), rAUaMax);
      }
  }
  MTOdiag::logScalarFieldHealth(runTime, opt, uaStage, "rAUa", rAUa.primitiveField());
  volVectorField HbyAa("HbyAa", Ua);
  HbyAa = rAUa*UaEqn.H();
  surfaceScalarField phiHbyAa("phiHbyAa", fvc::flux(HbyAa));
  adjustPhi(phiHbyAa, Ua, pa);
  
  tmp<volScalarField> rAtUa(rAUa);
  bool useConsistentBranchUa = simple.consistent();
  if (simple.consistent())
  {
      rAtUa = 1.0/(1.0/rAUa - UaEqn.H1());
      label nBadRAtUa = 0;
      volScalarField& rAtUaRef = rAtUa.ref();
      forAll(rAtUaRef, cellI)
      {
          scalar& rv = rAtUaRef[cellI];
          if (!std::isfinite(rv))
          {
              rv = rAUaMin;
              nBadRAtUa++;
          }
          else
          {
              rv = Foam::min(Foam::max(rv, rAUaMin), rAUaMax);
          }
      }
      if (nBadRAtUa > 0)
      {
          useConsistentBranchUa = false;
          rAtUaRef = rAUa;
      }
  }
  if (useConsistentBranchUa)
  {
      phiHbyAa += fvc::interpolate(rAtUa() - rAUa)*fvc::snGrad(pa)*mesh.magSf();
      HbyAa -= (rAUa - rAtUa())*fvc::grad(pa);
  }
  MTOdiag::logScalarFieldHealth(runTime, opt, uaStage, "rAtUa", rAtUa().primitiveField());

  const scalar phiAaAvg = Foam::max(gAverage(mag(phiHbyAa.primitiveField())), 1e-12);
  const scalar phiAaMaxAbs = 1e6*phiAaAvg;
  forAll(phiHbyAa, faceI)
  {
      scalar& phiv = phiHbyAa[faceI];
      if (!std::isfinite(phiv))
      {
          phiv = 0.0;
      }
      else
      {
          phiv = Foam::min(Foam::max(phiv, -phiAaMaxAbs), phiAaMaxAbs);
      }
  }
  MTOdiag::logScalarFieldHealth(runTime, opt, uaStage, "phiHbyAa", phiHbyAa.primitiveField());
  tUaEqn.clear();

  // Update the pressure BCs to ensure flux consistency
  constrainPressure(pa, Ua, phiHbyAa, rAtUa(), MRF);           
  // Non-orthogonal pressure corrector loop
  bool badPaResidual = false;
  scalar paFirstInitRes = 1e300;
  scalar paLastFinalRes = 1e300;
  label nonOrthSolveI = 0;
  while (simple.correctNonOrthogonal())
  {
      fvScalarMatrix paEqn
      (
          fvm::laplacian(rAtUa(), pa) == fvc::div(phiHbyAa)
      );

      paEqn.setReference(paRefCell, paRefValue);
      SolverPerformance<scalar> sppa = paEqn.solve();
      logSolverStats(runTime, "AdjointFlow_Ua", i, "pa", sppa, opt, runTime.timeName());

      const scalar paInitRes = sppa.initialResidual();
      paLastFinalRes = sppa.finalResidual();
      if (nonOrthSolveI == 0)
      {
          paFirstInitRes = paInitRes;
      }
      nonOrthSolveI++;
      if (!std::isfinite(paInitRes) || !std::isfinite(paLastFinalRes))
      {
          badPaResidual = true;
          break;
      }
      MTOdiag::logMetric(runTime, opt, uaStage, "pa_initRes", paInitRes);
      MTOdiag::logMetric(runTime, opt, uaStage, "pa_finalRes", paLastFinalRes);

      if (simple.finalNonOrthogonalIter())
      {
          phia = phiHbyAa - paEqn.flux();
      }
  }

  if (badPaResidual)
  {
      Info<< "Warning: non-finite pa residual detected; restoring last safe Ua/pa state." << endl;
      Ua = UaSafe;
      pa = paSafe;
      phia = phiaSafe;
      break;
  }

  if (i > 0 && paFirstInitRes > 1.15*prevPaInitRes)
  {
      paDivergeCount++;
  }
  else
  {
      paDivergeCount = 0;
  }
  prevPaInitRes = paFirstInitRes;
  MTOdiag::logMetric(runTime, opt, uaStage, "pa_firstInitRes", paFirstInitRes);
  MTOdiag::logMetric(runTime, opt, uaStage, "pa_driftCount", scalar(paDivergeCount));
  if (paDivergeCount >= 2)
  {
      Info<< "Warning: pa first-pass residual drifting up; restoring last safe Ua/pa state and exiting inner loop." << endl;
      Ua = UaSafe;
      pa = paSafe;
      phia = phiaSafe;
      break;
  }

  //#include "adjointContinuityErrs.H"

  // Explicitly relax pressure for adjoint momentum corrector
  pa.relax();

  // Adjoint momentum corrector
  Ua = HbyAa - rAtUa()*fvc::grad(pa);
  Ua.correctBoundaryConditions();
  fvOptions.correct(Ua);

  // Guard against non-finite adjoint fields before next sub-iteration
  forAll(Ua, cellI)
  {
      vector& uv = Ua[cellI];
      if (!(std::isfinite(uv.x()) && std::isfinite(uv.y()) && std::isfinite(uv.z())))
      {
          uv = vector::zero;
      }
      else
      {
          const scalar um = mag(uv);
          if (um > UaHardCap)
          {
              uv *= (UaHardCap/(um + VSMALL));
          }
      }
  }
  forAll(pa, cellI)
  {
      scalar& pv = pa[cellI];
      if (!std::isfinite(pv))
      {
          pv = 0.0;
      }
      else if (pv > paHardCap)
      {
          pv = paHardCap;
      }
      else if (pv < -paHardCap)
      {
          pv = -paHardCap;
      }
  }
  MTOdiag::logVectorFieldHealth(runTime, opt, uaStage, "Ua", Ua.primitiveField());
  MTOdiag::logScalarFieldHealth(runTime, opt, uaStage, "pa", pa.primitiveField());
  const volScalarField magUaNow("magUaNow", mag(Ua));
  const volScalarField magPaNow("magPaNow", mag(pa));
  const scalar UaMaxNow = Foam::max(magUaNow.primitiveField());
  const scalar paAbsMaxNow = Foam::max(magPaNow.primitiveField());
  MTOdiag::logMetric(runTime, opt, uaStage, "Ua_hardCap", UaHardCap);
  MTOdiag::logMetric(runTime, opt, uaStage, "pa_hardCap", paHardCap);
  MTOdiag::logMetric(runTime, opt, uaStage, "Ua_absMax", UaMaxNow);
  MTOdiag::logMetric(runTime, opt, uaStage, "pa_absMax", paAbsMaxNow);
  if (!std::isfinite(UaMaxNow) || !std::isfinite(paAbsMaxNow) || UaMaxNow > 1.05*UaHardCap || paAbsMaxNow > 1.05*paHardCap)
  {
      Info<< "Warning: Ua/pa exceeded hard cap; restoring last safe Ua/pa state and exiting inner loop." << endl;
      Ua = UaSafe;
      pa = paSafe;
      phia = phiaSafe;
      break;
  }
  
  Ua.storePrevIter();
  pa.storePrevIter();
  phia.storePrevIter();
  UaSafe = Ua;
  paSafe = pa;
  phiaSafe = phia;

  if (paLastFinalRes < 5e-5 && i >= 2)
  {
      break;
  }
}