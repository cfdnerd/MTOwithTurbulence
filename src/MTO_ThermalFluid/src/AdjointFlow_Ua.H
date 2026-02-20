scalar prevPaInitRes = 1e300;
label paDivergeCount = 0;
const label maxAdjUaIter = 15;
volVectorField UaSafe(Ua);
volScalarField paSafe(pa);
surfaceScalarField phiaSafe(phia);
for(i=0;i<maxAdjUaIter;i++)
{ 
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
      const dimensionedScalar kSmall("kSmall", dimArea/dimTime/dimTime, 1e-6);
      const dimensionedScalar omegaSmall("omegaSmall", dimless/dimTime, 1e-6);
      // Adjoint momentum source: scale to [0 1 -2]. ka-term: div/(k); omegaa-term: div/(omega)*nut
      turbSourceUa = -fvc::div(ka*2*nut*twoSymm(fvc::grad(U)))/(k + kSmall)
                     -nut*fvc::div(omegaa*2*gamma1*twoSymm(fvc::grad(U)))/(omega + omegaSmall);

      // Clamp turbSourceUa to avoid local spikes that can destabilize Ua/pa solves
      const scalar srcAvg = Foam::max(gAverage(mag(turbSourceUa.primitiveField())), 1e-12);
      const scalar srcMaxAbs = 1e5*srcAvg;
      forAll(turbSourceUa, cellI)
      {
          vector& sv = turbSourceUa[cellI];
          if (!(std::isfinite(sv.x()) && std::isfinite(sv.y()) && std::isfinite(sv.z())))
          {
              sv = vector::zero;
          }
          else
          {
              sv.x() = Foam::min(Foam::max(sv.x(), -srcMaxAbs), srcMaxAbs);
              sv.y() = Foam::min(Foam::max(sv.y(), -srcMaxAbs), srcMaxAbs);
              sv.z() = Foam::min(Foam::max(sv.z(), -srcMaxAbs), srcMaxAbs);
          }
      }
  }

  tmp<fvVectorMatrix> tUaEqn
  (
      fvm::div(-phi, Ua)
    + adjointTransposeConvectiona
    - fvm::laplacian(turbulence->nuEff(), Ua)
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
  }
  forAll(pa, cellI)
  {
      scalar& pv = pa[cellI];
      if (!std::isfinite(pv))
      {
          pv = 0.0;
      }
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