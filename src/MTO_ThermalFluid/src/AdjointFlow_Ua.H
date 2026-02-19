for(i=0;i<40;i++)
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
      const dimensionedScalar kSmall("kSmall", dimArea/dimTime/dimTime, 1e-10);
      const dimensionedScalar omegaSmall("omegaSmall", dimless/dimTime, 1e-10);
      // Adjoint momentum source: raw terms have [0 3 -4] and [0 1 -3]; scale by 1/k and 1/omega to get [0 1 -2]
      turbSourceUa = -fvc::div(ka*2*nut*twoSymm(fvc::grad(U)))/(k + kSmall)
                     -fvc::div(omegaa*2*gamma1*twoSymm(fvc::grad(U)))/(omega + omegaSmall);
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
  volVectorField HbyAa("HbyAa", Ua);
  HbyAa = rAUa*UaEqn.H();
  tUaEqn.clear();          
  surfaceScalarField phiHbyAa("phiHbyAa", fvc::flux(HbyAa));
  adjustPhi(phiHbyAa, Ua, pa);
  
  tmp<volScalarField> rAtUa(rAUa);
  
   if (simple.consistent())
  {
      rAtUa = 1.0/(1.0/rAUa - UaEqn.H1());
      phiHbyAa +=
          fvc::interpolate(rAtUa() - rAUa)*fvc::snGrad(pa)*mesh.magSf();
      HbyAa -= (rAUa - rAtUa())*fvc::grad(pa);
  }

  tUaEqn.clear();

  // Update the pressure BCs to ensure flux consistency
  constrainPressure(pa, Ua, phiHbyAa, rAtUa(), MRF);           
  // Non-orthogonal pressure corrector loop
  while (simple.correctNonOrthogonal())
  {
      fvScalarMatrix paEqn
      (
          fvm::laplacian(rAUa, pa) == fvc::div(phiHbyAa)
      );

      paEqn.setReference(paRefCell, paRefValue);
      SolverPerformance<scalar> sppa = paEqn.solve();
      logSolverStats(runTime, "AdjointFlow_Ua", i, "pa", sppa, opt, runTime.timeName());

      if (simple.finalNonOrthogonalIter())
      {
          phia = phiHbyAa - paEqn.flux();
      }
  }

  //#include "adjointContinuityErrs.H"

  // Explicitly relax pressure for adjoint momentum corrector
  pa.relax();

  // Adjoint momentum corrector
  Ua = HbyAa - rAUa*fvc::grad(pa);
  Ua.correctBoundaryConditions();
  fvOptions.correct(Ua);
  
  Ua.storePrevIter();
  pa.storePrevIter();
  phia.storePrevIter();
}