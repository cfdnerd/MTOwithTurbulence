T=T2;
// Turbulent thermal diffusivity: alphat = nut/Prt (from thermalProperties)
volScalarField alphaEff(DT + turbulence->nut()/Prt);
// Floor alphaEff to bound thermal matrix condition number (avoids FPE in iterative solvers)
{
    scalar aEffAvg = gAverage(alphaEff.primitiveField());
    scalar aEffFloor = Foam::max(1e-6*aEffAvg, 1e-12);
    alphaEff = max(alphaEff, dimensionedScalar("aEffFloor", alphaEff.dimensions(), aEffFloor));
}
fvScalarMatrix TEqn
(
    fvm::div(phi,T)
   -fvm::laplacian(alphaEff, T)
);
if(DDHS)
{
 TEqn-=heatSource/rhoc*(1-xh);
}
else
{
 TEqn-=heatSource/rhoc;
}
SolverPerformance<scalar> spT = TEqn.solve();
logSolverStats(runTime, "Primal_T", 0, "T", spT, opt, runTime.timeName());

