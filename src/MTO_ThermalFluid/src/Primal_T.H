T=T2;
// Turbulent thermal diffusivity: alphat = nut/Prt (IncompressibleTurbulenceModel has no alphat())
const dimensionedScalar Prt("Prt", dimless, 0.85);
volScalarField alphaEff(DT + turbulence->nut()/Prt);
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

