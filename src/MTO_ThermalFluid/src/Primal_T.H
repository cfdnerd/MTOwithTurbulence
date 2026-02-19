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
 TEqn-=100000000.0/4180000*Q*(1-xh);
}
else
{
 TEqn-=100000000.0/4180000*Q;
}
TEqn.solve();

