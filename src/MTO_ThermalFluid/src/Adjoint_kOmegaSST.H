// Adjoint k-omega SST turbulence equations (sensitivity-consistent)
// Solves adjoint ka, omegaa when RAS turbulence is active
// See Kavvadias et al. (2014) for full derivation
{
    if (mesh.thisDb().foundObject<volScalarField>("k") && mesh.thisDb().foundObject<volScalarField>("omega"))
    {
        const volScalarField& k = mesh.thisDb().lookupObject<volScalarField>("k");
        const volScalarField& omega = mesh.thisDb().lookupObject<volScalarField>("omega");
        const volScalarField& nut = turbulence->nut();

        // SST constants (must match primal kOmegaSST)
        const dimensionedScalar alphaK1("alphaK1", dimless, 0.85);
        const dimensionedScalar alphaOmega1("alphaOmega1", dimless, 0.5);
        const dimensionedScalar betaStar("betaStar", dimless, 0.09);
        const dimensionedScalar beta2("beta2", dimless, 0.0828);
        const dimensionedScalar Prt("Prt", dimless, 0.85);
        const dimensionedScalar kSmall("kSmall", dimArea/dimTime/dimTime, 1e-10);
        const dimensionedScalar omegaSmall("omegaSmall", dimless/dimTime, 1e-10);

        volScalarField dNutdk(nut/(k + kSmall));
        volScalarField dNutdOmega(-nut/(omega + omegaSmall));

        // Effective diffusivities for k and omega
        volScalarField DkEff("DkEff", nu + alphaK1*nut);
        volScalarField DomegaEff("DomegaEff", nu + alphaOmega1*nut);

        // RHS: source from momentum (d(nut)/dk * dev2(twoSymm(grad(U))) : grad(Ua))
        volScalarField kaRhsMomentum(-dNutdk*(dev2(twoSymm(fvc::grad(U))) && fvc::grad(Ua)));
        volScalarField omegaaRhsMomentum(-dNutdOmega*(dev2(twoSymm(fvc::grad(U))) && fvc::grad(Ua)));

        // RHS: source from thermal (objFunction 1) - alphat = nut/Prt, d(alphaEff)/dk = dNutdk/Prt
        volScalarField kaRhsThermal(kaRhsMomentum * 0.0);
        volScalarField omegaaRhsThermal(omegaaRhsMomentum * 0.0);
        if (objFunction == 1)
        {
            kaRhsThermal = (dNutdk/Prt)*(fvc::grad(T) & fvc::grad(Tb));
            omegaaRhsThermal = (dNutdOmega/Prt)*(fvc::grad(T) & fvc::grad(Tb));
        }

        // Adjoint k equation: -div(phi,ka) - div(DkEff*grad(ka)) + beta*omega*ka = RHS
        fvScalarMatrix kaEqn
        (
            - fvm::div(phi, ka)
            - fvm::laplacian(DkEff, ka)
            + fvm::Sp(betaStar*omega, ka)
         ==
            kaRhsMomentum + kaRhsThermal
        );
        kaEqn.relax();
        SolverPerformance<scalar> spKa = kaEqn.solve();
        logSolverStats(runTime, "Adjoint_kOmegaSST", 0, "ka", spKa, opt, runTime.timeName());

        // Adjoint omega equation: -div(phi,omegaa) - div(DomegaEff*grad(omegaa)) + 2*beta2*omega*omegaa = RHS
        fvScalarMatrix omegaaEqn
        (
            - fvm::div(phi, omegaa)
            - fvm::laplacian(DomegaEff, omegaa)
            + fvm::Sp(scalar(2)*beta2.value()*omega, omegaa)
         ==
            omegaaRhsMomentum + omegaaRhsThermal
        );
        omegaaEqn.relax();
        SolverPerformance<scalar> spOmegaa = omegaaEqn.solve();
        logSolverStats(runTime, "Adjoint_kOmegaSST", 0, "omegaa", spOmegaa, opt, runTime.timeName());
    }
}
