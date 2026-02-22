// Adjoint k-omega SST turbulence equations (sensitivity-consistent)
// Solves adjoint ka, omegaa when RAS turbulence is active
// See Kavvadias et al. (2014) for full derivation
{
    if (mesh.thisDb().foundObject<volScalarField>("k") && mesh.thisDb().foundObject<volScalarField>("omega"))
    {
        const scalar rampDen = Foam::max(scalar(adjTurbRampEnd - adjTurbRampStart), scalar(1));
        const scalar rampXi = Foam::min
        (
            scalar(1),
            Foam::max(scalar(0), scalar(opt - adjTurbRampStart)/rampDen)
        );
        const scalar turbAdjBeta = Foam::min(scalar(1), Foam::max(adjTurbBetaMin + (1.0 - adjTurbBetaMin)*rampXi, scalar(0)));
        const dimensionedScalar betaAdj("betaAdj", dimless, turbAdjBeta);
        const dimensionedScalar adjPseudoTurb("adjPseudoTurb", dimless/dimTime, adjPseudoInvDtTurb);
        const volScalarField& k = mesh.thisDb().lookupObject<volScalarField>("k");
        const volScalarField& omega = mesh.thisDb().lookupObject<volScalarField>("omega");
        const volScalarField& nut = turbulence->nut();
        const volScalarField nutSafe(max(nut, dimensionedScalar("zero", nut.dimensions(), 0.0)));

        // SST constants (must match primal kOmegaSST)
        const dimensionedScalar alphaK1("alphaK1", dimless, 0.85);
        const dimensionedScalar alphaOmega1("alphaOmega1", dimless, 0.5);
        const dimensionedScalar betaStar("betaStar", dimless, 0.09);
        const dimensionedScalar beta2("beta2", dimless, 0.0828);
        const dimensionedScalar kSmall("kSmall", dimArea/dimTime/dimTime, 1e-10);
        const dimensionedScalar omegaSmall("omegaSmall", dimless/dimTime, 1e-10);

        volScalarField dNutdk(nutSafe/(k + kSmall));
        volScalarField dNutdOmega(-nutSafe/(omega + omegaSmall));

        // Effective diffusivities for k and omega (use nutSafe to avoid negative diffusivity)
        volScalarField DkEff("DkEff", nu + alphaK1*nutSafe);
        volScalarField DomegaEff("DomegaEff", nu + alphaOmega1*nutSafe);

        // RHS: source from momentum (d(nut)/dk * dev2(twoSymm(grad(U))) : grad(Ua))
        // Raw term has dim [0 0 -1]; ka eqn needs [0 2 -3], omegaa eqn needs [0 0 -2]
        // Scale: ka by k, omegaa by omega (to match equation dimensions)
        volScalarField kaRhsMomentum
        (
            (k + kSmall)*(-dNutdk*(dev2(twoSymm(fvc::grad(U))) && fvc::grad(Ua)))
        );
        volScalarField omegaaRhsMomentum
        (
            (omega + omegaSmall)*(omega + omegaSmall)/(k + kSmall)
           *(-dNutdOmega*(dev2(twoSymm(fvc::grad(U))) && fvc::grad(Ua)))
        );

        // RHS: thermal source (objFunction 1) - alphat = nut/Prt affects T via diffusion
        // Scale by k²/Tref² for ka ([0 2 -3]) and omega²/Tref² for omegaa ([0 0 -2])
        // Cap gradTgradTb to avoid blow-up when Tb is clipped (crash-fix opt313)
        volScalarField kaRhsThermal(kaRhsMomentum * 0.0);
        volScalarField omegaaRhsThermal(omegaaRhsMomentum * 0.0);
        if (objFunction == 1)
        {
            const dimensionedScalar Tref2("Tref2", dimTemperature*dimTemperature, max(1.0, gMax(T)*gMax(T)));
            volScalarField gradTgradTb("gradTgradTb", fvc::grad(T) & fvc::grad(Tb));
            forAll(gradTgradTb, cellI)
            {
                scalar& g = gradTgradTb[cellI];
                if (!std::isfinite(g))
                    g = 0.0;
                else if (Foam::mag(g) > adjTbGradTgradTbCapMax)
                    g = (g > 0 ? adjTbGradTgradTbCapMax : -adjTbGradTgradTbCapMax);
            }
            kaRhsThermal = (k + kSmall)*(k + kSmall)/Tref2 * (dNutdk/Prt) * gradTgradTb;
            omegaaRhsThermal = (omega + omegaSmall)*(omega + omegaSmall)/Tref2
                             * (dNutdOmega/Prt) * gradTgradTb;
        }

        // Adjoint k equation: -div(phi,ka) - div(DkEff*grad(ka)) + beta*omega*ka = RHS
        fvScalarMatrix kaEqn
        (
            - fvm::div(phi, ka)
            - fvm::laplacian(DkEff, ka)
            + fvm::Sp(adjPseudoTurb, ka)
            + fvm::Sp(betaStar*omega, ka)
         ==
            betaAdj*(kaRhsMomentum + kaRhsThermal)
        );
        kaEqn.relax();
        SolverPerformance<scalar> spKa = kaEqn.solve();
        logSolverStats(runTime, "Adjoint_kOmegaSST", 0, "ka", spKa, opt, runTime.timeName());

        // Adjoint omega equation: -div(phi,omegaa) - div(DomegaEff*grad(omegaa)) + 2*beta2*omega*omegaa = RHS
        fvScalarMatrix omegaaEqn
        (
            - fvm::div(phi, omegaa)
            - fvm::laplacian(DomegaEff, omegaa)
            + fvm::Sp(adjPseudoTurb, omegaa)
            + fvm::Sp(scalar(2)*beta2.value()*omega, omegaa)
         ==
            betaAdj*(omegaaRhsMomentum + omegaaRhsThermal)
        );
        omegaaEqn.relax();
        SolverPerformance<scalar> spOmegaa = omegaaEqn.solve();
        logSolverStats(runTime, "Adjoint_kOmegaSST", 0, "omegaa", spOmegaa, opt, runTime.timeName());
        MTOdiag::logMetric(runTime, opt, "Adjoint_kOmegaSST", "turbAdjBeta", turbAdjBeta);
        // Field health for adjoint turbulence (diagLevel 2)
        MTOdiag::logScalarFieldHealth(runTime, opt, "Adjoint_kOmegaSST", "ka", ka.primitiveField());
        MTOdiag::logScalarFieldHealth(runTime, opt, "Adjoint_kOmegaSST", "omegaa", omegaa.primitiveField());
    }
}
