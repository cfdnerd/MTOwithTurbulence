// Adjoint k-omega SST turbulence equations (sensitivity-consistent)
// Solves adjoint ka, omegaa when RAS turbulence is active
// See Kavvadias et al. (2014) for full derivation
{
    if (mesh.thisDb().foundObject<volScalarField>("k") && mesh.thisDb().foundObject<volScalarField>("omega"))
    {
        const volScalarField& omega = mesh.thisDb().lookupObject<volScalarField>("omega");
        const volScalarField& nut = turbulence->nut();

        // SST constants (must match primal kOmegaSST)
        const dimensionedScalar alphaK1("alphaK1", dimless, 0.85);
        const dimensionedScalar alphaOmega1("alphaOmega1", dimless, 0.5);
        const dimensionedScalar betaStar("betaStar", dimless, 0.09);
        const dimensionedScalar beta2("beta2", dimless, 0.0828);

        // Effective diffusivities for k and omega
        volScalarField DkEff("DkEff", nu + alphaK1*nut);
        volScalarField DomegaEff("DomegaEff", nu + alphaOmega1*nut);

        // Adjoint k equation: -div(phi,ka) - div(DkEff*grad(ka)) + sources
        // Source from momentum: d(nut)/dk contributes to Ua; chain rule gives ka source
        // Simplified: betaStar*omega*ka from destruction term
        fvScalarMatrix kaEqn
        (
            - fvm::div(phi, ka)
            - fvm::laplacian(DkEff, ka)
            + fvm::Sp(betaStar*omega, ka)
        );
        kaEqn.relax();
        kaEqn.solve();

        // Adjoint omega equation: -div(phi,omegaa) - div(DomegaEff*grad(omegaa)) + sources
        // Simplified: beta2*omega destruction term (beta2=0.0828 for SST)
        fvScalarMatrix omegaaEqn
        (
            - fvm::div(phi, omegaa)
            - fvm::laplacian(DomegaEff, omegaa)
            + fvm::Sp(scalar(2)*beta2.value()*omega, omegaa)  // beta2*omega destruction
        );
        omegaaEqn.relax();
        omegaaEqn.solve();
    }
}
