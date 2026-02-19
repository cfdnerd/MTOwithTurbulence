// Brinkman k/omega bounds: clip k and omega in high-alpha (porous/solid) regions for stability
// Called after turbulence->correct()
{
    if (mesh.thisDb().foundObject<volScalarField>("k") && mesh.thisDb().foundObject<volScalarField>("omega"))
    {
        volScalarField& kRef = const_cast<volScalarField&>(mesh.thisDb().lookupObject<volScalarField>("k"));
        volScalarField& omegaRef = const_cast<volScalarField&>(mesh.thisDb().lookupObject<volScalarField>("omega"));
        const scalar kMin = 1e-10;
        const scalar omegaMin = 1e-6;
        const scalar omegaMax = 1e10;
        const scalar alphaThreshold = 0.1 * alphaMax.value();
        forAll(alpha, celli)
        {
            if (alpha[celli] > alphaThreshold)
            {
                kRef[celli] = Foam::max(kRef[celli], kMin);
                omegaRef[celli] = Foam::max(Foam::min(omegaRef[celli], omegaMax), omegaMin);
            }
        }
        kRef.correctBoundaryConditions();
        omegaRef.correctBoundaryConditions();
    }
}
