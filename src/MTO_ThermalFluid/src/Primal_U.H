for(i=0;i<80;i++)
{
MRF.correctBoundaryVelocity(U);
tmp<fvVectorMatrix> tUEqn
(
    fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
  + fvm::Sp(alpha, U)
 ==
    fvOptions(U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (simple.momentumPredictor())
{
    SolverPerformance<vector> spU = solve(UEqn == -fvc::grad(p));
    logSolverStats(runTime, "Primal_U", i, "U", spU, opt, runTime.timeName());

    fvOptions.correct(U);
}
//****************************************//
volScalarField rAU(1.0/UEqn.A());
// Sanitize and bound rAU to keep pressure matrix scaling stable
const scalar rAUavg = Foam::max(gAverage(rAU.primitiveField()), 1e-12);
const scalar rAUmin = 1e-6*rAUavg;
const scalar rAUmax = 1e6*rAUavg;
forAll(rAU, cellI)
{
    scalar& rv = rAU[cellI];
    if (!std::isfinite(rv))
    {
        rv = rAUmin;
    }
    else
    {
        rv = Foam::min(Foam::max(rv, rAUmin), rAUmax);
    }
}
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
MRF.makeRelative(phiHbyA);
adjustPhi(phiHbyA, U, p);

tmp<volScalarField> rAtU(rAU);

bool useConsistentBranch = simple.consistent();
if (simple.consistent())
{
    rAtU = 1.0/(1.0/rAU - UEqn.H1());
    label nBadRAtU = 0;
    volScalarField& rAtURef = rAtU.ref();
    forAll(rAtURef, cellI)
    {
        scalar& rv = rAtURef[cellI];
        if (!std::isfinite(rv))
        {
            rv = rAUmin;
            nBadRAtU++;
        }
        else
        {
            rv = Foam::min(Foam::max(rv, rAUmin), rAUmax);
        }
    }
    if (nBadRAtU > 0)
    {
        useConsistentBranch = false;
        rAtURef = rAU;
    }
}
if (useConsistentBranch)
{
    phiHbyA += fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
    HbyA -= (rAU - rAtU())*fvc::grad(p);
}

// Sanitize and clamp phiHbyA before pressure solve
const scalar phiAvg = Foam::max(gAverage(mag(phiHbyA.primitiveField())), 1e-12);
const scalar phiMaxAbs = 1e6*phiAvg;
forAll(phiHbyA, faceI)
{
    scalar& phiv = phiHbyA[faceI];
    if (!std::isfinite(phiv))
    {
        phiv = 0.0;
    }
    else
    {
        phiv = Foam::min(Foam::max(phiv, -phiMaxAbs), phiMaxAbs);
    }
}

tUEqn.clear();

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAtU(), MRF);

// Non-orthogonal pressure corrector loop
while (simple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);

    SolverPerformance<scalar> spp = pEqn.solve();
    logSolverStats(runTime, "Primal_U", i, "p", spp, opt, runTime.timeName());

    if (simple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

// Momentum corrector
U = HbyA - rAtU()*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);

U.storePrevIter();
p.storePrevIter();
phi.storePrevIter();

}



