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
// Floor rAU to bound pressure matrix condition number (avoids FPE in PBiCGStab)
const scalar rAUfloor = 1e-6 * gAverage(rAU);
rAU = max(rAU, dimensionedScalar("rAUfloor", rAU.dimensions(), max(rAUfloor, 1e-12)));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
MRF.makeRelative(phiHbyA);
adjustPhi(phiHbyA, U, p);

tmp<volScalarField> rAtU(rAU);

if (simple.consistent())
{
    rAtU = 1.0/(1.0/rAU - UEqn.H1());
    rAtU.ref() = max(rAtU(), dimensionedScalar("rAtUfloor", rAtU().dimensions(), max(rAUfloor, 1e-12)));
    phiHbyA +=
        fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
    HbyA -= (rAU - rAtU())*fvc::grad(p);
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



