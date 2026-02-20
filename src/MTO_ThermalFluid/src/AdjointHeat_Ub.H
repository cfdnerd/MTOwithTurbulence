for(i=0;i<50;i++)
{    
volVectorField adjointTransposeConvectionb((fvc::grad(U) & Ub));

tmp<fvVectorMatrix> tUbEqn
(
    fvm::div(-phi, Ub)
  + adjointTransposeConvectionb
  - fvm::laplacian(turbulence->nuEff(), Ub)
  + fvm::Sp(alpha, Ub)
  + b1*Tb*fvc::grad(T)
 ==
    fvOptions(Ub)
);
fvVectorMatrix& UbEqn = tUbEqn.ref();

UbEqn.relax();

fvOptions.constrain(UbEqn);

SolverPerformance<vector> spUb = solve(UbEqn == -fvc::grad(pb));
logSolverStats(runTime, "AdjointHeat_Ub", i, "Ub", spUb, opt, runTime.timeName());

fvOptions.correct(Ub);
//****************************************
volScalarField rAUb(1.0/UbEqn.A());
const scalar rAUbAvg = Foam::max(gAverage(rAUb.primitiveField()), 1e-12);
const scalar rAUbMin = 1e-6*rAUbAvg;
const scalar rAUbMax = 1e6*rAUbAvg;
forAll(rAUb, cellI)
{
    scalar& rv = rAUb[cellI];
    if (!Foam::isfinite(rv))
    {
        rv = rAUbMin;
    }
    else
    {
        rv = Foam::min(Foam::max(rv, rAUbMin), rAUbMax);
    }
}
volVectorField HbyAb("HbyAb", Ub);
HbyAb = rAUb*UbEqn.H();
surfaceScalarField phiHbyAb("phiHbyAb", fvc::flux(HbyAb));
adjustPhi(phiHbyAb, Ub, pb);

tmp<volScalarField> rAtUb(rAUb);

bool useConsistentBranchUb = simple.consistent();
if (simple.consistent())
{
    rAtUb = 1.0/(1.0/rAUb - UbEqn.H1());
    label nBadRAtUb = 0;
    volScalarField& rAtUbRef = rAtUb.ref();
    forAll(rAtUbRef, cellI)
    {
        scalar& rv = rAtUbRef[cellI];
        if (!Foam::isfinite(rv))
        {
            rv = rAUbMin;
            nBadRAtUb++;
        }
        else
        {
            rv = Foam::min(Foam::max(rv, rAUbMin), rAUbMax);
        }
    }
    if (nBadRAtUb > 0)
    {
        useConsistentBranchUb = false;
        rAtUbRef = rAUb;
    }
}
if (useConsistentBranchUb)
{
    phiHbyAb += fvc::interpolate(rAtUb() - rAUb)*fvc::snGrad(pb)*mesh.magSf();
    HbyAb -= (rAUb - rAtUb())*fvc::grad(pb);
}
const scalar phiAbAvg = Foam::max(gAverage(mag(phiHbyAb)), 1e-12);
const scalar phiAbMaxAbs = 1e6*phiAbAvg;
forAll(phiHbyAb, faceI)
{
    scalar& phiv = phiHbyAb[faceI];
    if (!Foam::isfinite(phiv))
    {
        phiv = 0.0;
    }
    else
    {
        phiv = Foam::min(Foam::max(phiv, -phiAbMaxAbs), phiAbMaxAbs);
    }
}
tUbEqn.clear();

// Update the pressure BCs to ensure flux consistency
constrainPressure(pb, Ub, phiHbyAb, rAtUb(), MRF);  
// Non-orthogonal pressure corrector loop
while (simple.correctNonOrthogonal())
{
    fvScalarMatrix pbEqn
    (
        fvm::laplacian(rAtUb(), pb) == fvc::div(phiHbyAb)
    );

    pbEqn.setReference(pbRefCell, pbRefValue);
    SolverPerformance<scalar> sppb = pbEqn.solve();
    logSolverStats(runTime, "AdjointHeat_Ub", i, "pb", sppb, opt, runTime.timeName());

    if (simple.finalNonOrthogonalIter())
    {
        phib = phiHbyAb - pbEqn.flux();
    }
}

// Guard against non-finite adjoint fields before correction
forAll(Ub, cellI)
{
    vector& uv = Ub[cellI];
    if (!(Foam::isfinite(uv.x()) && Foam::isfinite(uv.y()) && Foam::isfinite(uv.z())))
    {
        uv = vector::zero;
    }
}
forAll(pb, cellI)
{
    scalar& pv = pb[cellI];
    if (!Foam::isfinite(pv))
    {
        pv = 0.0;
    }
}

//#include "adjointContinuityErrs.H"

// Explicitly relax pressure for adjoint momentum corrector
pb.relax();

// Adjoint momentum corrector
Ub = HbyAb - rAtUb()*fvc::grad(pb);
Ub.correctBoundaryConditions();
fvOptions.correct(Ub);

Ub.storePrevIter();
pb.storePrevIter();
phib.storePrevIter();

if
(
    i >= 2
 && spUb.finalResidual() < 1e-4
)
{
    break;
}
}