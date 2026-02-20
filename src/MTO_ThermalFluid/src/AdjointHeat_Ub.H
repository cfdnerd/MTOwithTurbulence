scalar prevPbInitRes = 1e300;
label pbDivergeCount = 0;
const label maxAdjUbIter = 15;
scalar prevPhiHbyAbMax = -1.0;
volVectorField UbSafe(Ub);
volScalarField pbSafe(pb);
surfaceScalarField phibSafe(phib);
for(i=0;i<maxAdjUbIter;i++)
{    
const word ubStage(word("UbLoop_") + Foam::name(i));
const scalar heatRampDen = Foam::max(scalar(adjHeatRampEnd - adjHeatRampStart), scalar(1));
const scalar heatRampXi = Foam::min
(
    scalar(1),
    Foam::max(scalar(0), scalar(opt - adjHeatRampStart)/heatRampDen)
);
const scalar heatAdjBeta = Foam::min(scalar(1), Foam::max(adjHeatBetaMin + (1.0 - adjHeatBetaMin)*heatRampXi, scalar(0)));
const dimensionedScalar betaHeatAdj("betaHeatAdj", dimless, heatAdjBeta);
const dimensionedScalar adjPseudoUb("adjPseudoUb", dimless/dimTime, adjPseudoInvDtUb);
const volScalarField magUStageUb("magUStageUb", mag(U));
const scalar UavgStageUb = Foam::max(gAverage(magUStageUb.primitiveField()), 1e-8);
const scalar UmaxStageUb = Foam::max(gMax(magUStageUb.primitiveField()), UavgStageUb);
const volScalarField magUbSafe("magUbSafe", mag(UbSafe));
const volScalarField magPbSafe("magPbSafe", mag(pbSafe));
const scalar UbSafeAvg = Foam::max(gAverage(magUbSafe.primitiveField()), 1e-8);
const scalar pbSafeAvg = Foam::max(gAverage(magPbSafe.primitiveField()), 1e-8);
const scalar UbHardCap = Foam::min(Foam::max(50.0*UbSafeAvg, 1e3), Foam::max(1e6*UmaxStageUb, 1e3));
const scalar pbHardCap = Foam::min(Foam::max(50.0*pbSafeAvg, 1e4), Foam::max(1e6*UmaxStageUb*UmaxStageUb, 1e5));
volVectorField adjointTransposeConvectionb((fvc::grad(U) & Ub));

tmp<fvVectorMatrix> tUbEqn
(
    fvm::div(-phi, Ub)
  + adjointTransposeConvectionb
  - fvm::laplacian(turbulence->nuEff(), Ub)
  + fvm::Sp(adjPseudoUb, Ub)
  + fvm::Sp(alpha, Ub)
  + betaHeatAdj*b1*Tb*fvc::grad(T)
 ==
    fvOptions(Ub)
);
fvVectorMatrix& UbEqn = tUbEqn.ref();

UbEqn.relax();

fvOptions.constrain(UbEqn);

SolverPerformance<vector> spUb = solve(UbEqn == -fvc::grad(pb));
logSolverStats(runTime, "AdjointHeat_Ub", i, "Ub", spUb, opt, runTime.timeName());

fvOptions.correct(Ub);
forAll(Ub, cellI)
{
    vector& uv = Ub[cellI];
    if (!(std::isfinite(uv.x()) && std::isfinite(uv.y()) && std::isfinite(uv.z())))
    {
        uv = vector::zero;
    }
    else
    {
        const scalar um = mag(uv);
        if (um > UbHardCap)
        {
            uv *= (UbHardCap/(um + VSMALL));
        }
    }
}
//****************************************
volScalarField rAUb(1.0/UbEqn.A());
const scalar rAUbAvg = Foam::max(gAverage(rAUb.primitiveField()), 1e-12);
const scalar rAUbMin = 1e-6*rAUbAvg;
const scalar rAUbMax = 1e6*rAUbAvg;
forAll(rAUb, cellI)
{
    scalar& rv = rAUb[cellI];
    if (!std::isfinite(rv))
    {
        rv = rAUbMin;
    }
    else
    {
        rv = Foam::min(Foam::max(rv, rAUbMin), rAUbMax);
    }
}
MTOdiag::logScalarFieldHealth(runTime, opt, ubStage, "rAUb", rAUb.primitiveField());
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
        if (!std::isfinite(rv))
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
MTOdiag::logScalarFieldHealth(runTime, opt, ubStage, "rAtUb", rAtUb().primitiveField());
const surfaceScalarField magPhiBase("magPhiBase", mag(phi));
const scalar phiBaseAvg = Foam::max(gAverage(magPhiBase.primitiveField()), 1e-12);
const scalar phiCapPhys = Foam::max(1e3*phiBaseAvg, 1e-8);
const scalar phiCapHist = (prevPhiHbyAbMax > 0.0)
    ? Foam::max(5.0*prevPhiHbyAbMax, 1e-8)
    : phiCapPhys;
const scalar phiAbMaxAbs = Foam::min(phiCapPhys, phiCapHist);
scalar phiUsedMax = 0.0;
forAll(phiHbyAb, faceI)
{
    scalar& phiv = phiHbyAb[faceI];
    if (!std::isfinite(phiv))
    {
        phiv = 0.0;
    }
    else
    {
        phiv = Foam::min(Foam::max(phiv, -phiAbMaxAbs), phiAbMaxAbs);
    }
    phiUsedMax = Foam::max(phiUsedMax, Foam::mag(phiv));
}
prevPhiHbyAbMax = phiUsedMax;
MTOdiag::logScalarFieldHealth(runTime, opt, ubStage, "phiHbyAb", phiHbyAb.primitiveField());
MTOdiag::logMetric(runTime, opt, ubStage, "phiHbyAb_capPhys", phiCapPhys);
MTOdiag::logMetric(runTime, opt, ubStage, "phiHbyAb_capHist", phiCapHist);
MTOdiag::logMetric(runTime, opt, ubStage, "phiHbyAb_capUsed", phiAbMaxAbs);
MTOdiag::logMetric(runTime, opt, ubStage, "phiHbyAb_usedMax", phiUsedMax);
MTOdiag::logMetric(runTime, opt, ubStage, "heatAdjBeta", heatAdjBeta);
tUbEqn.clear();

// Update the pressure BCs to ensure flux consistency
constrainPressure(pb, Ub, phiHbyAb, rAtUb(), MRF);  
// Non-orthogonal pressure corrector loop
bool badPbResidual = false;
scalar pbFirstInitRes = 1e300;
scalar pbLastFinalRes = 1e300;
label nonOrthSolveI = 0;
while (simple.correctNonOrthogonal())
{
    fvScalarMatrix pbEqn
    (
        fvm::laplacian(rAtUb(), pb) == fvc::div(phiHbyAb)
    );

    pbEqn.setReference(pbRefCell, pbRefValue);
    SolverPerformance<scalar> sppb = pbEqn.solve();
    logSolverStats(runTime, "AdjointHeat_Ub", i, "pb", sppb, opt, runTime.timeName());
    const scalar pbInitRes = sppb.initialResidual();
    pbLastFinalRes = sppb.finalResidual();
    if (nonOrthSolveI == 0)
    {
        pbFirstInitRes = pbInitRes;
    }
    nonOrthSolveI++;
    if (!std::isfinite(pbInitRes) || !std::isfinite(pbLastFinalRes))
    {
        badPbResidual = true;
        break;
    }
    MTOdiag::logMetric(runTime, opt, ubStage, "pb_initRes", pbInitRes);
    MTOdiag::logMetric(runTime, opt, ubStage, "pb_finalRes", pbLastFinalRes);

    if (simple.finalNonOrthogonalIter())
    {
        phib = phiHbyAb - pbEqn.flux();
    }
}

if (badPbResidual)
{
    Info<< "Warning: non-finite pb residual detected; restoring last safe Ub/pb state." << endl;
    Ub = UbSafe;
    pb = pbSafe;
    phib = phibSafe;
    break;
}

if (i > 0 && pbFirstInitRes > 1.15*prevPbInitRes)
{
    pbDivergeCount++;
}
else
{
    pbDivergeCount = 0;
}
prevPbInitRes = pbFirstInitRes;
MTOdiag::logMetric(runTime, opt, ubStage, "pb_firstInitRes", pbFirstInitRes);
MTOdiag::logMetric(runTime, opt, ubStage, "pb_driftCount", scalar(pbDivergeCount));
if (pbDivergeCount >= 2)
{
    Info<< "Warning: pb first-pass residual drifting up; restoring last safe Ub/pb state and exiting inner loop." << endl;
    Ub = UbSafe;
    pb = pbSafe;
    phib = phibSafe;
    break;
}

// Guard against non-finite adjoint fields before correction
forAll(Ub, cellI)
{
    vector& uv = Ub[cellI];
    if (!(std::isfinite(uv.x()) && std::isfinite(uv.y()) && std::isfinite(uv.z())))
    {
        uv = vector::zero;
    }
}
forAll(pb, cellI)
{
    scalar& pv = pb[cellI];
    if (!std::isfinite(pv))
    {
        pv = 0.0;
    }
    else if (pv > pbHardCap)
    {
        pv = pbHardCap;
    }
    else if (pv < -pbHardCap)
    {
        pv = -pbHardCap;
    }
}
MTOdiag::logVectorFieldHealth(runTime, opt, ubStage, "Ub", Ub.primitiveField());
MTOdiag::logScalarFieldHealth(runTime, opt, ubStage, "pb", pb.primitiveField());
const volScalarField magUbNow("magUbNow", mag(Ub));
const volScalarField magPbNow("magPbNow", mag(pb));
const scalar UbMaxNow = Foam::max(magUbNow.primitiveField());
const scalar pbAbsMaxNow = Foam::max(magPbNow.primitiveField());
MTOdiag::logMetric(runTime, opt, ubStage, "Ub_hardCap", UbHardCap);
MTOdiag::logMetric(runTime, opt, ubStage, "pb_hardCap", pbHardCap);
MTOdiag::logMetric(runTime, opt, ubStage, "Ub_absMax", UbMaxNow);
MTOdiag::logMetric(runTime, opt, ubStage, "pb_absMax", pbAbsMaxNow);
if (!std::isfinite(UbMaxNow) || !std::isfinite(pbAbsMaxNow) || UbMaxNow > 1.05*UbHardCap || pbAbsMaxNow > 1.05*pbHardCap)
{
    Info<< "Warning: Ub/pb exceeded hard cap; restoring last safe Ub/pb state and exiting inner loop." << endl;
    Ub = UbSafe;
    pb = pbSafe;
    phib = phibSafe;
    break;
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
UbSafe = Ub;
pbSafe = pb;
phibSafe = phib;

if (pbLastFinalRes < 5e-5 && i >= 2)
{
    break;
}
}