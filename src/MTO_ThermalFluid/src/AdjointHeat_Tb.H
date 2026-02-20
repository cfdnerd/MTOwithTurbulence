// Tb equation robustness: avoid unconditional p-norm normalization and guard solve state.
volScalarField TbSafe(Tb);
sumT = 0;
const scalar maxTAbs = Foam::max(gMax(mag(T).primitiveField()), 1e-12);
if (objFunction == 2)
{
    for (i=0; i<n; i++)
    {
        const scalar tNorm = Foam::max(Foam::mag(T[i])/maxTAbs, 1e-12);
        sumT += std::pow(tNorm, Pnorm);
    }
    reduce(sumT, sumOp<scalar>());
    sumT = std::pow(Foam::max(sumT, scalar(1e-30)), 1.0/Pnorm - 1.0);
}

// alphaEff floored in Primal_T; use same field for consistency
const dimensionedScalar adjPseudoTb("adjPseudoTb", dimless/dimTime, adjPseudoInvDtTb);
fvScalarMatrix TEqnb
(
   -fvm::div(-phi,Tb) 
   +fvm::laplacian(alphaEff, Tb)
   +fvm::Sp(adjPseudoTb, Tb)
);
if (objFunction==1)
{
   TEqnb -= Q;
}
if (objFunction==2)
{
   const volScalarField tNormSafe
   (
       "tNormSafe",
       Foam::max(mag(T)/maxTAbs, scalar(1e-12))
   );
   TEqnb -= b2*sumT*Foam::pow(tNormSafe, Pnorm-1);
}
SolverPerformance<scalar> spTb = TEqnb.solve();
logSolverStats(runTime, "AdjointHeat_Tb", 0, "Tb", spTb, opt, runTime.timeName());
const scalar tbInitRes = spTb.initialResidual();
const scalar tbFinalRes = spTb.finalResidual();
MTOdiag::logMetric(runTime, opt, "AdjointHeat_Tb", "tb_initRes", tbInitRes);
MTOdiag::logMetric(runTime, opt, "AdjointHeat_Tb", "tb_finalRes", tbFinalRes);
if (!std::isfinite(tbInitRes) || !std::isfinite(tbFinalRes))
{
    Info<< "Warning: non-finite Tb residual; restoring last safe Tb state." << endl;
    Tb = TbSafe;
}
const volScalarField magTbSafe("magTbSafe", mag(TbSafe));
const scalar tbSafeAvg = Foam::max(gAverage(magTbSafe.primitiveField()), 1e-8);
const scalar tbHardCap = Foam::max(50.0*tbSafeAvg, 1e3);
forAll(Tb, cellI)
{
    scalar& tv = Tb[cellI];
    if (!std::isfinite(tv))
    {
        tv = 0.0;
    }
    else if (tv > tbHardCap)
    {
        tv = tbHardCap;
    }
    else if (tv < -tbHardCap)
    {
        tv = -tbHardCap;
    }
}
MTOdiag::logScalarFieldHealth(runTime, opt, "AdjointHeat_Tb", "Tb", Tb.primitiveField());
const volScalarField magTbNow("magTbNow", mag(Tb));
const scalar tbAbsMax = gMax(magTbNow.primitiveField());
MTOdiag::logMetric(runTime, opt, "AdjointHeat_Tb", "tb_hardCap", tbHardCap);
MTOdiag::logMetric(runTime, opt, "AdjointHeat_Tb", "tb_absMax", tbAbsMax);
