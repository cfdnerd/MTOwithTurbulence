// Tb equation robustness: avoid unconditional p-norm normalization and guard solve state.
const word tbStage("AdjointHeat_Tb");
MTOdiag::logMetric(runTime, opt, tbStage, "tb_phase_enter", 1.0);
volScalarField TbSafe(Tb);
// Reset Tb before solve if excessive (avoids ill-conditioned matrix and FPE in PBiCGStab)
{
    const volScalarField magTbIn("magTbIn", mag(Tb));
    const scalar tbAbsMaxIn = gMax(magTbIn.primitiveField());
    if (!std::isfinite(tbAbsMaxIn) || tbAbsMaxIn > tbPreSolveReset)
    {
        Tb = dimensionedScalar("zero", Tb.dimensions(), 0.0);
        TbSafe = Tb;
        if (std::isfinite(tbAbsMaxIn) && tbAbsMaxIn > tbPreSolveReset)
        {
            Info<< "AdjointHeat_Tb: reset Tb to zero (max|Tb|=" << tbAbsMaxIn << " > " << tbPreSolveReset << ")" << endl;
        }
        MTOdiag::logMetric(runTime, opt, tbStage, "tb_preSolve_reset", 1.0);
    }
}
sumT = 0;
const volScalarField magTNow("magTNow", mag(T));
scalar maxTAbs = Foam::max(gMax(magTNow.primitiveField()), 1e-12);
if (!std::isfinite(maxTAbs))
{
    maxTAbs = 1e-12;
    MTOdiag::logMetric(runTime, opt, tbStage, "tb_maxTAbs_reset", 1.0);
}
MTOdiag::logMetric(runTime, opt, tbStage, "tb_maxTAbs", maxTAbs);
MTOdiag::logScalarFieldHealth(runTime, opt, tbStage, "T", T.primitiveField());
if (objFunction == 2)
{
    for (i=0; i<n; i++)
    {
        const scalar tNorm = Foam::max(Foam::mag(T[i])/maxTAbs, 1e-12);
        sumT += std::pow(tNorm, Pnorm);
    }
    reduce(sumT, sumOp<scalar>());
    sumT = std::pow(Foam::max(sumT, scalar(1e-30)), 1.0/Pnorm - 1.0);
    if (!std::isfinite(sumT))
    {
        sumT = 0.0;
        MTOdiag::logMetric(runTime, opt, tbStage, "tb_sumT_reset", 1.0);
    }
}
MTOdiag::logMetric(runTime, opt, tbStage, "tb_sumT", sumT);
MTOdiag::logMetric(runTime, opt, tbStage, "tb_phase_after_norm", 1.0);

// alphaEff floored in Primal_T; use same field for consistency
const dimensionedScalar adjPseudoTb("adjPseudoTb", dimless/dimTime, adjPseudoInvDtTb);
MTOdiag::logScalarFieldHealth(runTime, opt, tbStage, "alphaEff", alphaEff.primitiveField());
MTOdiag::logMetric(runTime, opt, tbStage, "tb_phase_pre_assembly", 1.0);
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
// Cap RHS source to avoid FPE in PBiCGStab from extreme residual values (crash-fix opt313)
{
    scalarField& src = TEqnb.source();
    forAll(src, i)
    {
        if (!std::isfinite(src[i]))
            src[i] = 0.0;
    }
    scalar srcMax = gMax(mag(src));
    if (!std::isfinite(srcMax))
        srcMax = 0.0;
    if (srcMax > tbRhsCapMax && srcMax > 1e-30)
    {
        src *= tbRhsCapMax / srcMax;
        MTOdiag::logMetric(runTime, opt, tbStage, "tb_rhs_capped", 1.0);
    }
}
MTOdiag::logMetric(runTime, opt, tbStage, "tb_phase_pre_solve", 1.0);
const bool primalTSingular =
(
    spT.singular()
 || !std::isfinite(spT.initialResidual())
 || !std::isfinite(spT.finalResidual())
 || spT.finalResidual() > 1.0
);
SolverPerformance<scalar> spTb;
if (primalTSingular)
{
    Tb = TbSafe;
    Info<< "AdjointHeat_Tb: skipping solve (primal T singular or bad residual); keeping previous Tb." << endl;
    MTOdiag::logMetric(runTime, opt, tbStage, "tb_skip_primalT_bad", 1.0);
}
else
{
    spTb = TEqnb.solve();
    logSolverStats(runTime, "AdjointHeat_Tb", 0, "Tb", spTb, opt, runTime.timeName());
}
const scalar tbInitRes = primalTSingular ? 0.0 : spTb.initialResidual();
const scalar tbFinalRes = primalTSingular ? 0.0 : spTb.finalResidual();
MTOdiag::logMetric(runTime, opt, tbStage, "tb_initRes", tbInitRes);
MTOdiag::logMetric(runTime, opt, tbStage, "tb_finalRes", tbFinalRes);
MTOdiag::logMetric(runTime, opt, tbStage, "tb_phase_post_solve", 1.0);
if (!primalTSingular && (!std::isfinite(tbInitRes) || !std::isfinite(tbFinalRes)))
{
    Info<< "Warning: non-finite Tb residual; restoring last safe Tb state." << endl;
    Tb = TbSafe;
}
const volScalarField magTbSafe("magTbSafe", mag(TbSafe));
const scalar tbSafeAvg = Foam::max(gAverage(magTbSafe.primitiveField()), 1e-8);
const scalar tbHardCap = Foam::min(Foam::max(10.0*tbSafeAvg, tbHardCapMin), tbHardCapMax);
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
MTOdiag::logMetric(runTime, opt, tbStage, "tb_phase_post_clip", 1.0);
MTOdiag::logScalarFieldHealth(runTime, opt, tbStage, "Tb", Tb.primitiveField());
const volScalarField magTbNow("magTbNow", mag(Tb));
const scalar tbAbsMax = gMax(magTbNow.primitiveField());
MTOdiag::logMetric(runTime, opt, tbStage, "tb_hardCap", tbHardCap);
MTOdiag::logMetric(runTime, opt, tbStage, "tb_absMax", tbAbsMax);
MTOdiag::logMetric(runTime, opt, tbStage, "tb_phase_exit", 1.0);
