//Tb=T2;
sumT=0;
maxT=gMax(T);
for(i=0;i<n;i++)
{
   sumT=sumT+std::pow(T[i]/maxT,Pnorm);
}      
reduce(sumT, sumOp<scalar>());
sumT=std::pow(sumT,1.0/Pnorm-1);

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
   TEqnb -= b2*sumT*Foam::pow(T/maxT,Pnorm-1);
}
SolverPerformance<scalar> spTb = TEqnb.solve();
logSolverStats(runTime, "AdjointHeat_Tb", 0, "Tb", spTb, opt, runTime.timeName());
