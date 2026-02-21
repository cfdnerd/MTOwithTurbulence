//***********************************// Cascaded PDE filter (Lazarov & Wang 2016)
// Two-stage Helmholtz filter for gray suppression in turbulent channel flows.
// Stage 1: x -> xp1  |  Stage 2: xp1 -> xp  |  Heaviside: xp -> xh
//***********************************//Stage 1: PDE filter x -> xp1
{
    fvScalarMatrix xp1Eqn(fvm::laplacian(xp1)-fvm::Sp(b,xp1)+x*b);
    SolverPerformance<scalar> spXp1 = xp1Eqn.solve();
    logSolverStats(runTime, "filter_x", 0, "xp1", spXp1, opt, runTime.timeName());
}
//***********************************//Stage 2: PDE filter xp1 -> xp
{
    fvScalarMatrix xpEqn(fvm::laplacian(xp)-fvm::Sp(bFilter2,xp)+xp1*bFilter2);
    SolverPerformance<scalar> spXp = xpEqn.solve();
    logSolverStats(runTime, "filter_x", 1, "xp", spXp, opt, runTime.timeName());
}
//***********************************//Heaviside projection (sharper ramp for cascaded)
del=Foam::min(0.4*opt, 200);
eta0=0.0001;
eta1=1.0;
y0=diff(xp,mesh.V(),del,eta0,n);
reduce(y0, sumOp<scalar>());
do
{
  eta5=(eta0+eta1)/2.0;
  y5=diff(xp,mesh.V(),del,eta5,n);
  reduce(y5, sumOp<scalar>());
  if(y0*y5<0)
  { 
     eta1=eta5;
  }
  else
  { 
     eta0=eta5;
     y0=y5;
  }
} while ((eta1-eta0)>0.0001);

for(i=0;i<n;i++)
{
   if(xp[i]<=eta5)
   {
     xh[i]=eta5*(Foam::exp(-del*(1-xp[i]/eta5))-(1-xp[i]/eta5)*Foam::exp(-del)); 
   }
   else
   {
     xh[i]=eta5+(1-eta5)*(1-Foam::exp(-del*(xp[i]-eta5)/(1-eta5))+(xp[i]-eta5)*Foam::exp(-del)/(1-eta5));
   }   
}
