//***********************************// LEGACY: Single-stage PDE filter (set aside for cascaded)
// Preserved from original implementation. Use filter_x.H (cascaded) as default.
//***********************************//PDE filter
{
    fvScalarMatrix xpEqn(fvm::laplacian(xp)-fvm::Sp(b,xp)+x*b);
    SolverPerformance<scalar> spXp = xpEqn.solve();
    logSolverStats(runTime, "filter_x", 0, "xp", spXp, opt, runTime.timeName());
}
//***********************************//Heaviside_filter
// Turbulent tuning: gentler ramp, lower cap, delayed sharpening (heavisideDelFactor/Cap/Start/Init)
del = (opt > heavisideDelStart)
    ? Foam::min(heavisideDelFactor*(opt - heavisideDelStart), heavisideDelCap)
    : heavisideDelInit;
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
