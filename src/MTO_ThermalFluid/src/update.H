//***********************************//update mu and lambda
opt++;
nuEff = turbulence->nuEff();
if(solid_area)
{
   forAll(cells_solid, i)
   {
      xh[cells_solid[i]]=0;
   }
}
if(fluid_area)
{
   forAll(cells_fluid, i)
   {
      xh[cells_fluid[i]]=1.0;
   }
}
if(test_area)
{
   forAll(cells_test, i)
   {
      xh[cells_test[i]]=1.0;
   }
}
xh.correctBoundaryConditions();
xh = Foam::max(0.0, Foam::min(1.0, xh));
xh.correctBoundaryConditions();
if(runTime.writeTime())
{
   xh.write();  
}

if(opt<=63)
{
  alphaMax=alphamax*(opt/7.0+1.0);
}
// Freeze alphaMax when design is extreme (solid-heavy) to avoid ill-conditioning
scalar xhMinVal = gMin(xh);
if(opt>100 && xhMinVal >= 0.15)
{
  alphaMax=alphaMax*1.02;
}
alphaMax.value()=Foam::min(alphaMax.value(),alphaMAX.value());
// Mesh-scaled cap: alpha ∝ nuEff/h^2 keeps pressure matrix condition number bounded (Brinkman lit.)
// Use nuEff (not nu) for turbulent flow consistency
scalar alphaMaxMesh = 1e30;
{
  scalar nCell = returnReduce(mesh.nCells(), sumOp<label>());
  scalar volTot = gSum(mesh.V());
  scalar hChar = Foam::pow(max(volTot/nCell, 1e-30), 1.0/3.0);
  scalar nuEffVal = gAverage(nuEff.primitiveField());
  alphaMaxMesh = alphaMaxMeshFactor*max(nuEffVal, nu.value())/max(hChar*hChar, 1e-30);
  alphaMax.value() = Foam::min(alphaMax.value(), alphaMaxMesh);
}
if(opt>=80)
{
  qu=0.005+(opt-80)*1e-4;
}
qu=Foam::min(qu,0.01);
// Adaptive xh floor: when design is solid-heavy (min xh low), raise floor to avoid alpha spikes and FPE
scalar xhFloorVal = (xhMinVal < 0.15) ? 0.03 : ((xhMinVal < 0.05) ? 0.02 : 1e-3);
const volScalarField xhSafe(Foam::max(xh, scalar(xhFloorVal)));
volScalarField alphaRaw(alphaMax*qu*(1-xhSafe)/(qu+xhSafe));
scalar alphaCapVal = (xhMinVal < 0.15) ? Foam::min(scalar(3e3), alphaMaxMesh) : Foam::min(scalar(1e4), alphaMaxMesh);
alpha=Foam::min(alphaRaw, dimensionedScalar("alphaCap", alphaRaw.dimensions(), alphaCapVal));
// Use xhSafe in DT to limit thermal diffusivity range and improve matrix conditioning
DT=(ks+(kf-ks)*xhSafe*(1+qu)/(qu+xhSafe))/rhoc;

// Diagnostics checkpoint: track coefficient health right after update
MTOdiag::logScalarFieldHealth(runTime, opt, "update", "xh", xh.primitiveField());
MTOdiag::logScalarFieldHealth(runTime, opt, "update", "alpha", alpha.primitiveField());
MTOdiag::logScalarFieldHealth(runTime, opt, "update", "DT", DT.primitiveField());
MTOdiag::logMetric(runTime, opt, "update", "alphaMax_value", alphaMax.value());
MTOdiag::logMetric(runTime, opt, "update", "alphaCap", alphaCapVal);
MTOdiag::logMetric(runTime, opt, "update", "xhFloor", xhFloorVal);
