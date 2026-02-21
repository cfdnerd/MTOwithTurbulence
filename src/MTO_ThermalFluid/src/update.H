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
// Design-boundary blend: blend design cells adjacent to fluid_area toward xh=1 (smooth interface)
if(fluid_area && designBoundaryBlend > 1e-12)
{
   const labelList& own = mesh.owner();
   const labelList& nei = mesh.neighbour();
   HashSet<label> fluidSet(cells_fluid.size());
   forAll(cells_fluid, i) fluidSet.insert(cells_fluid[i]);
   HashSet<label> solidSet;
   if(solid_area) { solidSet.resize(cells_solid.size()); forAll(cells_solid, i) solidSet.insert(cells_solid[i]); }
   HashSet<label> testSet;
   if(test_area)  { testSet.resize(cells_test.size());  forAll(cells_test, i)  testSet.insert(cells_test[i]);  }
   forAll(own, facei)
   {
      label o = own[facei];
      label n = nei[facei];
      bool oFluid = fluidSet.found(o);
      bool nFluid = fluidSet.found(n);
      bool nDesign = !fluidSet.found(n) && (!solid_area || !solidSet.found(n)) && (!test_area || !testSet.found(n));
      bool oDesign = !fluidSet.found(o) && (!solid_area || !solidSet.found(o)) && (!test_area || !testSet.found(o));
      if(oFluid && nDesign)
      {
         scalar& xhn = xh[n];
         xhn = (1.0-designBoundaryBlend)*xhn + designBoundaryBlend*1.0;
      }
      else if(!oFluid && nFluid && oDesign)
      {
         scalar& xho = xh[o];
         xho = (1.0-designBoundaryBlend)*xho + designBoundaryBlend*1.0;
      }
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
if(opt>=quRampStart)
{
  qu=Foam::min(quInit+(opt-quRampStart)*quRampRate, quCap);
}
else
{
  qu=quInit;
}
// Adaptive xh floor: when design is solid-heavy (min xh low), raise floor to avoid alpha spikes and FPE
// xhFloorTurbMult raises floor for turbulent flow (reduces Brinkman stiffness)
scalar xhFloorBase = (xhMinVal < 0.15) ? 0.03 : ((xhMinVal < 0.05) ? 0.02 : 1e-3);
scalar xhFloorVal = Foam::min(xhFloorBase * xhFloorTurbMult, 0.05);
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
