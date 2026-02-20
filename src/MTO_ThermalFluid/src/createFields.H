Info<< "Reading optimization parameters\n" << endl;
IOdictionary optProperties
(
   IOobject
   (
       "optProperties",
       runTime.constant(),
       mesh,
       IOobject::MUST_READ_IF_MODIFIED,
       IOobject::NO_WRITE
   )
);
//optimization parameters
scalar raa0(readScalar(optProperties.lookup("raa0")));
scalar mma_init(readScalar(optProperties.lookup("mma_init")));
scalar mma_dec(readScalar(optProperties.lookup("mma_dec")));
scalar mma_inc(readScalar(optProperties.lookup("mma_inc")));
scalar movlim(readScalar(optProperties.lookup("movlim")));
scalar voluse(readScalar(optProperties.lookup("voluse")));
scalar GeoDim(readScalar(optProperties.lookup("GeoDim")));
scalar PowerDiss0(readScalar(optProperties.lookup("PowerDiss0")));
scalar PowerDissMax(readScalar(optProperties.lookup("PowerDissMax")));
scalar PowerDissRelax(readScalar(optProperties.lookup("PowerDissRelax")));
scalar filterR(readScalar(optProperties.lookup("filterR")));
dimensionedScalar b("b", dimensionSet(0,-2,0,0,0,0,0),1.0);//filter radius of the PDE filter
dimensionedScalar b1("b1", dimensionSet(0,2,-2,-2,0,0,0),1.0);
Switch solid_area(optProperties.lookup("solid_area"));
Switch fluid_area(optProperties.lookup("fluid_area"));
Switch test_area(optProperties.lookup("test_area"));
scalar Pnorm(readScalar(optProperties.lookup("Pnorm")));
scalar objFunction(readScalar(optProperties.lookup("objFunction")));
scalar nAdjTurbPasses(optProperties.lookupOrDefault<scalar>("nAdjTurbPasses", 1));
scalar turbSensDamp(optProperties.lookupOrDefault<scalar>("turbSensDamp", 0));
// Robust adjoint framework controls
scalar adjPseudoInvDtUa(optProperties.lookupOrDefault<scalar>("adjPseudoInvDtUa", 50.0));
scalar adjPseudoInvDtUb(optProperties.lookupOrDefault<scalar>("adjPseudoInvDtUb", 25.0));
scalar adjPseudoInvDtTb(optProperties.lookupOrDefault<scalar>("adjPseudoInvDtTb", 10.0));
scalar adjPseudoInvDtTurb(optProperties.lookupOrDefault<scalar>("adjPseudoInvDtTurb", 100.0));
scalar adjTurbBetaMin(optProperties.lookupOrDefault<scalar>("adjTurbBetaMin", 0.2));
label adjTurbRampStart(optProperties.lookupOrDefault<label>("adjTurbRampStart", 1));
label adjTurbRampEnd(optProperties.lookupOrDefault<label>("adjTurbRampEnd", 80));
label adjTurbLagEvery(optProperties.lookupOrDefault<label>("adjTurbLagEvery", 2));
scalar adjHeatBetaMin(optProperties.lookupOrDefault<scalar>("adjHeatBetaMin", 0.2));
label adjHeatRampStart(optProperties.lookupOrDefault<label>("adjHeatRampStart", 1));
label adjHeatRampEnd(optProperties.lookupOrDefault<label>("adjHeatRampEnd", 120));
scalar adjUaPassDecay(optProperties.lookupOrDefault<scalar>("adjUaPassDecay", 0.35));
scalar adjUaPassMinGain(optProperties.lookupOrDefault<scalar>("adjUaPassMinGain", 0.10));
scalar adjUaLocalCapCoeff(optProperties.lookupOrDefault<scalar>("adjUaLocalCapCoeff", 200.0));
scalar adjUaRhsCapCoeff(optProperties.lookupOrDefault<scalar>("adjUaRhsCapCoeff", 200.0));
scalar adjUaAllPassTbMax(optProperties.lookupOrDefault<scalar>("adjUaAllPassTbMax", 1e74));
scalar adjUaHigherPassTbMax(optProperties.lookupOrDefault<scalar>("adjUaHigherPassTbMax", 1e60));
scalar tbPreSolveReset(optProperties.lookupOrDefault<scalar>("tbPreSolveReset", 1e6));
scalar tbHardCapMax(optProperties.lookupOrDefault<scalar>("tbHardCapMax", 1e6));
scalar tbRhsCapMax(optProperties.lookupOrDefault<scalar>("tbRhsCapMax", 1e12));
label diagLevel(optProperties.lookupOrDefault<label>("diagLevel", 1));
label diagEvery(optProperties.lookupOrDefault<label>("diagEvery", 5));
label diagRank(optProperties.lookupOrDefault<label>("diagRank", -1));
label diagOptStart(optProperties.lookupOrDefault<label>("diagOptStart", 120));
Info<< "Reading pseudo-density field\n" << endl;
volScalarField xp
(
    IOobject
    (
        "xp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    scalar(voluse),
    zeroGradientFvPatchScalarField::typeName
);
volScalarField xh
(
    IOobject
    (
        "xh",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    scalar(voluse),
    zeroGradientFvPatchScalarField::typeName
);
volScalarField x(xp);
volScalarField drho(x);

// find the cells in the non-design area
labelList cells_solid,cells_fluid,cells_test;
if(solid_area)
{
  word zoneName="solid_area";
  label zoneID=mesh.cellZones().findZoneID(zoneName);
  cells_solid=mesh.cellZones()[zoneID];
  forAll(cells_solid, i)
  {
     x[cells_solid[i]]=0;
  }
}
if(fluid_area)
{
  word zoneName="fluid_area";
  label zoneID=mesh.cellZones().findZoneID(zoneName);
  cells_fluid=mesh.cellZones()[zoneID];
  forAll(cells_fluid, i)
  {
     x[cells_fluid[i]]=1.0;
  }
}
if(test_area)
{
  word zoneName="test_area";
  label zoneID=mesh.cellZones().findZoneID(zoneName);
  cells_test=mesh.cellZones()[zoneID];
  forAll(cells_test, i)
  {
     x[cells_test[i]]=1.0;
     xh[cells_test[i]]=1.0;
  }
}
volScalarField gsensVol
(
    IOobject
    (
        "gsensVol",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    xh,
    zeroGradientFvPatchScalarField::typeName
);
volScalarField gsenshVol(xh);

// find the cells on the boundary
dictionary conFunc = mesh.solutionDict().subDict("constraintFunctionDict_flow");
int nObjPatch =conFunc.lookupOrDefault<scalar>("numberConstraintPatches",2);
wordList conPatchNames = conFunc.lookup("constraintPatchesNames");
label conPatchList[nObjPatch];//patch list
int iLoop;
for(iLoop=0; iLoop<nObjPatch; iLoop++)
{
    conPatchList[iLoop]=mesh.boundaryMesh().findPatchID(conPatchNames[iLoop]);
}   
