Info<< "Initialization\n" << endl;
if (Pstream::master())
{
    const fileName caseRoot(runTime.rootPath()/runTime.globalCaseName());
    const fileName histFile(caseRoot/"optimization.hst");
    std::remove(histFile.c_str());
    const fileName solverLogFile(caseRoot/"solver_status.log");
    std::remove(solverLogFile.c_str());
} 

int i=0,opt=0;
double del=0;
scalar n(0),N(0),eta0,eta5,eta1,y0,y5;
n=x.size();//n:length of x in the i-th core
N=x.size();//N:design freedom
reduce(N, sumOp<scalar>()); 

scalar V(0),area=gSum(mesh.V()),MeanT(0),PowerDiss(0),maxT(0),sumT(0);
double len = filterLenAniso
    ? std::pow(Foam::max(area/N, 1e-30), 1.0/Foam::max(GeoDim, 1))
    : std::pow(Foam::max(area/N, 1e-30), 1.0/3);
if(filterTurbReScale)
{
    scalar Uref = Foam::max(gAverage(mag(U.primitiveField())), 1e-12);
    scalar Re = Foam::max(Uref*len/nu.value(), 1e-6);
    scalar turbScale = Foam::max(Foam::min(Foam::sqrt(Re_ref/Re), 2.0), 0.5);
    len *= turbScale;
}
double denomSq1 = Foam::max(filterR*len/3.464, 1e-30);
double denomSq2 = Foam::max(filterR2*len/3.464, 1e-30);
b.value()=1.0/(denomSq1*denomSq1);
bFilter2.value()=1.0/(denomSq2*denomSq2);

// Write mesh/filter estimates at top of solver_status.log (master only)
writeMeshFilterHeader(runTime, mesh, area, N, len, filterR, filterR2, GeoDim, nu, U);

//***********************************//MMA
int m=2;//number of constraint
MMA mma(n,m);
mma.asyminit = mma_init;
mma.asymdec = mma_dec;
mma.asyminc = mma_inc;
mma.raa0 = raa0;
mma.movelimit = movlim;
std::vector<double> xmma(n);
std::vector<double> dfdx(n);
std::vector<double> gx(m);
std::vector<std::vector<double>> dgdx(m);
for (int i = 0; i < m; i++)
{
  dgdx[i].resize(n);
}

