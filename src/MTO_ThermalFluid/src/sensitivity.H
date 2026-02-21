Info<<"sensitivity analysis"<<endl;

fsenshMeanT=-alphaMax*(1+qu)*qu/((qu+xh)*(qu+xh))*(U&Ub)+b1*(kf-ks)*(1+qu)*qu/((qu+xh)*(qu+xh))*(fvc::grad(T)&fvc::grad(Tb))/rhoc;
if(DDHS)//design dependent heat source
{
 fsenshMeanT += b1*Tb*heatSource/rhoc;
}
gsenshPowerDiss=-alphaMax*(1+qu)*qu/((qu+xh)*(qu+xh))*(U&Ua);

// Turbulence sensitivity: ka*dRk/dxh + omegaa*dRomega/dxh when turbSensDamp>0 (Rk~-c*alpha*k, Romega~-c*alpha*omega in solid)
if (turbSensDamp > 0 && mesh.thisDb().foundObject<volScalarField>("k") && mesh.thisDb().foundObject<volScalarField>("omega"))
{
    const volScalarField& ka = mesh.thisDb().lookupObject<volScalarField>("ka");
    const volScalarField& omegaa = mesh.thisDb().lookupObject<volScalarField>("omegaa");
    const volScalarField& k = mesh.thisDb().lookupObject<volScalarField>("k");
    const volScalarField& omega = mesh.thisDb().lookupObject<volScalarField>("omega");
    const volScalarField dAlphaDxh(alphaMax*(1+qu)*qu/((qu+xh)*(qu+xh)));
    gsenshPowerDiss -= turbSensDamp * (ka * k + omegaa * omega) * dAlphaDxh;
    if (objFunction == 1)
    {
        fsenshMeanT -= turbSensDamp * (ka * k + omegaa * omega) * dAlphaDxh;
    }
}

if(useCascadedFilter)
{
    #include "filter_chainrule.H"
}
else
{
    #include "filter_chainrule_legacy.H"
}

gx[0]=(PowerDiss/Foam::max(PowerDissRelax-opt*0.2,PowerDissMax)-1.0);
gx[1]=V;

for(i=0;i<n;i++)
{
 xmma[i] = x[i];
 dfdx[i]=fsensMeanT[i]/N;//sensitivities of objective function
 dgdx[0][i]=gsensPowerDiss[i]/N;//sensitivities of constraint function
 dgdx[1][i]=gsensVol[i]/N;//sensitivities of constraint function
}


Info<< "\n run MMA \n" << endl;
mma.MMAsolver(xmma, dfdx, gx, dgdx);
// Backtracked MMA update: accept only finite, bounded design states
std::vector<double> xOld(n);
for(i=0;i<n;i++)
{
    xOld[i] = x[i];
}
scalar lambda = 1.0;
bool accepted = false;
for (label attempt = 0; attempt < 6; attempt++)
{
    bool finiteAndBounded = true;
    for(i=0;i<n;i++)
    {
        const scalar xi = xOld[i] + lambda*(xmma[i] - xOld[i]);
        if (!std::isfinite(xi) || xi < -1e-6 || xi > 1.0 + 1e-6)
        {
            finiteAndBounded = false;
            break;
        }
    }

    if (finiteAndBounded)
    {
        for(i=0;i<n;i++)
        {
            x[i] = Foam::min(Foam::max(xOld[i] + lambda*(xmma[i] - xOld[i]), scalar(0.0)), scalar(1.0));
        }
        accepted = true;
        break;
    }

    lambda *= 0.5;
}

if (!accepted)
{
    Info<< "Warning: MMA trial update rejected after backtracking; keeping previous design." << endl;
    for(i=0;i<n;i++)
    {
        x[i] = xOld[i];
    }
}
// Restore design variable in fixed zones (MMA updates all cells; filter needs correct x in solid/fluid/test)
if(solid_area)
{
    forAll(cells_solid, i) { x[cells_solid[i]] = 0; }
}
if(fluid_area)
{
    forAll(cells_fluid, i) { x[cells_fluid[i]] = 1.0; }
}
if(test_area)
{
    forAll(cells_test, i) { x[cells_test[i]] = 1.0; }
}

if(useCascadedFilter)
{
    #include "filter_x.H"
}
else
{
    #include "filter_x_legacy.H"
}
