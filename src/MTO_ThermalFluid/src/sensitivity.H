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

#include "filter_chainrule.H"

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
for(i=0;i<n;i++)
{
   x[i]=xmma[i];
}

#include "filter_x.H"
