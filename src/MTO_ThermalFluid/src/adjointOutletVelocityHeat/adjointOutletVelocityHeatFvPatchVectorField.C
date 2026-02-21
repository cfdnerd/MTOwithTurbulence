/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "adjointOutletVelocityHeatFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvPatchFieldMapper.H"
#include "RASModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointOutletVelocityHeatFvPatchVectorField::
adjointOutletVelocityHeatFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


Foam::adjointOutletVelocityHeatFvPatchVectorField::
adjointOutletVelocityHeatFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::adjointOutletVelocityHeatFvPatchVectorField::
adjointOutletVelocityHeatFvPatchVectorField
(
    const adjointOutletVelocityHeatFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::adjointOutletVelocityHeatFvPatchVectorField::
adjointOutletVelocityHeatFvPatchVectorField
(
    const adjointOutletVelocityHeatFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void Foam::adjointOutletVelocityHeatFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& phiap =
        patch().lookupPatchField<surfaceScalarField, scalar>("phib");

    const fvPatchField<vector>& Up =
        patch().lookupPatchField<volVectorField, vector>("U");

    const fvPatchField<vector>& Uap =
    	patch().lookupPatchField<volVectorField, vector>("Ub");

    const fvsPatchField<scalar>& phip =
    	patch().lookupPatchField<surfaceScalarField, scalar>("phi");
    	
    const volScalarField& nu = patch().boundaryMesh().mesh().lookupObject<volScalarField>("nu");
    scalarField nueff = nu.boundaryField()[patch().index()];
    if (patch().boundaryMesh().mesh().foundObject<volScalarField>("nut"))
    {
        const volScalarField& nut = patch().boundaryMesh().mesh().lookupObject<volScalarField>("nut");
        nueff += nut.boundaryField()[patch().index()];
    }

    const scalarField& deltainv =
        patch().deltaCoeffs(); // dist^(-1)

    const scalarField magSfSafe(max(patch().magSf(), Foam::VSMALL));
//Primal velocity, mag of normal component and tangential component
    scalarField Up_ns = phip/magSfSafe;

    vectorField Up_t = Up - (phip * patch().Sf())/(magSfSafe*magSfSafe);

//Tangential component of adjoint velocity in neighbouring node
    vectorField Uaneigh = Uap.patchInternalField();
    vectorField Uaneigh_n = (Uaneigh & patch().nf())*patch().nf();
    vectorField Uaneigh_t = Uaneigh - Uaneigh_n;

    vectorField Uap_t = (nueff*deltainv*Uaneigh_t) / (Up_ns+nueff*deltainv+Foam::VSMALL);

    vectorField Uap_n = (phiap * patch().Sf())/(magSfSafe*magSfSafe);

    operator==(Uap_t+Uap_n);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::adjointOutletVelocityHeatFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        adjointOutletVelocityHeatFvPatchVectorField
    );
}


// ************************************************************************* //
