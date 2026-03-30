/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
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

#include "haenschSurfaceTension.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceTensionModels
{
    defineTypeNameAndDebug(haenschSurfaceTension, 0);
    addToRunTimeSelectionTable
    (
        surfaceTensionModel,
        haenschSurfaceTension,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceTensionModels::haenschSurfaceTension::
haenschSurfaceTension
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    surfaceTensionModel(dict, pair, registerObject),
    sigma_("sigma", dimSigma, dict),
    haensch_ab_(dict.lookupOrDefault<scalar>("haenschAb", 100)),
    haensch_alphaMin_
    (
	dict.lookupOrDefault<scalar>
	(
	    IOobject::groupName("haenschAlphaClustMin", this->pair_.phase1().name()),
	    0.1
	)
    ),
    haensch_alphaMax_
    (
	dict.lookupOrDefault<scalar>
	(
	    IOobject::groupName("haenschAlphaClustMax", this->pair_.phase1().name()),
	    0.9
	)
    )
{
    Info << "Selected Haensch surface tension blending with a factor ab " << haensch_ab_ << " and critical cluster boundaries from " << haensch_alphaMin_ << " to " << haensch_alphaMax_ << ".\n";
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceTensionModels::haenschSurfaceTension::
~haenschSurfaceTension()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::surfaceTensionModels::haenschSurfaceTension::sigma() const
{
    tmp<volScalarField> gamma =
	(0.5*tanh(haensch_ab_*(this->pair_.phase1()-haensch_alphaMin_))+0.5) *
	(0.5*tanh(haensch_ab_*(haensch_alphaMax_-this->pair_.phase1()))+0.5)*sigma_;

    return gamma;
}

Foam::tmp<Foam::scalarField>
Foam::surfaceTensionModels::haenschSurfaceTension::sigma
(
    label patchi
) const
{
    volScalarField gamma =
	(0.5*tanh(haensch_ab_*(this->pair_.phase1()-haensch_alphaMin_))+0.5) *
	(0.5*tanh(haensch_ab_*(haensch_alphaMax_-this->pair_.phase1()))+0.5)*sigma_;

    return tmp<scalarField>
    (
	gamma.boundaryField()[patchi]
    );
}


// ************************************************************************* //
