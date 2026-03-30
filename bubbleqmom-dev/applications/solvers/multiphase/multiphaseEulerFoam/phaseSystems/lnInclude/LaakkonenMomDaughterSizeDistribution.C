/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

#include "LaakkonenMomDaughterSizeDistribution.H"
#include "addToRunTimeSelectionTable.H"
#include "breakupMomModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace daughterSizeDistributionMomModels
{
    defineTypeNameAndDebug(LaakkonenMomDaughterSizeDistribution, 0);
    addToRunTimeSelectionTable
    (
        daughterSizeDistributionMomModel,
        LaakkonenMomDaughterSizeDistribution,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::daughterSizeDistributionMomModels::
LaakkonenMomDaughterSizeDistribution::LaakkonenMomDaughterSizeDistribution
(
    const breakupMomModel& breakup,
    const dictionary& dict
)
:
    daughterSizeDistributionMomModel(breakup, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::daughterSizeDistributionMomModels::
LaakkonenMomDaughterSizeDistribution::~LaakkonenMomDaughterSizeDistribution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::diameterModels::daughterSizeDistributionMomModels::
LaakkonenMomDaughterSizeDistribution::calcIntegrald
(
    volScalarField& daughterIntegral,
    const volScalarField& di,
    const label i,
    const label k
) const
{
	daughterIntegral =
	        (3240. * pow(di,k))/
	        ((k+9.)*(k+12.)*(k+15.));
}

void
Foam::diameterModels::daughterSizeDistributionMomModels::
LaakkonenMomDaughterSizeDistribution::calcIntegralm
(
    volScalarField& daughterIntegral,
    const volScalarField& mi,
    const label i,
    const label k
) const
{
	daughterIntegral =
	        (120. * pow(mi,k))/
	        ((k+3.)*(k+4.)*(k+5.));
}


// ************************************************************************* //
