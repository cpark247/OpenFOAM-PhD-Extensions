/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "uniformBinaryMom.H"
#include "addToRunTimeSelectionTable.H"
#include "breakupMomModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace daughterSizeDistributionMomModels
{
    defineTypeNameAndDebug(uniformBinaryMom, 0);
    addToRunTimeSelectionTable
    (
        daughterSizeDistributionMomModel,
        uniformBinaryMom,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::daughterSizeDistributionMomModels::uniformBinaryMom::
uniformBinaryMom
(
    const breakupMomModel& breakup,
    const dictionary& dict
)
:
    daughterSizeDistributionMomModel(breakup, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::daughterSizeDistributionMomModels::uniformBinaryMom::
~uniformBinaryMom()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void 
Foam::diameterModels::daughterSizeDistributionMomModels::
uniformBinaryMom::calcIntegrald
(
    volScalarField& daughterIntegral,
    const volScalarField& di,
    const label i,
    const label k
) const
{
	daughterIntegral =
	                   (pow(cbrt(0.5)*di, k+1))/
	                   (k+1.);
}


void 
Foam::diameterModels::daughterSizeDistributionMomModels::
uniformBinaryMom::calcIntegralm
(
    volScalarField& daughterIntegral,
    const volScalarField& mi,
    const label i,
    const label k
) const
{
	daughterIntegral = 6.0*pow(mi, k)/(k + 3.0);
}


// ************************************************************************* //
