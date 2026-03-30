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

#include "hydrodynamicMom.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace coalescenceMomModels
{
    defineTypeNameAndDebug(hydrodynamicMom, 0);
    addToRunTimeSelectionTable(coalescenceMomModel, hydrodynamicMom, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::coalescenceMomModels::hydrodynamicMom::hydrodynamicMom
(
    const populationBalanceMomModel& popBal,
    const dictionary& dict
)
:
    coalescenceMomModel(popBal, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::coalescenceMomModels::hydrodynamicMom::addToCoalescenceRate
(
    volScalarField& coalescenceRate,
    const volScalarField& di,
    const volScalarField& dj,
    const label i,
    const label j
)
{

	const dimensionedScalar oneSecond("oneSecond", dimTime, 1.);

    coalescenceRate +=
        pow3((di + dj))/oneSecond;
}


// ************************************************************************* //
