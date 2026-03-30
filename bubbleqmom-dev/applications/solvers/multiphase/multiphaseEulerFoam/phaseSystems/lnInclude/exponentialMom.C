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

#include "exponentialMom.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace breakupMomModels
{
    defineTypeNameAndDebug(exponentialMom, 0);
    addToRunTimeSelectionTable(breakupMomModel, exponentialMom, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::breakupMomModels::exponentialMom::exponentialMom
(
    const populationBalanceMomModel& popBal,
    const dictionary& dict
)
:
    breakupMomModel(popBal, dict),
    exponent_(readScalar(dict.lookup("exponent"))),
    C_(readScalar(dict.lookup("C")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::breakupMomModels::exponentialMom::setBreakupRate
(
    volScalarField& breakupRate,
    const volScalarField& di,
    const label i
)
{
	const dimensionedScalar oneSecond("oneSecond", dimTime, 1.);
	const dimensionedScalar cubDim("cubDim", dimLength*dimLength*dimLength, 1.);
    breakupRate =
        C_*exp(exponent_*(
						  pow3(di)*constant::mathematical::pi/6.0/cubDim
					     )
				)/oneSecond;
}


// ************************************************************************* //
