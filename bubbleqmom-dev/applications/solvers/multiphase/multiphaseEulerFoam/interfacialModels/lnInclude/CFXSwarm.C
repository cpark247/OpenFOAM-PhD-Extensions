/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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

#include "CFXSwarm.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace swarmCorrections
{
    defineTypeNameAndDebug(CFXSwarm, 0);
    addToRunTimeSelectionTable
    (
        swarmCorrection,
        CFXSwarm,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swarmCorrections::CFXSwarm::CFXSwarm
(
    const dictionary& dict,
    const phasePair& pair
)
:
    swarmCorrection(dict, pair),
    residualSwarmFactor_("residualSwarmFactor", dimless, dict),
    limit_("limit", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::swarmCorrections::CFXSwarm::~CFXSwarm()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::swarmCorrections::CFXSwarm::Cs() const
{
    Info<<"Computing drag correction factor in CFXSwarm"<<endl;
    return 
	neg(this->pair_.dispersed()-limit_)
	+pos0(this->pair_.dispersed()-limit_)
	*max
	(
		(1.0-this->pair_.dispersed())
		/(1.0-limit_)
		*limit_
		/max(this->pair_.dispersed(),this->pair_.dispersed().residualAlpha()),
		residualSwarmFactor_
	);
}


// ************************************************************************* //
