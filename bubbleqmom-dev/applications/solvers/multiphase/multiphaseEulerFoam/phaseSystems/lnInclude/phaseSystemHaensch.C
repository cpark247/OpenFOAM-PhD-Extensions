/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2021 OpenFOAM Foundation
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

#include "phaseSystemHaensch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseSystemHaensch, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseSystemHaensch::phaseSystemHaensch
(
    const fvMesh& mesh
)
:
    phaseSystem(mesh),
    haensch_ab_(lookupOrDefault<scalar>("haenschAb", 100))
{

    forAll(phases(), phasei)
    {
        const volScalarField& alphai = phases()[phasei];
        mesh_.setFluxRequired(alphai.name());

	phaseModel& phase = phases()[phasei];

	haensch_alphaMin_.insert
	(
	    phase.name(),
	    lookupOrDefault<scalar>
	    (
		IOobject::groupName("haenschAlphaClustMin", phase.name()),
		0.1
	    )
	);

	haensch_alphaMax_.insert
	(
	    phase.name(),
	    lookupOrDefault<scalar>
	    (
		IOobject::groupName("haenschAlphaClustMax", phase.name()),
		0.9
	    )
	);

	gamma_.insert
	(
	    phase.name(),
	    new volScalarField
	    (
		IOobject
		(
		    IOobject::groupName("gammaHaensch", phase.name()),
		    this->mesh().time().timeName(),
		    this->mesh(),
		    IOobject::NO_READ,
		    IOobject::AUTO_WRITE
		),
		this->mesh(),
		dimensionedScalar
		(
		    IOobject::groupName("gammaHaensch", phase.name()),
		    dimless,
		    0.0
		)
	    )
	);
    }

    Info << "Selected phaseSystemHaensch with a factor ab " << haensch_ab_ << " and critical cluster boundaries from " << haensch_alphaMin_ << " to " << haensch_alphaMax_ << ".\n";
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseSystemHaensch::~phaseSystemHaensch()
{}

// ************************************************************************* //
