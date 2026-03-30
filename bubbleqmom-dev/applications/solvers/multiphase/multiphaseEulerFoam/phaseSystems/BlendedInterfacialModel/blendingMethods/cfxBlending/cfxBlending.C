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

#include "cfxBlending.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{
    defineTypeNameAndDebug(cfxBlending, 0);

    addToRunTimeSelectionTable
    (
        blendingMethod,
        cfxBlending,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendingMethods::cfxBlending::cfxBlending
(
    const dictionary& dict,
    const wordList& phaseNames
)
:
    blendingMethod(dict),
    alphaMin_
        (
            "alphaMin_",
            dimless,
            dict.lookupOrDefault("alphaMin",1e-3)
        ),
    alphaMax_
        (
            "alphaMax_",
            dimless,
            dict.lookupOrDefault("alphaMax",0.8)
        ),
	continuousPhase_(dict.lookup("continuousPhase"))
{
      
    Info << "alphaMin: "
	 << alphaMin_.value() << endl
	 << "alphaMax: "
	 << alphaMax_.value() << endl;

    if
    (
       	alphaMin_
        > alphaMax_
    )
    {
        FatalErrorInFunction
            << "The supplied maximum volume fraction"
            << " is less than the minimum value."
            << endl << exit(FatalError);
    }
    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blendingMethods::cfxBlending::~cfxBlending()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::blendingMethods::cfxBlending::f1
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    if (continuousPhase_ == "")
    {
        FatalErrorInFunction
            << "The CFX phase  blending requires at least one continuous "
            << "phase, specified by keyword 'continuous' in phaseProperties."
            << exit(FatalError);
    }

    const fvMesh& mesh(phase1.mesh());
    
    volScalarField alpha_ = max (phase1,alphaMin_);
//    Info << "alphaMin: "
//         << alphaMin_.value() << endl
//         << "alphaMax: "
//         << alphaMax_.value() << endl;

    if (phase2.name() == continuousPhase_)
    {
        return min
            (
                max
                (
                    ((1 - alpha_)*alphaMax_)/((1 - alphaMax_)*alpha_),
                    alphaMin_/alpha_ 
                ),
                scalar (1)
            );
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "f1",
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar
                (
                    "f1",
                    dimless,
                    phase2.name() == continuousPhase_
                )
            )
        );
    }
}


Foam::tmp<Foam::volScalarField> Foam::blendingMethods::cfxBlending::f2
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    if (continuousPhase_ == "")
    {
        FatalErrorInFunction
            << "The CFX phase  blending requires at least one continuous "
            << "phase, specified by keyword 'continuous' in phaseProperties."
            << exit(FatalError);
    }
    const fvMesh& mesh(phase1.mesh());

    volScalarField alpha_ = max (phase2,alphaMin_);
	
    if (phase1.name() == continuousPhase_)
    {
        return min
            (
                max
                (
                    ((1 - alpha_)*alphaMax_)/((1 - alphaMax_)*alpha_ ),
                    alphaMin_/alpha_
                ), 
                scalar (1)
            );
    }
    else
    {
        return tmp<volScalarField> 
	    (
            new volScalarField
            (
                IOobject
                (
                    "f2",
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar
                (
                    "f2",
                    dimless,
                    phase1.name() == continuousPhase_
                )
            )
        );

    }

}

// ************************************************************************* //
