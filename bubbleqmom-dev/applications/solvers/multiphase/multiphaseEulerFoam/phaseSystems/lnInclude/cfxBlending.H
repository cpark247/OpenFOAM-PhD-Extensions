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

Class
    Foam::blendingMethods::cfxBlending

Description
    Blending method based on CFX blending for interfacial area density of 
    the particle  model. 

    r =
    max (r, r_{min}), for r < r_{max}
    max (\fraq{1-r}{1-r_{max}}*r_{max},r_{min}) for r >= r_{max}
    
    

SourceFiles
    cfxBlending.C

\*---------------------------------------------------------------------------*/

#ifndef cfxBlending_H
#define cfxBlending_H

#include "blendingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{

/*---------------------------------------------------------------------------*\
                         Class cfxBlending Declaration
\*---------------------------------------------------------------------------*/

class cfxBlending
:
    public blendingMethod
{
    // Private data

        //- Minimum volume fraction, CFX 1e-6
        const dimensionedScalar alphaMin_;

        //- Maximum packing, CFX 0.8
        const dimensionedScalar alphaMax_;

	//- Name of the continuous phase
        const word continuousPhase_;

public:

    //- Runtime type information
    TypeName("cfxBlending");


    // Constructors

        //- Construct from a dictionary and two phases
        cfxBlending
        (
            const dictionary& dict,
            const wordList& phaseNames
        );


    //- Destructor
    ~cfxBlending();


    // Member Functions

        //- Factor for primary phase
        virtual tmp<volScalarField> f1
        (
            const phaseModel& phase1,
            const phaseModel& phase2
        ) const;

        //- Factor for secondary phase
        virtual tmp<volScalarField> f2
        (
            const phaseModel& phase1,
            const phaseModel& phase2
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blendingMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
