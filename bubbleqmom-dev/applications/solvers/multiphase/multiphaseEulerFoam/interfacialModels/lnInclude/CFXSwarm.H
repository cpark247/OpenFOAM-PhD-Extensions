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
    Foam::swarmCorrections::CFXSwarm

Description
    scales the drag force in the same way as in Ansys CFX

SourceFiles
    CFXSwarm.C

\*---------------------------------------------------------------------------*/

#ifndef CFXSwarm_H
#define CFXSwarm_H

#include "swarmCorrection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace swarmCorrections
{

/*---------------------------------------------------------------------------*\
                         Class CFXSwarm Declaration
\*---------------------------------------------------------------------------*/

class CFXSwarm
:
    public swarmCorrection
{
    // Private data

	//- Residual correction factor
	const dimensionedScalar residualSwarmFactor_;

        //- Upper volume fraction limit for drag force
        const dimensionedScalar limit_;
	
public:

    //- Runtime type information
    TypeName("CFXDrag");


    // Constructors

        //- Construct from a dictionary and a phase pair
        CFXSwarm
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~CFXSwarm();


    // Member Functions

        //- Swarm correction coefficient
        virtual tmp<volScalarField> Cs() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace swarmCorrections
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
