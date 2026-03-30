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

Class
    Foam::diameterModels::coalescenceMomModels::constantCoalescenceMom

Description
    Constant coalescence kernel. Used for verification and validation of the
    coalescence formulation implemented in the populationBalanceMomModel class.

SourceFiles
    constantCoalescenceMom.C

\*---------------------------------------------------------------------------*/

#ifndef CONSTANTCOALESCENCEMOM_H
#define CONSTANTCOALESCENCEMOM_H

#include "coalescenceMomModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace coalescenceMomModels
{

/*---------------------------------------------------------------------------*\
                     Class constantCoalescenceMom Declaration
\*---------------------------------------------------------------------------*/

class constantCoalescenceMom
:
    public coalescenceMomModel
{
    // Private data

        //- constant coalescenceRate
        const dimensionedScalar rate_;


public:

    //- Runtime type information
    TypeName("constantMom");

    // Constructor

        constantCoalescenceMom
        (
            const populationBalanceMomModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~constantCoalescenceMom()
    {}


    // Member Functions

        //- Add to coalescenceRate
        virtual void addToCoalescenceRate
        (
            volScalarField& coalescenceRate,
            const volScalarField& di,
            const volScalarField& dj,
            const label i,
            const label j
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coalescenceMomModels
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
