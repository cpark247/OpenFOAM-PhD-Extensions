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
    Foam::diameterModels::breakupMomModels::exponentialMom

Description
    Exponential kernel. Used for verification and validation of the breakup
    formulation implemented in the populationBalanceMomModel class.

SourceFiles
    exponentialMom.C

\*---------------------------------------------------------------------------*/

#ifndef exponentialMom_H
#define exponentialMom_H

#include "breakupMomModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace breakupMomModels
{

/*---------------------------------------------------------------------------*\
                         Class exponentialMom Declaration
\*---------------------------------------------------------------------------*/

class exponentialMom
:
    public breakupMomModel
{
    // Private data

        //- Exponent
        scalar exponent_;

        //- Coefficient
        scalar C_;


public:

    //- Runtime type information
    TypeName("exponentialMom");

    // Constructor

        exponentialMom
        (
            const populationBalanceMomModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~exponentialMom()
    {}


    // Member Functions

        //- Set total breakupRate
        virtual void setBreakupRate
        (
            volScalarField& breakupRate,
            const volScalarField& di,
            const label i
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace breakupMomModels
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
