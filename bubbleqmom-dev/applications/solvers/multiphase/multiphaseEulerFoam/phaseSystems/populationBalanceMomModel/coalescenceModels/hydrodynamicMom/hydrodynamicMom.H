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
    Foam::diameterModels::coalescenceMomModels::hydrodynamicMom

Description
    Hydrodynamic kernel. Used for verification and validation of the coalescence
    formulation implemented in the populationBalanceMomModel class.

SourceFiles
    hydrodynamicMom.C

\*---------------------------------------------------------------------------*/

#ifndef HYDRODYNAMICMOM_H
#define HYDRODYNAMICMOM_H

#include "coalescenceMomModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace coalescenceMomModels
{

/*---------------------------------------------------------------------------*\
                        Class hydrodynamicMom Declaration
\*---------------------------------------------------------------------------*/

class hydrodynamicMom
:
    public coalescenceMomModel
{
public:

    //- Runtime type information
    TypeName("hydrodynamicMom");

    // Constructor

        hydrodynamicMom
        (
            const populationBalanceMomModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~hydrodynamicMom()
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
