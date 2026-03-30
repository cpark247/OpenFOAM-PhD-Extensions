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
    Foam::diameterModels::daughterSizeDistributionMomModels::uniformBinaryMom

Description
    Uniform binary daughter size distribution. Used for verification and
    validation of the breakup formulation implemented in the
    populationBalanceMomModel class.

SourceFiles
    uniformBinaryMom.C

\*---------------------------------------------------------------------------*/

#ifndef uniformBinaryMom_H
#define uniformBinaryMom_H

#include "daughterSizeDistributionMomModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace daughterSizeDistributionMomModels
{

/*---------------------------------------------------------------------------*\
                      Class uniformBinaryMom Declaration
\*---------------------------------------------------------------------------*/

class uniformBinaryMom
:
    public daughterSizeDistributionMomModel
{
public:

    //- Runtime type information
    TypeName("uniformBinaryMom");


    // Constructor

        uniformBinaryMom
        (
            const breakupMomModel& breakup,
            const dictionary& dict
        );


    //- Destructor
    virtual ~uniformBinaryMom();


    // Member Functions

        //- Return total number of particles assigned to class i when a particle
        //  of class k breaks
        virtual void calcIntegrald
        (
            volScalarField& daughterIntegral,
            const volScalarField& di,
            const label i,
            const label k
        ) const;

        //- Return total number of particles assigned to class i when a particle
        //  of class k breaks
        virtual void calcIntegralm
        (
            volScalarField& daughterIntegral,
            const volScalarField& mi,
            const label i,
            const label k
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace daughterSizeDistributionMomModels
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
