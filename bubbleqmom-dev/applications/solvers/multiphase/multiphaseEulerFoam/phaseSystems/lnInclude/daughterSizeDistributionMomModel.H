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
    Foam::diameterModels::daughterSizeDistributionMomModel

Description
    Base class for daughter size distribution models.

SourceFiles
    daughterSizeDistributionMomModel.C

\*---------------------------------------------------------------------------*/

#ifndef daughterSizeDistributionMomModel_H
#define daughterSizeDistributionMomModel_H

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "fvMesh.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{

class breakupMomModel;

/*---------------------------------------------------------------------------*\
                Class daughterSizeDistributionMomModel Declaration
\*---------------------------------------------------------------------------*/

class daughterSizeDistributionMomModel
{
protected:

    // Protected data

        //- Reference to the breakupMomModel
        const breakupMomModel& breakup_;


public:

    //- Runtime type information
    TypeName("daughterSizeDistributionMomModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            daughterSizeDistributionMomModel,
            dictionary,
            (
                const breakupMomModel& breakup,
                const dictionary& dict
            ),
            (breakup, dict)
        );


    // Constructor

        daughterSizeDistributionMomModel
        (
            const breakupMomModel& breakup,
            const dictionary& dict
        );


    //- Destructor
    virtual ~daughterSizeDistributionMomModel();


    // Selectors

        static autoPtr<daughterSizeDistributionMomModel> New
        (
            const breakupMomModel& breakup,
            const dictionary& dict
        );


    // Member Functions

        //- Correct field independent expressions
        void precompute();

        //- Calculate the rate \f$ \int_0^{d_i} d^k \beta(d|d_i) dd \f$ where k is the moment order
        virtual void calcIntegrald
        (
            volScalarField& daughterIntegral,
            const volScalarField& di,
            const label i,
            const label k
        ) const = 0;

        //- Calculate the rate \f$ \int_0^{m_i} m^k \beta(m|m_i) dm \f$ where k is the moment order
        virtual void calcIntegralm
        (
            volScalarField& daughterIntegral,
            const volScalarField& mi,
            const label i,
            const label k
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
