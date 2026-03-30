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
    Foam::diameterModels::breakupMomModel

Description
    Base class for breakup models which give a total breakup rate and a separate
    daughter size distribution function.

SourceFiles
    breakupMomModel.C

\*---------------------------------------------------------------------------*/

#ifndef breakupMomModel_H
#define breakupMomModel_H

#include "populationBalanceMomModel.H"
#include "daughterSizeDistributionMomModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{

/*---------------------------------------------------------------------------*\
                        Class breakupMomModel Declaration
\*---------------------------------------------------------------------------*/

class breakupMomModel
{
protected:

    // Protected data

        //- Reference to the populationBalanceMomModel
        const populationBalanceMomModel& popBal_;

        //- Dictionary
        dictionary dict_;

        //- Daughter size distribution model
        autoPtr<daughterSizeDistributionMomModel> dsd_;


public:

    //- Runtime type information
    TypeName("breakupMomModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            breakupMomModel,
            dictionary,
            (
                const populationBalanceMomModel& popBal,
                const dictionary& dict
            ),
            (popBal, dict)
        );


    //- Class used for the read-construction of
    //  PtrLists of breakup models
    class iNew
    {
        const populationBalanceMomModel& popBal_;

    public:

        iNew(const populationBalanceMomModel& popBal)
        :
            popBal_(popBal)
        {}

        autoPtr<breakupMomModel> operator()(Istream& is) const
        {
            word type(is);
            dictionary dict(is);
            return breakupMomModel::New(type, popBal_, dict);
        }
    };


    // Constructor

        breakupMomModel
        (
            const populationBalanceMomModel& popBal,
            const dictionary& dict
        );

        autoPtr<breakupMomModel> clone() const
        {
            NotImplemented;
            return autoPtr<breakupMomModel>(nullptr);
        }


    // Selector

        static autoPtr<breakupMomModel> New
        (
            const word& type,
            const populationBalanceMomModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~breakupMomModel()
    {}


    // Member Functions

        //- Return reference to the populationBalanceMomModel
        const populationBalanceMomModel& popBal() const
        {
            return popBal_;
        }

        //- Return const-reference to daughter size distribution pointer
        const autoPtr<daughterSizeDistributionMomModel>& dsdPtr() const
        {
            return dsd_;
        }

        //- Return reference to daughter size distribution pointer
        autoPtr<daughterSizeDistributionMomModel>& dsdPtr()
        {
            return dsd_;
        }

        //- Correct diameter independent expressions
        virtual void precompute();

        //- Set total breakupRate
        virtual void setBreakupRate
        (
            volScalarField& breakupRate,
            const volScalarField& di,
            const label i
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
