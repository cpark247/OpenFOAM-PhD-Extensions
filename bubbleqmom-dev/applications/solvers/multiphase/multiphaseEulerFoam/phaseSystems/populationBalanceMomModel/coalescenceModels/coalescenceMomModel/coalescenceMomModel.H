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
    Foam::diameterModels::coalescenceMomModel

Description
    Base class for coalescence models.

SourceFiles
    coalescenceMomModel.C

\*---------------------------------------------------------------------------*/

#ifndef coalescenceMomModel_H
#define coalescenceMomModel_H

#include "populationBalanceMomModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{

/*---------------------------------------------------------------------------*\
                      Class coalescenceMomModel Declaration
\*---------------------------------------------------------------------------*/

class coalescenceMomModel
{
protected:

    // Protected data

        //- Reference to the populationBalanceMomModel
        const populationBalanceMomModel& popBal_;


public:

    //- Runtime type information
    TypeName("coalescenceMomModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            coalescenceMomModel,
            dictionary,
            (
                const populationBalanceMomModel& popBal,
                const dictionary& dict
            ),
            (popBal, dict)
        );


    //- Class used for the read-construction of
    //  PtrLists of coalescence models
    class iNew
    {
        const populationBalanceMomModel& popBal_;

    public:

        iNew(const populationBalanceMomModel& popBal)
        :
            popBal_(popBal)
        {}

        autoPtr<coalescenceMomModel> operator()(Istream& is) const
        {
            word type(is);
            dictionary dict(is);
            return coalescenceMomModel::New(type, popBal_, dict);
        }
    };


    // Constructor

        coalescenceMomModel
        (
            const populationBalanceMomModel& popBal,
            const dictionary& dict
        );

        autoPtr<coalescenceMomModel> clone() const
        {
            NotImplemented;
            return autoPtr<coalescenceMomModel>(nullptr);
        }


    // Selector

        static autoPtr<coalescenceMomModel> New
        (
            const word& type,
            const populationBalanceMomModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~coalescenceMomModel()
    {}


    // Member Functions

        //- Correct diameter independent expressions
        virtual void precompute();

        //- Add to coalescenceRate
        virtual void addToCoalescenceRate
        (
            volScalarField& coalescenceRate,
            const volScalarField& di,
            const volScalarField& dj,
            const label i,
            const label j
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
