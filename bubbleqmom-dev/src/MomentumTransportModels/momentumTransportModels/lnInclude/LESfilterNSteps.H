/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::LESfilterNSteps

Description
    Abstract class for LES filters which allows multiple filtering steps.

SourceFiles
    LESfilterNSteps.C

\*---------------------------------------------------------------------------*/

#ifndef LESfilterNSteps_H
#define LESfilterNSteps_H

#include "volFields.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class LESfilterNSteps Declaration
\*---------------------------------------------------------------------------*/

class LESfilterNSteps
{
    // Private data

        const fvMesh& mesh_;


    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        LESfilterNSteps(const LESfilterNSteps&);
        void operator=(const LESfilterNSteps&);


public:

    //- Runtime type information
    TypeName("LESfilterNSteps");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            LESfilterNSteps,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& LESfilterNStepsDict
            ),
            (mesh, LESfilterNStepsDict)
        );


    // Constructors

        //- Construct from components
        LESfilterNSteps(const fvMesh& mesh)
        :
            mesh_(mesh)
        {}


    // Selectors

        //- Return a reference to the selected LES filter
        static autoPtr<LESfilterNSteps> New
        (
            const fvMesh&,
            const dictionary&
        );


    //- Destructor
    virtual ~LESfilterNSteps()
    {}


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Read the LESfilterNSteps dictionary
        virtual void read(const dictionary&) = 0;


    // Member Operators

        virtual tmp<volScalarField> operator()
        (
            const tmp<volScalarField>&,
            const tmp<volScalarField>&
        ) = 0;

        virtual tmp<volVectorField> operator()
        (
            const tmp<volVectorField>&,
            const tmp<volScalarField>&
        ) = 0;

        virtual tmp<volScalarField> operator()
        (
            const tmp<volScalarField>&,
            const tmp<volScalarField>&,
            const label
        ) = 0;

        virtual tmp<volVectorField> operator()
        (
            const tmp<volVectorField>&,
            const tmp<volScalarField>&,
            const label
        ) = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
