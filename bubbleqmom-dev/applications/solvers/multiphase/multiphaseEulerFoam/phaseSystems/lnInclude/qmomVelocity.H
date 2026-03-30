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
    Foam::diameterModels::qmomVelocity

Description
    This diameterModel is intended for use with a populationBalanceModel in
    order to simulate polydispersed bubbly or particulate flows. It can hold any
    number of sizeGroups from which the Sauter mean diameter is calculated. It
    can also be used as a diameterModel without a populationBalance and would
    then behave like a constantDiameter model. In this case, some arbitrary name
    must be entered for the populationBalance keyword.

Usage
    \table
        Property          | Description
        populationBalance | Name of the corresponding populationBalance
        sizeGroups        | List of sizeGroups
    \endtable

    Example
    \verbatim
    diameterModel   qmomVelocity;
    qmomVelocityCoeffs
    {
        populationBalance    bubbles;

        shapeModel	spherical;

    }
    \endverbatim

See also
    Foam::diameterModels::sizeGroup
    Foam::diameterModels::populationBalanceModel

SourceFiles
    qmomVelocity.C

\*---------------------------------------------------------------------------*/

#ifndef QMOMVELOCITY_H
#define QMOMVELOCITY_H

#include "diameterModel.H"
#include "multivariateScheme.H"
#include "convectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{

class sizeGroup;

/*---------------------------------------------------------------------------*\
                        Class qmomVelocity Declaration
\*---------------------------------------------------------------------------*/

class qmomVelocity
:
    public diameterModel
{
    // Private data

        //- Name of the populationBalance this qmomVelocity belongs to
        word popBalName_;

        //- Number-based Sauter-mean diameter of the phase
        volScalarField d_;

        //- Mass transfer rate
        volScalarField dmdt_;




public:

    //- Runtime type information
    TypeName("qmomVelocity");


    // Constructors

    //- Construct from components
    qmomVelocity
    (
        const dictionary& diameterProperties,
        const phaseModel& phase
    );


    //- Destructor
    virtual ~qmomVelocity();


    // Member Functions

        //- Return name of populationBalance this qmomVelocity belongs to
        inline const word& popBalName() const;

        //- Return const-reference to mass transfer rate
        inline const volScalarField& dmdt() const;

        //- Return reference to mass transfer rate
        inline volScalarField& dmdt();

        //- Read diameterProperties dictionary
        virtual bool read(const dictionary& diameterProperties);

        //- Return the Sauter-mean diameter
        virtual tmp<volScalarField> d() const;

        //- Return the surface area per unit volume
        virtual tmp<volScalarField> a() const;

        //- update the Sauter-mean diameter
        virtual void updateSauterMeanDiam(const volScalarField& inputMeanDiameter);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "qmomVelocityI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
