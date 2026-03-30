/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
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
    Foam::functionObjects::filterField

Description
    Gives back a filtered version of a field.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    filterField.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_filterField_H
#define functionObjects_filterField_H

#include "fieldExpression.H"
#include "LESfilterNSteps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                             Class filterField Declaration
\*---------------------------------------------------------------------------*/

class filterField
:
    public fieldExpression
{
    // Private data

	//- using LESfilter with multiple substeps
	autoPtr<LESfilterNSteps> filterPtr_;
	LESfilterNSteps& filter_;

	//- number of substeps
	label steps_;

	//- filter width
	scalar filterWidth_;
	volScalarField deltaFilter_;


    // Private Member Functions

        //- Calculate the filtered field and register the result
        template<class Type>
        bool calcFilter();

        //- Calculate the filtered field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("filterField");


    // Constructors

        //- Construct from Time and dictionary
        filterField
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~filterField();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "filterFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
