/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "filterField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(filterField, 0);
    addToRunTimeSelectionTable(functionObject, filterField, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::filterField::calc()
{
    bool processed = false;

    #define processType(fieldType, none)                                       \
        processed = processed || calcFilter<fieldType>();
    //FOR_ALL_FIELD_TYPES(processType)	tensors cannot be filtered
    processType(scalar, none) \
    processType(vector, none)

    if (!processed)
    {
        cannotFindObject(fieldName_);
    }

    return processed;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::filterField::filterField
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, typeName),
    filterPtr_(LESfilterNSteps::New(mesh_,dict)),
    filter_(filterPtr_()),
    steps_(dict.lookupOrDefault<label>("steps",1)),
    filterWidth_(dict.lookup<scalar>("width")),
    deltaFilter_
    (
        IOobject
        (
            "deltaFilter",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimLength, 0.0),
        "zeroGradient"
    )
{
    volScalarField gridSize
    (
            IOobject
            (
                "gridSize",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimLength, 0.0),
            "zeroGradient"
    );
    gridSize.ref() = cbrt(mesh_.V());
    gridSize.correctBoundaryConditions();

    dimensionedScalar finalWidth("finalWidth", dimLength, filterWidth_);
    volScalarField localFilterFactor = finalWidth/(sqrt((scalar)steps_) * gridSize);

    Info<< "sub-filter factor min/max: "
        << min(localFilterFactor.internalField())
        << " / "
        << max(localFilterFactor.internalField())
        << endl;

    localFilterFactor.correctBoundaryConditions();
    deltaFilter_ = localFilterFactor * gridSize;

    Info<< "final filter width min/max: "
	<< min(deltaFilter_.internalField())*sqrt((scalar)steps_)
	<< " / "
	<< max(deltaFilter_.internalField())*sqrt((scalar)steps_)
	<< endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::filterField::~filterField()
{}


// ************************************************************************* //
