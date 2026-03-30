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

\*---------------------------------------------------------------------------*/

#include "momentsTruncFilter.H"
#include "addToRunTimeSelectionTable.H"
#include "calculatedFvPatchFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "bound.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(momentsTruncFilter, 0);
    addToRunTimeSelectionTable(LESfilterNSteps, momentsTruncFilter, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::momentsTruncFilter::momentsTruncFilter(const fvMesh& mesh, scalar secondMCoeff)
:
    LESfilterNSteps(mesh),
    secondMCoeff_(secondMCoeff)
{
    
}


Foam::momentsTruncFilter::momentsTruncFilter(const fvMesh& mesh, const dictionary& bd)
:
    LESfilterNSteps(mesh),
    secondMCoeff_(readScalar(bd.subDict(type() + "Coeffs").lookup("secondMCoeff")))//,
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::momentsTruncFilter::read(const dictionary& bd)
{
    bd.subDict(type() + "Coeffs").lookup("secondMCoeff") >> secondMCoeff_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::momentsTruncFilter::operator()
(
    const tmp<volScalarField>& unFilteredField,
    const tmp<volScalarField>& filterDelta
)
{
    tmp<volScalarField> filteredField =
            unFilteredField() + fvc::laplacian(secondMCoeff_ * pow(filterDelta, 2.0), unFilteredField()); 
    
    unFilteredField.clear();
    filterDelta.clear();

    return filteredField;
}

Foam::tmp<Foam::volVectorField> Foam::momentsTruncFilter::operator()
(
    const tmp<volVectorField>& unFilteredField,
    const tmp<volScalarField>& filterDelta
)
{
    tmp<volVectorField> filteredField =
            unFilteredField() + fvc::laplacian(secondMCoeff_ * pow(filterDelta, 2.0), unFilteredField()); 
    
    unFilteredField.clear();
    filterDelta.clear();

    return filteredField;
}  

Foam::tmp<Foam::volScalarField> Foam::momentsTruncFilter::operator()
(
    const tmp<volScalarField>& unFilteredField,
    const tmp<volScalarField>& filterDelta,
    const label n
)
{
    volScalarField tmpField = unFilteredField();
    volScalarField delta = filterDelta(); 

    volScalarField subFilteredField = tmpField; 
    word gradScheme = "grad(laplaceFilter)";
    word divScheme = "div(laplaceFilter)"; 
    Info << "Number of filter sub steps: " << n << endl;
    for(int i=0; i<n; i++)
    { 
        //filtered = tmpField + fvc::laplacian(secondMCoeff_ * pow(delta, 2.0), tmpField);
        // div(grad) more robust than laplacian
        subFilteredField = tmpField + fvc::div(secondMCoeff_ * pow(delta, 2.0) * fvc::grad(tmpField,gradScheme),divScheme);
        subFilteredField.correctBoundaryConditions();
        tmpField = subFilteredField;
    }
  	
    scalar maxFiltered = max(subFilteredField.primitiveField());
	scalar minFiltered = min(subFilteredField.primitiveField());

    reduce(maxFiltered, maxOp<scalar>());
    reduce(minFiltered, minOp<scalar>());
 
    Info << "filtered min/max: " << minFiltered << "/" << maxFiltered << endl;
    //tmp<volScalarField> filteredField = subFilteredField + fvc::laplacian(secondMCoeff_ * pow(delta, 2.0), subFilteredField);
    tmp<volScalarField> filteredField = subFilteredField + fvc::div(secondMCoeff_ * pow(delta, 2.0) * fvc::grad(subFilteredField,gradScheme),divScheme);
    
    unFilteredField.clear();
    filterDelta.clear();
    
    return filteredField;
}

Foam::tmp<Foam::volVectorField> Foam::momentsTruncFilter::operator()
(
    const tmp<volVectorField>& unFilteredField,
    const tmp<volScalarField>& filterDelta,
    const label n
)
{
    volVectorField tmpField = unFilteredField();
    volScalarField delta = filterDelta(); 

    volVectorField subFilteredField = tmpField; 
    const word gradScheme("grad(laplaceFilter)");
    const word divScheme("div(laplaceFilter)"); 
    Info << "Number of filter sub steps: " << n << endl;
    for(int i=0; i<n; i++)
    { 
        //subFilteredField = tmpField + fvc::laplacian(secondMCoeff_ * pow(delta, 2.0), tmpField);
        // div(grad) more robust than laplacian
        subFilteredField = tmpField + fvc::div(secondMCoeff_ * pow(delta, 2.0) * fvc::grad(tmpField,gradScheme),divScheme);
        subFilteredField.correctBoundaryConditions();
        tmpField = subFilteredField;
    }
  	
    scalar maxFiltered = max(mag(subFilteredField.primitiveField()));
	scalar minFiltered = min(mag(subFilteredField.primitiveField()));

    reduce(maxFiltered, maxOp<scalar>());
    reduce(minFiltered, minOp<scalar>());
 
    Info << "filtered min/max: " << minFiltered << "/" << maxFiltered << endl;
    //tmp<volVectorField> filteredField = subFilteredField + fvc::laplacian(secondMCoeff_ * pow(delta, 2.0), subFilteredField);
    tmp<volVectorField> filteredField = subFilteredField + fvc::div(secondMCoeff_ * pow(delta, 2.0) * fvc::grad(subFilteredField,gradScheme),divScheme);
    unFilteredField.clear();
    filterDelta.clear();
    
    return filteredField;
}
// ************************************************************************* //
