/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "alphaSum.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "cellSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(alphaSum, 0);
    addToRunTimeSelectionTable(functionObject, alphaSum, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::alphaSum::alphaSum
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    OS_("alphaSum"),
    cellSetName_("evalAlpha"),
    gasVol_(Zero)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::alphaSum::~alphaSum()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::alphaSum::read(const dictionary& dict)
{
    alphaName_ = IOobject::groupName("alpha",dict.lookup("phase"));

    Info<< indent
        << "- selecting cells using cellSet " << cellSetName_ << endl;
 
    cellSet selectedCells(mesh_, cellSetName_);
    cells_ = selectedCells.toc();
    return true;
}


bool Foam::functionObjects::alphaSum::execute()
{
    const volScalarField& alpha = mesh_.lookupObject<volScalarField>(alphaName_);
    const volScalarField::Internal alphaV = alpha.internalField()*mesh_.V();

    gasVol_ = 0.0;
    forAll(cells_,i)
    {
        const label celli = cells_[i];
        gasVol_ += alphaV[celli];
    }
    reduce(gasVol_, sumOp<scalar>());

    return true;
}


bool Foam::functionObjects::alphaSum::end()
{
    return true;
}


bool Foam::functionObjects::alphaSum::write()
{
    if (Pstream::master())
    {
	if(OS_.opened())
	{
	    OS_ << mesh_.time().value() << tab << gasVol_ << endl;
	}
    }

    return true;
}


// ************************************************************************* //
