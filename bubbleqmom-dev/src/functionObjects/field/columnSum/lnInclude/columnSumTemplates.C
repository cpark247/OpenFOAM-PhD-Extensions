/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "volFields.H"
#include "meshStructureGI.H"
#include "globalIndex.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::columnSum::columnSumField
(
    const word& fieldName
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (foundObject<fieldType>(fieldName))
    {
	const fieldType& fld = lookupObject<fieldType>(fieldName);
	const typename fieldType::Internal gradFld = mag(fvc::grad(fld));

        const word resultName(averageName(fieldName));

        if (!obr_.foundObject<fieldType>(resultName))
        {
	    IOobject resHeader
	    (
	        resultName,
	        fld.mesh().time().timeName(),
	        fld.mesh(),
	        IOobject::NO_READ,
	        IOobject::NO_WRITE
	    );
	    fieldType* resPtr(new fieldType(resHeader, fld));
	    obr_.objectRegistry::store(resPtr);
        }
        fieldType& res = obr_.lookupObjectRef<fieldType>(resultName);

        const meshStructureGI& ms = meshAddressing(fld.mesh());
        if (globalFaces_.empty())
        {
            return false;
        }

        const labelList& cellToPatchFace = ms.cellToPatchFaceAddressing();

        // Brute force: collect per-global-patchface on all processors
        Field<Type> regionField(globalFaces_().size(), Zero);
        scalarList regionCount(globalFaces_().size(), Zero);

	// Use cell volume as weighting factor so this can be used on
	// structured but non-equidistant mesh 
        forAll(cellToPatchFace, celli)
        {
            const label regioni = cellToPatchFace[celli];
            regionField[regioni] += gradFld[celli]*fld.mesh().V()[celli];
            regionCount[regioni] += fld.mesh().V()[celli];
        }

        // Global sum
        Pstream::listCombineGather(regionField, plusEqOp<Type>());
        Pstream::listCombineScatter(regionField);
        Pstream::listCombineGather(regionCount, plusEqOp<scalar>());
        Pstream::listCombineScatter(regionCount);

        forAll(regionField, regioni)
        {
            regionField[regioni] /= regionCount[regioni];
        }

        // And send result back
        forAll(cellToPatchFace, celli)
        {
            const label regioni = cellToPatchFace[celli];
            res[celli] = regionField[regioni];
        }
        res.correctBoundaryConditions();

        return true;
    }

    return false;
}


// ************************************************************************* //
