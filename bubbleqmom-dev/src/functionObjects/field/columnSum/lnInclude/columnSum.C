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

#include "columnSum.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "meshStructureGI.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(columnSum, 0);
    addToRunTimeSelectionTable(functionObject, columnSum, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::meshStructureGI&
Foam::functionObjects::columnSum::meshAddressing(const polyMesh& mesh) const
{
    if (meshStructureGIPtr_.empty())
    {
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();

        // Count
        label sz = 0;
        forAll (patchIDs_, patchi)
        {
            sz += pbm[patchIDs_[patchi]].size();
        }

        // Fill
        labelList meshFaces(sz);
        sz = 0;
        forAll (patchIDs_, patchi)
        {
            label start = pbm[patchIDs_[patchi]].start();
            label size = pbm[patchIDs_[patchi]].size();
            for (label i = 0; i < size; ++i)
            {
                meshFaces[sz++] = start+i;
            }
        }

        if (sz == 0)
        {
            // TODO: If source patch is a cyclic it may have have been
            // converted to a processorCyclic for parallel runs

            WarningInFunction
                << "Requested patches have zero faces"
                << endl;
        }

        uindirectPrimitivePatch uip
        (
            UIndirectList<face>(mesh.faces(), meshFaces),
            mesh.points()
        );

        globalFaces_.reset(new globalIndex(uip.size()));
        globalEdges_.reset(new globalIndex(uip.nEdges()));
        globalPoints_.reset(new globalIndex(uip.nPoints()));
        meshStructureGIPtr_.reset
        (
            new meshStructureGI
            (
                mesh,
                uip,
                globalFaces_(),
                globalEdges_(),
                globalPoints_()
            )
        );
    }

    return *meshStructureGIPtr_;
}


const Foam::word Foam::functionObjects::columnSum::averageName
(
    const word& fieldName
) const
{
    return name() + ":columnSum(" + fieldName + ")";
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::columnSum::columnSum
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    patchIDs_(),
    fieldSet_(mesh_)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::columnSum::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    patchIDs_ =
        mesh_.boundaryMesh().patchSet
        (
            wordReList(dict.lookup("patches"))
        ).sortedToc();

    fieldSet_.read(dict);

    return true;
}


bool Foam::functionObjects::columnSum::execute()
{
    // Make fields up to date with current selection
    fieldSet_.updateSelection();
    const wordHashSet& fieldNameSet = fieldSet_.selectionNames();

    forAllConstIter (wordHashSet, fieldNameSet , nameIter)
    {
        columnSumField<scalar>(nameIter.key());
        //columnSumField<vector>(nameIter.key());
        //columnSumField<sphericalTensor>(nameIter.key());
        //columnSumField<symmTensor>(nameIter.key());
        //columnSumField<tensor>(nameIter.key());
    }

    return true;
}


bool Foam::functionObjects::columnSum::write()
{
    const wordHashSet& fieldNameSet = fieldSet_.selectionNames();
    forAllConstIter (wordHashSet, fieldNameSet, nameIter)
    {
        const word resultName("columnSum(" + nameIter.key() + ")");
        const regIOobject& obj =
            obr_.lookupObject<regIOobject>(averageName(nameIter.key()));

        if (obr_.foundObject<regIOobject>(averageName(nameIter.key())))
        {
            obj.write();
        }
    }

    return true;
}


// ************************************************************************* //
