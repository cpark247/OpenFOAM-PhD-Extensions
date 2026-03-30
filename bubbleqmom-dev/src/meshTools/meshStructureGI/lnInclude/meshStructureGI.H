/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    Foam::meshStructureGI

Description
    Detect extruded mesh structure given a set of patch faces.
    Can use global indexing like OF-2006.

SourceFiles
    meshStructureGI.C

\*---------------------------------------------------------------------------*/

#ifndef meshStructureGI_H
#define meshStructureGI_H

#include "labelList.H"
#include "uindirectPrimitivePatch.H"
#include "className.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;


/*---------------------------------------------------------------------------*\
                         Class meshStructureGI Declaration
\*---------------------------------------------------------------------------*/

class meshStructureGI
{
    // Private data

        //- Cell to patch face
        labelList cellToPatchFaceAddressing_;

        //- Cell to layer
        labelList cellLayer_;

        //- Face to patch face
        labelList faceToPatchFaceAddressing_;

        //- Face to patch edge
        labelList faceToPatchEdgeAddressing_;

        //- Face to layer
        labelList faceLayer_;

        //- Point to patch point
        labelList pointToPatchPointAddressing_;

        //- Point to layer
        labelList pointLayer_;

        //- Is mesh structured?
        bool structured_;


   // Private Member Functions

        //- Is cell structured
        bool isStructuredCell
        (
            const polyMesh& mesh,
            const label layerI,
            const label celli
        ) const;

        //- Calculate all maps.
        void correct
        (
            const polyMesh& mesh,
            const uindirectPrimitivePatch& pp,
            const globalIndex& globalFaces,
            const globalIndex& globalEdges,
            const globalIndex& globalPoints
        );


public:

    // Declare name of the class and its debug switch
    ClassName("meshStructureGI");

    // Constructors

	//- Construct from mesh and faces in mesh and global indexing
        meshStructureGI
        (
            const polyMesh& mesh,
            const uindirectPrimitivePatch&,
            const globalIndex& globalFaces,
            const globalIndex& globalEdges,
            const globalIndex& globalPoints
        );


    // Member Functions

        //- Is mesh structured?
        inline bool structured() const;

        //- Cell to patch face
        inline const labelList& cellToPatchFaceAddressing() const;

        //- Cell to patch face
        inline labelList& cellToPatchFaceAddressing();

        //- Cell to layer
        inline const labelList& cellLayer() const;

        //- Cell to layer
        inline labelList& cellLayer();

        //- Face to patch face
        inline const labelList& faceToPatchFaceAddressing() const;

        //- Face to patch face
        inline labelList& faceToPatchFaceAddressing();

        //- Face to patch edge
        inline const labelList& faceToPatchEdgeAddressing() const;

        //- Face to patch edge
        inline labelList& faceToPatchEdgeAddressing();

        //- Face to layer
        inline const labelList& faceLayer() const;

        //- Face to layer
        inline labelList& faceLayer();

        //- Point to patch point
        inline const labelList& pointToPatchPointAddressing() const;

        //- Point to patch point
        inline labelList& pointToPatchPointAddressing();

        //- Point to layer
        inline const labelList& pointLayer() const;

        //- Point to layer
        inline labelList& pointLayer();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshStructureGII.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
