/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::patchInterpolateFE40

Description
    Foam::patchInterpolateFE40

    Interpolates field values from a set of supplied source points in space to supplied target points.
    Two interpolation strategies can be chosen:
    - standard OpenFOAM interpolation (newInterpolation=false, explanation taken from timeVarriyngMappedFixedValue):
        - Points need to be more or less on a plane since get triangulated in 2D.

        -At startup this bc does the triangulation and determines linear
        interpolation (triangle it is in and weights to the 3 vertices)
        for every face centre. Interpolates linearly inbetween times.

        - Switch on debug flag to have it dump the triangulation (in transformed
        space) and transform face centres.
    - Taking values from nearest neighbour point (newInterpolation=true)

SourceFiles
    patchInterpolateFE40.C

\*---------------------------------------------------------------------------*/

#ifndef patchInterpolateFE40_H
#define patchInterpolateFE40_H

#include "FixedList.H"
#include "instantList.H"
#include "coordinateSystem.H"
//#include "AverageIOField.H"
#include "parseTable.h"
#include "nn.h"

#define POINT_SIZE 3

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class patchInterpolateFE40 Declaration
\*---------------------------------------------------------------------------*/

template <class Type> 
class patchInterpolateFE40
{

protected:
    // Protected data
        
        // Informations about source and target mesh
        const pointField& pSource_;
        const pointField& pTarget_;
        label nPointsSource_;
        label nPointsTarget_;
        
        //- use value of nearest neighbour instead of interpolating
        vector maxP_;
        vector minP_;
        double min_[POINT_SIZE];
        double max_[POINT_SIZE];
        double req_[POINT_SIZE];
        double *sample_;
        double *sol_;
        int splits_[POINT_SIZE];
        int dims_;
        label fieldDim_;
        
        //- Fraction of perturbation (fraction of bounding box) to add
        scalar perturb_;

    // Private Member Functions
        // save informations about 
        // source (sample) and target (request) points
        // (needed for nearest neighbour search)
        void setSampleCoord(label pointI);
        void setRequestCoord(label pointI);
        void setSampleValue(label i, label pointI, const Field<Type>& sourceFld );
        void setSolutionField(label i, label pointI, Field<Type>& fld );
        
public:

    //- Runtime type information
        TypeName("patchInterpolateFE40");

    // Constructors

        //- Construct from patch and internal field
        patchInterpolateFE40
        (
            const pointField& pSource,
            const pointField& pTarget,
            const label& fieldDim
        );

        //- Construct as copy
        patchInterpolateFE40
        (
            const patchInterpolateFE40&
        );
    
    // Destructor
        ~patchInterpolateFE40();

    // Member functions
    
        // Interpolation function
        tmp<Field<Type> > interpolate(const Field<Type>& sourceFld);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
