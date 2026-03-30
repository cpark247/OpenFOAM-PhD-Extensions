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

\*---------------------------------------------------------------------------*/

#include "patchInterpolateFE40.H"
//#include "foamTime.H"
#include "triSurfaceTools.H"
//#include "triSurface.H"
#include "vector2D.H"
#include "OFstream.H"
//#include "long.H"
// * * * *  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// force instantiation of class with type scalar and vector
// (otherwise some functions are maybe not defined 
// when called by executable) 
template class patchInterpolateFE40<vector>;
template class patchInterpolateFE40<scalar>;

defineTemplateTypeNameAndDebug(patchInterpolateFE40<vector>, 0);
defineTemplateTypeNameAndDebug(patchInterpolateFE40<scalar>, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type> 
void patchInterpolateFE40<Type>::setSampleCoord(label pointI)
{
    if (debug)
    {
        printf("\t Set sample coordinate for point %i out of %i\n",pointI,nPointsSource_);
    }
    
    for (int i=0; i<POINT_SIZE; i++)
    {
        sample_[i] = pSource_[pointI][i];
    }
    if (debug)
    {
        Info << "Sample Coordinates of point " << pointI << ": ( "
        << sample_[0] << " " << sample_[1] << " " << sample_[2] << " )"  << endl;
    }
}

template<class Type> 
void patchInterpolateFE40<Type>::setRequestCoord(label pointI)
{
    if (debug)
    {
        printf("\t Set request coordinate for point %i out of %i\n",pointI,nPointsTarget_);
    }
    
    for (int i=0; i<POINT_SIZE; i++)
    {
        req_[i] = pTarget_[pointI][i];
    }

    if (debug)
    {
        Info << "Request Coordinates of point " << pointI << ": ( "
        << req_[0] << " " << req_[1] << " " << req_[2] << " )"  << endl;
    }
}

template<class Type> 
void patchInterpolateFE40<Type>::setSampleValue(label i, label pointI, const Field<Type>& sourceFld )
{
    FatalErrorIn
    (
        "patchInterpolateFE40<Type>::setSampleValues"
    )   << "Only implemented for vector and scalar class"
        << exit(FatalError);
}

template<> 
void patchInterpolateFE40<scalar>::setSampleValue(label i, label pointI, const Field<scalar>& sourceFld )
{
    sample_[POINT_SIZE+i] = sourceFld[pointI]; 
}

template<> 
void patchInterpolateFE40<vector>::setSampleValue(label i, label pointI, const Field<vector>& sourceFld )
{
    sample_[POINT_SIZE+i] = sourceFld[pointI].component(i); 
}

template<class Type> 
void patchInterpolateFE40<Type>::setSolutionField(label i, label pointI, Field<Type>& sourceFld )
{
    FatalErrorIn
    (
        "patchInterpolateFE40<Type>::setSolutionField"
    )   << "Only implemented for vector and scalar class"
        << exit(FatalError);
}

template<> 
void patchInterpolateFE40<scalar>::setSolutionField(label i, label pointI, Field<scalar>& fld )
{
    fld[pointI] = sol_[0]; 
}

template<> 
void patchInterpolateFE40<vector>::setSolutionField(label i, label pointI, Field<vector>& fld )
{
    fld[pointI].component(i) = sol_[i]; 
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
patchInterpolateFE40<Type>::patchInterpolateFE40
(
    const pointField& pSource,
    const pointField& pTarget, 
    const label& fieldDim
)
:
    pSource_(pSource),
    pTarget_(pTarget),
    sample_(NULL),
    sol_(NULL),
    nPointsSource_(-1),
    nPointsTarget_(-1),
    fieldDim_(fieldDim),
    dims_(POINT_SIZE)
{    
    if (debug)
    { 
        Info << "Field has dimension: " << fieldDim_ << endl;
    }
    
    nPointsSource_ = pSource_.size();
    nPointsTarget_ = pTarget_.size();
    
    if(nPointsSource_ > 0 && nPointsTarget_ > 0) 
    { 
        // determine parameters for nn-Search
        
        // Allocate memory for sample and solution pointer
        sample_ = new double[POINT_SIZE+fieldDim_];
        sol_ = new double[fieldDim_];
        
        maxP_ = max(pSource_);
        minP_ = min(pSource_);
        forAll(maxP_,i)
        {
            max_[i]=maxP_[i];
            min_[i]=minP_[i];
            if(min_[i]==max_[i])
            {
                splits_[i]=0;
            //    max[i]+=1.0;
            }
            else
            {
                splits_[i]=0;
            }
        }
    }
    else if (nPointsSource_ <= 0)
    {
        FatalErrorIn
        (
            "patchInterpolateFE40::patchInterpolatFE40"
            "(const pointField&, const pointField&, const bool&, const label&, const label& )"
        )   << "No points in source given." << nl
            << "Hence, no interpolation possible"
            << exit(FatalError);
    }
    else if (nPointsTarget_ <= 0)
    {
        FatalErrorIn
        (
            "patchInterpolateFE40::patchInterpolatFE40"
            "(const pointField&, const pointField&, const bool&, const label&, const label& )"
        )   << "No points in target given." << nl
            << "Hence, no interpolation possible"
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class Type>
patchInterpolateFE40<Type>::~patchInterpolateFE40()
{
    if(sample_!=NULL)
    {
        delete(sample_);
        sample_ = NULL;
    }
    if(sol_!=NULL)
    {
        delete(sol_);
        sol_ = NULL;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class Type>
tmp<Field<Type> > patchInterpolateFE40<Type>::interpolate
(
    const Field<Type>& sourceFld
)
{
    // Create tmp field (tfld)
    // and set reference to this (fld)
    tmp<Field<Type> > tfld(new Field<Type>(pTarget_.size()));
    Field<Type>& fld = tfld.ref();

    // Creating new interpolation object for current time step
    nn *n = nn_create(dims_, fieldDim_, min_, max_, splits_);
    
    if (debug)
    {
        Info << "\tCreating samples for newInterpolation" << endl;
    }
    
    
    forAll(pSource_,pointI)
    {
        setSampleCoord(pointI);
        if (debug)
        {
            Info << "Field value of sample point " << pointI << ": ( ";
        }
        for(int i=0; i<fieldDim_; i++)
        {
            setSampleValue(i, pointI, sourceFld);
            if (debug)
            {
                Info << sample_[POINT_SIZE+i];
            }
        }
        if (debug)
        {
            Info << " )" << endl;
        }
        // append data to nearest neighbour object
        nn_addData(n, sample_);
    }
    
    if (debug)
    {
        Info << "\tCreating request and perform newInterpolation" << endl;
    }
    
    // Loop over all faces of target Patch
    forAll(pTarget_,pointI)
    {
        // Set request coordinate
        setRequestCoord(pointI);    
        
        // Determine field value of 
        // nearest neighbour in source mesh
        nn_getData(n, req_, sol_);
        
        // write field value to solution field
        forAll(minP_,coordI)
        {
            setSolutionField(coordI,pointI,fld);
        }
    }
    
    // free memory
    nn_destroy(n);
    
    if (debug)
    {
        Info << "Interpolated field:" << endl;
        Info << fld << endl;
    }
    
    //ToDo: Check performance for large fields
    return tfld;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
