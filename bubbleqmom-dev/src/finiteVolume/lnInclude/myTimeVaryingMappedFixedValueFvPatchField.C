/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "myTimeVaryingMappedFixedValueFvPatchField.H"
#include "Time.H"
#include "triSurfaceTools.H"
//#include "triSurface.H"
#include "vector2D.H"
#include "OFstream.H"
//#include "long.H"
//#include "rawIOField.H"

#define NDEBUG

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
myTimeVaryingMappedFixedValueFvPatchField<Type>::
myTimeVaryingMappedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    fieldTableName_(iF.name()),
    sampleTimes_(0),
    startSampleTime_(-1),
    startSampledValues_(0),
    endSampleTime_(-1),
    endSampledValues_(0),
    timeValPeriod_(0.0),
    binaryFile_("test.bin"),
    p_(NULL),
    endTime_(-1.0),
    nTimes_(-1),
    nPoints_(-1),
    nFacePoints_(-1),
    deltaT_(-1.0),
    fieldDim_(-1),
    referenceField_(p.size())
{
#ifdef DEBUG
    Pout<< "myTimeVaryingMappedFixedValue :"
        << " contruct from fvPatch and internalField" << endl;
#endif
}


template<class Type>
myTimeVaryingMappedFixedValueFvPatchField<Type>::
myTimeVaryingMappedFixedValueFvPatchField
(
    const myTimeVaryingMappedFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    fieldTableName_(ptf.fieldTableName_),
    sampleTimes_(0),
    startSampleTime_(-1),
    startSampledValues_(0),
    endSampleTime_(-1),
    endSampledValues_(0),
    timeValPeriod_(ptf.timeValPeriod_),
    binaryFile_(ptf.binaryFile_),
    p_(ptf.p_),
    endTime_(ptf.endTime_),
    nTimes_(ptf.nTimes_),
    nPoints_(ptf.nPoints_),
    nFacePoints_(ptf.nFacePoints_),
    deltaT_(ptf.deltaT_),
    fieldDim_(ptf.fieldDim_),
    referenceField_(mapper(ptf.referenceField_))
{
#ifdef DEBUG
    Pout<< "myTimeVaryingMappedFixedValue"
        << " : construct from mappedFixedValue and mapper" << endl;
#endif
}


template<class Type>
myTimeVaryingMappedFixedValueFvPatchField<Type>::
myTimeVaryingMappedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),
    fieldTableName_(iF.name()),
    sampleTimes_(0),
    startSampleTime_(-1),
    startSampledValues_(0),
    endSampleTime_(-1),
    endSampledValues_(0),
    timeValPeriod_(0.0),
    binaryFile_(dict.lookupOrDefault<fileName>("fileName","test.bin")),
    p_(NULL),
    nTimes_(-1),
    nPoints_(-1),
    nFacePoints_(-1),
    deltaT_(readScalar(dict.lookup("dtBinary"))),
    endTime_(readScalar(dict.lookup("endTimeBinary"))),
    fieldDim_(-1),
    referenceField_(p.size(),Foam::Zero)
    //referenceField_("referenceField", dict, p.size())
{
#ifdef DEBUG
    Pout<< "myTimeVaryingMappedFixedValue : construct from dictionary"
        << endl;
#endif
    // Loading binary for every processor does not lead to memory overhead
    // MMAP is used within flameletConfig to handle memory efficient data lookup
    binaryFile_.expand();
    p_ = pmmap_create_read(binaryFile_.c_str());

    nTimes_ = endTime_/deltaT_ + 0.5;
#ifdef DEBUG
    Info << "deltaT = " << deltaT_ <<  endl;
    Info << "endTime = " << endTime_ <<  endl;
    Info << "nTimes = " <<nTimes_ <<  endl;
#endif
    dict.readIfPresent("fieldTableName", fieldTableName_);
    
    if (dict.found("referenceField"))
    {
       referenceField_ = Field<Type>("referenceField", dict, p.size());
    }
     
    if (dict.found("value"))
    {
        fixedValueFvPatchField<Type>::operator==
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        //fixedValueFvPatchField<Type>::operator==(referenceField_);
        //old implementation:
        updateCoeffs();
    }
}

/*
template<class Type>
myTimeVaryingMappedFixedValueFvPatchField<Type>::
myTimeVaryingMappedFixedValueFvPatchField
(
    const myTimeVaryingMappedFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    fieldTableName_(ptf.fieldTableName_),
    sampleTimes_(ptf.sampleTimes_),
    startSampleTime_(ptf.startSampleTime_),
    startSampledValues_(ptf.startSampledValues_),
    endSampleTime_(ptf.endSampleTime_),
    endSampledValues_(ptf.endSampledValues_),
    timeValPeriod_(ptf.timeValPeriod_),
    binaryFile_(ptf.binaryFile_),
    p_(ptf.p_),
    endTime_(ptf.endTime_),
    nTimes_(ptf.nTimes_),
    nPoints_(ptf.nPoints_),
    nFacePoints_(ptf.nFacePoints_),
    deltaT_(ptf.deltaT_),
    fieldDim_(ptf.fieldDim_),
    referenceField_(ptf.referenceField_)
{
#ifdef DEBUG
    Pout<< "myTimeVaryingMappedFixedValue : copy construct"
        << endl;
#endif
}
*/

template<class Type>
myTimeVaryingMappedFixedValueFvPatchField<Type>::
myTimeVaryingMappedFixedValueFvPatchField
(
    const myTimeVaryingMappedFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    fieldTableName_(ptf.fieldTableName_),
    sampleTimes_(ptf.sampleTimes_),
    startSampleTime_(ptf.startSampleTime_),
    startSampledValues_(ptf.startSampledValues_),
    endSampleTime_(ptf.endSampleTime_),
    endSampledValues_(ptf.endSampledValues_),
    timeValPeriod_(ptf.timeValPeriod_),
    binaryFile_(ptf.binaryFile_),
    p_(ptf.p_),
    endTime_(ptf.endTime_),
    nTimes_(ptf.nTimes_),
    nPoints_(ptf.nPoints_),
    nFacePoints_(ptf.nFacePoints_),
    deltaT_(ptf.deltaT_),
    fieldDim_(ptf.fieldDim_),
    referenceField_(ptf.referenceField_)
{
#ifdef DEBUG
    Pout<< "myTimeVaryingMappedFixedValue :"
        << " copy construct, resetting internal field" << endl;
#endif
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void myTimeVaryingMappedFixedValueFvPatchField<Type>::setFieldDimension()
{
    FatalErrorIn
    (
        "myTimeVaryingMappedFixedValueFvPatchField<Type>::setFieldDimension()"
    )   << "Only implemented for vector and scalar class"
        << exit(FatalError);    
}
template<>
void myTimeVaryingMappedFixedValueFvPatchField<scalar>::setFieldDimension()
{ 
    fieldDim_ = FIELD_DIM_SCALAR;
}
template<>
void myTimeVaryingMappedFixedValueFvPatchField<vector>::setFieldDimension()
{ 
    fieldDim_ = FIELD_DIM_VECTOR;
}

template<class Type>
void myTimeVaryingMappedFixedValueFvPatchField<Type>::setSampledValues
(
    Field<Type>& sampledValues,
    const label & sampleTime
)
{
    FatalErrorIn
    (
        "myTimeVaryingMappedFixedValueFvPatchField<Type>::setFieldDimension()"
    )   << "Only implemented for vector and scalar class"
        << exit(FatalError);    
}
template<>
void myTimeVaryingMappedFixedValueFvPatchField<scalar>::setSampledValues
(
    Field<scalar>& sampledValues,
    const label & sampleTime
)
{
    int nFacePoints_=this->patch().patch().faceCentres().size();
    tmp<Field<scalar> > tfld(new Field<scalar>(nFacePoints_));
    Field<scalar>& fld = tfld.ref();

    // Create map for Source binary
    double (*map)[nTimes_][nFacePoints_];
    map=(double (*)[nTimes_][nFacePoints_]) pmmap_get_vptr(p_);

#ifdef DEBUG
        Info << ">>> Index in binary: " << sampleTime << " <<<" << endl;
#endif
        
    //Info << "Size of vals is " << vals.size() << endl; 
    for(int i=0;i<nFacePoints_;i++)
    { 
        fld[i]=(*map)[sampleTime][i];
    }
    
    // spatial interpolation at runt-time could be called here instead
    sampledValues = tfld;
}

template<>
void myTimeVaryingMappedFixedValueFvPatchField<vector>::setSampledValues
(
    Field<vector>& sampledValues,
    const label & sampleTime
)
{
    int nFacePoints_=this->patch().patch().faceCentres().size();
    tmp<Field<vector> > tfld(new Field<vector>(nFacePoints_));
    Field<vector>& fld = tfld.ref();
    
    // Create map for Source binary
    double (*map)[nTimes_][nFacePoints_][fieldDim_];
    map=(double (*)[nTimes_][nFacePoints_][fieldDim_]) pmmap_get_vptr(p_);
      
#ifdef DEBUG
        Info << ">>> Index in binary: " << sampleTime << " <<<" << endl;
#endif

    //Info << "Size of vals is " << vals_.size() << endl; 
    for(int i=0;i<nFacePoints_;i++)
    { 
        for (int j=0; j<fieldDim_; j++)
        {
            fld[i].component(j)=(*map)[sampleTime][i][j];
        }
    }
    
    // spatial interpolation at runt-time could be called here instead
    sampledValues = tfld;
}

template<class Type>
void myTimeVaryingMappedFixedValueFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
    m(referenceField_, referenceField_);
    if (startSampledValues_.size() > 0)
    {
        m(startSampledValues_, startSampledValues_);
        m(endSampledValues_, endSampledValues_);
    }
}


template<class Type>
void myTimeVaryingMappedFixedValueFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

    const myTimeVaryingMappedFixedValueFvPatchField<Type>& tiptf =
        refCast<const myTimeVaryingMappedFixedValueFvPatchField<Type> >(ptf);

    referenceField_.rmap(tiptf.referenceField_, addr);

    startSampledValues_.rmap(tiptf.startSampledValues_, addr);
    endSampledValues_.rmap(tiptf.endSampledValues_, addr);
}


template<class Type>
void myTimeVaryingMappedFixedValueFvPatchField<Type>::readSamplePoints()
{
    // Read the sample points
    Info << "Trying to read sample points" << endl;

    pointIOField samplePoints
    (
        IOobject
        (
            "points",
	    this->db().time().path().path()/this->db().time().constant()/"boundaryData"/this->patch().name(),
            "boundaryData"/this->patch().name(),
            this->db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            false
        )
    );
    const fileName samplePointsFile = samplePoints.filePath();

//   const Time& time = this->db().time();
//// Reread values and interpolate
//   const fileName samplePointsFile
//   (
//       time.globalPath()
//      /time.constant()
//      /"boundaryData"
//      /this->patch().name()
//      /"points"
//   );
// 
//   IOobject io
//   (
//       samplePointsFile,   // absolute path
//       time,
//       IOobject::MUST_READ,
//       IOobject::NO_WRITE,
//       false,              // no need to register
//       true                // is global object (currently not used)
//   );
//  
//    // Read data
//    const rawIOField<point> samplePoints(io, false);

#ifdef DEBUG
    Info<< "myTimeVaryingMappedFixedValueFvPatchField :"
        << " Read " << samplePoints.size() << " sample points from "
        << samplePointsFile << endl;
#endif
    nPoints_ = samplePoints.size();

    // Read the times for which data is available

    const fileName samplePointsDir = samplePointsFile.path();

#ifdef DEBUG
    Info << "nTimes = " << nTimes_ << endl;
    Info << "endTime = " << endTime_ <<  endl;
    Info << "deltaT = " << deltaT_ << endl;
#endif
    sampleTimes_.setSize(nTimes_);
    forAll(sampleTimes_,i)
    {
      //sampleTimes_[i]=instant(((i+1)*deltaT_));
      // Make sure, that there is an entry for t=0
      sampleTimes_[i]=instant(i*deltaT_);
    }

#ifdef DEBUG
    Info<< "myTimeVaryingMappedFixedValueFvPatchField : In directory "
            << samplePointsDir << " found times " << timeNames(sampleTimes_)
            << endl;
#endif

}

template<class Type>
wordList myTimeVaryingMappedFixedValueFvPatchField<Type>::timeNames
(
    const instantList& times
)
{
    wordList names(times.size());

    forAll(times, i)
    {
        names[i] = times[i].name();
    }
    return names;
}


template<class Type>
void myTimeVaryingMappedFixedValueFvPatchField<Type>::findTime
(
    const fileName& instance,
    const fileName& local,
    const scalar timeVal,
    label& lo,
    label& hi
)
{
    // Binary is read from the beginning if its end time is reached
    timeValPeriod_ = timeVal;
    lo = startSampleTime_;
    hi = -1;

    int full_periods = (int) (timeVal / (endTime_ - deltaT_));

#ifdef DEBUG
    Info << "sampleTime: " << startSampleTime_ << "/" << sampleTimes_.size() << endl;
    Info << "Full periods: " << full_periods << endl;
    Info << "Time value (period): " << timeVal - full_periods *(endTime_ - deltaT_) << endl; 
#endif

    //if (timeVal > (endTime_ - deltaT_))  
    if (timeVal > endTime_)  
    {
        //shift sampleTime to the first available binary time period
        //int full_periods = (int) (timeVal / (endTime_ - deltaT_)); 
        timeValPeriod_ = timeVal - full_periods * (endTime_ - deltaT_);
#ifdef DEBUG
    Pout<< "myTimeVaryingMappedFixedValueFvPatchField<Type>::findTime" << endl
        << " Current Simulation time exceeds end-time of binary File (" << (endTime_ - deltaT_) << ")" << ", number of period: " << full_periods << endl
        << " Set time for request from "<< timeVal << " to " << timeValPeriod_ << endl;
#endif
    } 

    // Estimate startSampleTime 
    startSampleTime_ = (int) (timeValPeriod_/deltaT_) - 1;   
#ifdef DEBUG
    Pout<< "myTimeVaryingMappedFixedValueFvPatchField<Type>::findTime" << endl
        << " Current startSampleTime_ estimate: "<< startSampleTime_ <<  endl;
#endif

    for (label i = startSampleTime_+1; i < sampleTimes_.size(); i++)
    {
        if (sampleTimes_[i].value() > timeValPeriod_)
        {
            break;
        }
        else
        {
            lo = i;
        }
    }

    if (lo == -1)
    {
        FatalErrorIn("findTime")
            << "Cannot find starting sampling values for current time "
            << timeValPeriod_ << nl
            << "Have sampling values for times "
            << timeNames(sampleTimes_) << nl
            << "In directory "
            <<  this->db().time().constant()/"boundaryData"/this->patch().name()
            << "\n    on patch " << this->patch().name()
            << " of field " << fieldTableName_
            << exit(FatalError);
    }

    if (lo < sampleTimes_.size()-1)
    {
        hi = lo+1;
    }


#ifdef DEBUG
    if (hi == -1)
    {
        Pout<< "findTime : Found time " << timeVal << " (" << timeValPeriod_ << ")" << " after"
            << " index:" << lo << " time:" << sampleTimes_[lo].value()
            << endl;
    }
    else
    {
        Pout<< "findTime : Found time " << timeVal << " (" << timeValPeriod_ << ")" << " inbetween"
            << " index:" << lo << " time:" << sampleTimes_[lo].value()
            << " and index:" << hi << " time:" << sampleTimes_[hi].value()
            << endl;
    }
#endif
}


template<class Type>
void myTimeVaryingMappedFixedValueFvPatchField<Type>::checkTable()
{
    // Initialise
    // Now done in constructor of patchInterpolate class
    if (startSampleTime_ == -1 && endSampleTime_ == -1)
    {
        readSamplePoints();
    }

    setFieldDimension();
    nFacePoints_ = this->patch().patch().faceCentres().size();

//TODO: Check needs to be adapted (not only on master)
//      Check on Processor if patchField is not empty
//      Size must match samplePoinst
//    if(noInterpolation_)
//    {
//        //Info << "Setting Interpolation to false" << endl;
//        //useNewInterpolation_=false;
//        if((Pstream::master())&&(nPoints_!=nFacePoints_))
//        {
//
//            FatalErrorIn("myTimeVaryingMappedFixedValueFE40::myTimeVaryingMappedFixedValueFE40")
//            << "Size of samplePoints is not equal patch size!"
//            << endl << abort(FatalError);
//        }
//    }
    
    // Find current time in sampleTimes
    label lo = -1;
    label hi = -1;

    findTime
    (
        this->db().time().constant(),
        "boundaryData"/this->patch().name(),
        this->db().time().value(),
        lo,
        hi
    );

    // Update sampled data fields.
#ifdef DEBUG
    Info << "lo: " << lo << endl; 
    Info << "hi: " << hi << endl; 
#endif 

    if (lo != startSampleTime_)
    {
        startSampleTime_ = lo;

        if (startSampleTime_ == endSampleTime_)
        {
            // No need to reread since are end values
 
#ifdef DEBUG
            Pout<< "checkTable : Setting startValues to (already read) "
                <<   "boundaryData"
                    /this->patch().name()
                    /sampleTimes_[startSampleTime_].name()
                << endl;
#endif            
            startSampledValues_ = endSampledValues_;
        }
        else
        {
#ifdef DEBUG
            Pout<< "checkTable : Reading startValues from "
                <<   "boundaryData"
                    /this->patch().name()
                    /sampleTimes_[lo].name()
                << endl;
#endif
            // set sampled values to field startSampledValues_ 
            setSampledValues
            (
                startSampledValues_,
                startSampleTime_
            );
        }
    }

    if (hi != endSampleTime_)
    {
        endSampleTime_ = hi;

        if (endSampleTime_ == -1)
        {
            // endTime no longer valid. Might as well clear endValues.
#ifdef DEBUG                
            Pout<< "checkTable : Clearing endValues" << endl;
#endif            
            endSampledValues_.clear();
        }
        else
        {
#ifdef DEBUG                
            Pout<< "checkTable : Reading endValues from "
                <<   "boundaryData"
                    /this->patch().name()
                    /sampleTimes_[endSampleTime_].name()
                << endl;
#endif        
        
            // set sampled values to field endSampledValues_
            setSampledValues
            (
                endSampledValues_,
                endSampleTime_
            );
        }
    }
}


template<class Type>
void myTimeVaryingMappedFixedValueFvPatchField<Type>::updateCoeffs()
{
    //TODO check if we get a performance incease when updating BC only on the corresponding processor
    //if (this->updated() || this->patch().size() == vals_.size())
    if (this->updated())
    {
        return;
    }

    checkTable();

    // Interpolate between the sampled data

    if (endSampleTime_ == -1)
    {
        // only start value
#ifdef DEBUG
        Pout<< "updateCoeffs : Sampled, non-interpolated values"
            << " from start time:"
            << sampleTimes_[startSampleTime_].name() << nl;
#endif
        this->operator==(startSampledValues_);
    }
    else
    {
        scalar start = sampleTimes_[startSampleTime_].value();
        scalar end = sampleTimes_[endSampleTime_].value();

        scalar s = (timeValPeriod_-start)/(end-start);

#ifdef DEBUG
        Pout<< "updateCoeffs : Sampled, interpolated values"
            << " between start time:"
            << sampleTimes_[startSampleTime_].name()
            << " and end time:" << sampleTimes_[endSampleTime_].name()
            << " with weight:" << s << endl;
#endif
        this->operator==((1-s)*startSampledValues_ + s*endSampledValues_);
    }

    // add reference Field
    const Field<Type>& fld = *this;
    this->operator==(fld+referenceField_);
    
    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void myTimeVaryingMappedFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("fileName") << binaryFile_ << token::END_STATEMENT << nl;
    os.writeKeyword("dtBinary") << deltaT_ << token::END_STATEMENT << nl;
    os.writeKeyword("endTimeBinary") << endTime_ << token::END_STATEMENT << nl;


    if (fieldTableName_ != this->internalField().name())
    {
        os.writeKeyword("fieldTableName") << fieldTableName_ << token::END_STATEMENT << nl;
    }
    writeEntry(os, "referenceField", referenceField_);
    writeEntry(os, "value", static_cast<const Field<Type>&>(*this));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
