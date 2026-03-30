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

Class
    Foam::myTimeVaryingMappedFixedValueFvPatchField

Description
    Foam::myTimeVaryingMappedFixedValueFvPatchField
    Example usage in U-field:
        boundaryField
        {
            Fuel
            {
                type            myTimeVaryingMappedFixedValue;
                setAverage      0;
                perturb         0;
                fileName        "$FOAM_CASE/constant/boundaryData/Fuel/dnsIgnition_INLET_dt_5.575e-08_endTime_5.715e-05_reduced_U.bin";
                dtBinary        5.575e-08;
                endTimeBinary   5.715e-05;
                readBinary      1;
                useNewInterpolation 0;
                noInterpolation 1;
                value           uniform ( 0 0 28 );
                referenceField  nonuniform List<vector> 2 ( (0 1 2) (0 1 2));
            }
        }

SourceFiles
    myTimeVaryingMappedFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef myTimeVaryingMappedFixedValueFvPatchField_H
#define myTimeVaryingMappedFixedValueFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "coordinateSystem.H"
#include "FixedList.H"
#include "instantList.H"
#include "pmmap.h"
//#include "AverageIOField.H"
#include "parseTable.h"

#define FIELD_DIM_SCALAR 1
#define FIELD_DIM_VECTOR 3


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class myTimeVaryingMappedFixedValueFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class myTimeVaryingMappedFixedValueFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private data

        //- Name of the field data table, defaults to the name of the field
        word fieldTableName_;

        //- List of boundaryData time directories
        instantList sampleTimes_;

        //- Current starting index in sampleTimes
        label startSampleTime_;

        //- Interpolated values from startSampleTime
        Field<Type> startSampledValues_;

        //- Current end index in sampleTimes
        label endSampleTime_;

        //- Interpolated values from endSampleTime
        Field<Type> endSampledValues_;

        //- Current time value shifted to lay in the time range of binary file
        // In this way, the binary file is read from the beginning again, when current simulation 
        // time is greater then end Time of binary
        // Set in findTime()
        scalar timeValPeriod_;

        //- new for binary read
        fileName binaryFile_;
        pmmap *p_;
        // End Time of the binary file, read from dictionary
        scalar endTime_;
        label nTimes_;
//        // Number of sample points
        label nPoints_;
//        // Number of points on Face
//        // (might differ from nPoints_ if patch is splitted 
//        // or wrong sample points file is given)
        label nFacePoints_;
        // Delta Time with which binary file was written, read from dictionary
        scalar deltaT_;
//
        //- Dimension of Field Types
        label fieldDim_;

        //- Reference field, to which the read field is added
        Field<Type> referenceField_;


    // Private Member Functions

        //- Get names of times
        static wordList timeNames(const instantList&);

        //- Find times around current time
        void findTime
        (
            const fileName& instance,
            const fileName& local,
            const scalar timeVal,
            label& lo,
            label& hi
        );

        //- Read boundary points and determine interpolation weights to patch
        //  faceCentres
        void readSamplePoints();

        //- Find boundary data inbetween current time and interpolate
        void checkTable();

        //- Do actual interpolation using current weights
        tmp<Field<Type> > interpolate(const Field<Type>&) const;

        //- Set the dimension of the field size
        void setFieldDimension();

        //- Set sampled values for timesteps before and after current time
        void setSampledValues
        (
            Field<Type>& sampledValues,
            const label & sampleTime
        );

public:

    //- Runtime type information
    TypeName("myTimeVaryingMappedFixedValue");


    // Constructors

        //- Construct from patch and internal field
        myTimeVaryingMappedFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        myTimeVaryingMappedFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given myTimeVaryingMappedFixedValueFvPatchField
        //  onto a new patch
        myTimeVaryingMappedFixedValueFvPatchField
        (
            const myTimeVaryingMappedFixedValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );


        //- Disallow copy without setting internal field reference
        myTimeVaryingMappedFixedValueFvPatchField
        (
            const myTimeVaryingMappedFixedValueFvPatchField<Type>&
        ) = delete;

        //- Construct and return a clone
        //virtual tmp<fvPatchField<Type> > clone() const
        //{
        //    return tmp<fvPatchField<Type> >
        //    (
        //        new myTimeVaryingMappedFixedValueFvPatchField<Type>(*this)
        //    );
        //}

        //- Construct as copy setting internal field reference
        myTimeVaryingMappedFixedValueFvPatchField
        (
            const myTimeVaryingMappedFixedValueFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new myTimeVaryingMappedFixedValueFvPatchField<Type>(*this, iF)
            );
        }
 

    // Member functions

        // Access

            //- Return the coordinateSystem
//            const coordinateSystem& referenceCS() const
//            {
//                return referenceCS_;
//            }
//
            //- Return the reference field
            const Field<Type>& referenceField() const
            {                                 
                return referenceField_;
            }

            //- Return reference to the reference field to allow adjustment
            Field<Type>& referenceField()
            {
                return referenceField_;
            }

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();
        
        // Utility functions



        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "myTimeVaryingMappedFixedValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
