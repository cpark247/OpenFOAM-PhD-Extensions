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

Class
    Foam::momentsTruncFilter

Description
    Moments truncation filter for LES.

	This method, which is based on a truncation of the moments of 
	the filter, consists in expanding any variable to be filtered
	into Taylor series before performing the convolution.

    \verbatim
    Kernel                 as filter          
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Box filter:            secondMCoeff = 1/24  
    Spherical box filter:  secondMCoeff = 1/64  
    Gaussian filter:       secondMCoeff = 1/24  
    \endverbatim

SourceFiles
    momentsTruncFilter.C

\*---------------------------------------------------------------------------*/

#ifndef momentsTruncFilter_H
#define momentsTruncFilter_H

#include "LESfilterNSteps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class momentsTruncFilter Declaration
\*---------------------------------------------------------------------------*/

class momentsTruncFilter
:
    public LESfilterNSteps
{
    // Private data

        scalar secondMCoeff_;

    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        momentsTruncFilter(const momentsTruncFilter&);
        void operator=(const momentsTruncFilter&);


public:

    //- Runtime type information
    TypeName("momentsTrunc");


    // Constructors

        //- Construct from components
        momentsTruncFilter(const fvMesh& mesh, scalar secondMCoeff);

        //- Construct from IOdictionary
        momentsTruncFilter(const fvMesh& mesh, const dictionary&);


    //- Destructor
    virtual ~momentsTruncFilter()
    {}


    // Member Functions

        //- Read the LESfilterNSteps dictionary
        virtual void read(const dictionary&);


    // Member Operators

        virtual tmp<volScalarField> operator()
        (
            const tmp<volScalarField>& unFilteredField,
            const tmp<volScalarField>& filterDelta
        );

        virtual tmp<volVectorField> operator()
        (
            const tmp<volVectorField>& unFilteredField,
            const tmp<volScalarField>& filterDelta
        );

        virtual tmp<volScalarField> operator()
        (
            const tmp<volScalarField>& unFilteredField,
            const tmp<volScalarField>& filterDelta,
            const label n
        );

        virtual tmp<volVectorField> operator()
        (
            const tmp<volVectorField>& unFilteredField,
            const tmp<volScalarField>& filterDelta,
            const label n
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
