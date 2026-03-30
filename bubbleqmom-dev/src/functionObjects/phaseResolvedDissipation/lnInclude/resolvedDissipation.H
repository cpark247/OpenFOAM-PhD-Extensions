/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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
    Foam::functionObjects::pseudoDissipationTotal

Description
    This function object calculates effective energy dissipation rate and performs
    averaging in time. The mean field is written into time folder according to
    output time settings. 

    At the moment energy dissipation rate function object can be used for 
    incompressible flows only. The determination of disispation rate for 
    compressible flows has to be implemented yet. 

    Example of function object specification:
    \verbatim
    pseudoDissipationTotal1
    {
        type           pseudoDissipationTotal;
        libs ("libpseudoDissipationTotalFunctionObject.so");
        ...
        fieldName      epsilon;
    }
    \endverbatim

Usage
    \table
        Property     | Description                          | Required | Default value
        type         | type name: pseudoDissipationTotal    | yes      |
        fieldName    | name of the energy dissipation field | yes      |
    \endtable

SourceFiles
    pseudoDissipationTotal.C

\*---------------------------------------------------------------------------*/

#ifndef resolvedDissipation_H
#define resolvedDissipation_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class resolvedDissipation Declaration
\*---------------------------------------------------------------------------*/

class resolvedDissipation
:
    public fvMeshFunctionObject
{
public:

        static const word turbulenceModelName;

protected:
    // Protected data

        //- Name of this set of resolvedDissipation object
        word fieldName_;
	
	//- Name of the current phase
	word phaseName_;

        //- on/off switch
        bool active_;

    // Protected Member Functions

        //- Disallow default bitwise copy construct
        resolvedDissipation(const resolvedDissipation&);

        //- Disallow default bitwise assignment
        void operator=(const resolvedDissipation&);

        //- Return true if compressible turbulence model is identified
        bool compressible();


public:

    //- Runtime type information
    TypeName("phaseResolvedDissipation");


    // Constructors

        //- Construct from Time and dictionary
        resolvedDissipation
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~resolvedDissipation();


    // Member Functions

        //- Read the resolvedDissipation data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Write the resolvedDissipation
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
