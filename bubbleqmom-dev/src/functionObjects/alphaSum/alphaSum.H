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

Class
    Foam::functionObjects::alphaSum

Description
    This function object calculates the volume of a specified phase inside a specified cellSet. cellSet has to be named evalAlpha.

    Example of function object specification:
    \verbatim
    alphaSum1
    {
        type           alphaSum;
        libs ("libalphaSumFunctionObject.so");
        phase		gas;
    }
    \endverbatim

Usage
    \table
        Property     | Description               | Required | Default value
        type         | type name: alphaSum | yes      |
        phase   | name of phase      | yes      |
    \endtable

SourceFiles
    alphaSum.C

\*---------------------------------------------------------------------------*/

#ifndef alphaSum_H
#define alphaSum_H

#include "fvMeshFunctionObject.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class alphaSum Declaration
\*---------------------------------------------------------------------------*/

class alphaSum
:
    public fvMeshFunctionObject
{
    // Private data

	OFstream OS_;
	labelList cells_;
	word cellSetName_;
	word alphaName_;
	scalar gasVol_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        alphaSum(const alphaSum&);

        //- Disallow default bitwise assignment
        void operator=(const alphaSum&);


public:

    //- Runtime type information
    TypeName("alphaSum");


    // Constructors

        //- Construct from Time and dictionary
        alphaSum
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~alphaSum();


    // Member Functions

        //- Read the alphaSum data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Write the alphaSum
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
