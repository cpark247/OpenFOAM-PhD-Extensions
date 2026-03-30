/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2021 OpenFOAM Foundation
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
    Foam::phaseSystemHaensch

Description
    Class to represent a system of phases and model interfacial transfers
    between them. Uses the equation to track the interface ttracking technique according to Haensch (2012) "A multi-field two-fluid concept for transitions between different scales of interfacial structures"

SourceFiles
    phaseSystemHaensch.C

\*---------------------------------------------------------------------------*/

#ifndef phaseSystemHaensch_H
#define phaseSystemHaensch_H

#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class phaseSystemHaensch Declaration
\*---------------------------------------------------------------------------*/

class phaseSystemHaensch
:
    public phaseSystem
{
protected:

    // Protected typedefs

	typedef HashPtrTable<volScalarField, word, word::hash>
	    gammaTable;

	typedef HashTable<scalar, word, word::hash> haensch_parameter;


    // Protected data

	//- factor a_B see Eq.(8) in Haensch
	scalar haensch_ab_;

	//- boundary \alpha_{clust,min} see Eq.(8) in Haensch
	haensch_parameter haensch_alphaMin_;

	//- boundary \alpha_{clust,max} see Eq.(8) in Haensch
	haensch_parameter haensch_alphaMax_;

	//- \phi_{clust} see Eq.(8) in Haensch
	gammaTable gamma_;

public:

    //- Runtime type information
    TypeName("phaseSystemHaensch");

    // Constructors
        //- Construct from fvMesh
        phaseSystemHaensch(const fvMesh& mesh);

    //- Destructor
    virtual ~phaseSystemHaensch();

    //- Member functions
	//- Solve for the phase fractions
	virtual void solve
	(
	    const PtrList<volScalarField>& rAUs,
	    const PtrList<surfaceScalarField>& rAUfs
	);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
