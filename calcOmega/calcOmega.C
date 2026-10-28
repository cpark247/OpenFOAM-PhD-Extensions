/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

Application
    calcOmega

Description
    Computes the specific dissipation rate omega = epsilon/(Cmu*k) (Cmu = 0.09)
    from existing k and epsilon fields and writes it, e.g. to initialise or
    switch a case from a k-epsilon to a k-omega type model.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "graph.H"
#include "OFstream.H"
#include "writeFile.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    //#include "readCalcOmegaDict.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    
    omega = epsilon/(max(k,dimensionedScalar("mink", k.dimensions(), SMALL)) * 0.09);
    omega.write();

    Info<< "end" << endl;

    return 0;
}


// ************************************************************************* //
