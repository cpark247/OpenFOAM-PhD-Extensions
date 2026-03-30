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
    boxTurb

Description
    Puts together resolved and modelled portion of RST to give the total RST

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


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // order of R components: Rxx Rxy Rxz Ryy Ryz Rzz
    // R_ = symmTensor(uPrime_uPrime_Mean_+2/3*kMod_, uPrime_vPrime_Mean_, uPrime_wPrime_Mean_, vPrime_vPrime_Mean_+2/3*kMod_, vPrime_wPrime_Mean_, wPrime_wPrime_Mean_+2/3*kMod_);
    
    // Order of R components: Rxx Rxy Rxz Ryy Ryz Rzz

    forAll(R_, cell)
    {
        R_[cell] = symmTensor(
            uPrime_uPrime_Mean_[cell] + 2.0/3.0 * kMod_[cell],
            uPrime_vPrime_Mean_[cell],
            uPrime_wPrime_Mean_[cell],
            vPrime_vPrime_Mean_[cell] + 2.0/3.0 * kMod_[cell],
            vPrime_wPrime_Mean_[cell],
            wPrime_wPrime_Mean_[cell] + 2.0/3.0 * kMod_[cell]
        );
    }

    R_.write();

    Info<< "end of total RST calculation" << endl;

    return 0;
}


// ************************************************************************* //
