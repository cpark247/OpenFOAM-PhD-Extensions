/*---------------------------------------------------------------------------*\
Date: June 24, 2005
Author: Eugene de Villiers
Source: http://www.cfd-online.com/Forums/openfoam-solving/58043-les-2.html#post187619
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Application
    perturbU

Description
    initialise channel velocity with superimposed streamwise streaks.
    To be used to force channelOodles to transition and reach a fully
    developed flow sooner.

    Reads in perturbUDict.

    EdV from paper:
        Schoppa, W. and Hussain, F.
        "Coherent structure dynamics in near wall turbulence",
        Fluid Dynamics Research, Vol 26, pp119-139, 2000.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

    // These values should be edited according to the case
    // In all cases, the bottom wall should coincide with y=0
    // h is the half-height of the channel
    // z is spanwise
    // x streamwise
    // Ubar should be set in transportProperties

    IOdictionary perturbDict
    (
        IOobject
        (
            "perturbUDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    const scalar Retau(readScalar(perturbDict.lookup("Retau")));
    //const scalar h(readScalar(perturbDict.lookup("h")));
    scalar a(readScalar(perturbDict.lookup("a"))); // y-axis length
    scalar b(readScalar(perturbDict.lookup("b"))); // z-axis length
    const bool setBulk(readBool(perturbDict.lookup("setBulk")));
    const bool perturb(readBool(perturbDict.lookup("perturb")));
    const direction streamDir(readLabel(perturbDict.lookup("streamwise")));
    const direction spanDir(readLabel(perturbDict.lookup("spanwise")));

    Info
        //<< "Duct half height       = " << h << nl
        << "Re(tau)                   = " << Retau << nl
        << "Set bulk flow             = " << Switch(setBulk) << nl
        << "Perturb flow              = " << Switch(perturb) << nl
        << "Streamwise flow component = " << streamDir << nl
        << "Spanwise flow component   = " << spanDir << nl
        << endl;


    if (!setBulk && !perturb)
    {
        FatalErrorIn(args.executable())
            << "At least one of setBulk or perturb needs to be set"
            << " to do anything to the velocity"
            << exit(FatalError);
    }

    if (streamDir > 2 || spanDir > 2 || streamDir == spanDir)
    {
        FatalErrorIn(args.executable())
            << "Spanwise and streamwise components have to be 0,1 or 2 and"
            << " differ from one another." << nl
            << "streamDir:" << streamDir
            << " spanDir:" << spanDir
            << exit(FatalError);
    }

    // Get component normal to streamDir and spanDir. This is the height
    // component.
    direction heightDir = 0;
    if (streamDir == heightDir)
    {
        heightDir++;
    }
    if (spanDir == heightDir)
    {
        heightDir++;
    }


    IOobject Uheader
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );
    Info<< "Reading U" << endl;
    volVectorField U(Uheader, mesh);

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );
    dimensionedVector Ubar
    (
        transportProperties.lookup("Ubar")
    );


    Info<< "nu      = " << nu << endl;
    Info<< "Ubar    = " << Ubar << endl;
    Info<< "Re(tau) = " << Retau << endl;
    //const scalar utau = Retau*nu.value()/h;
    const scalar utau = Retau*nu.value()/((a+b)/2);
    Info<< "u(tau)  = " << utau << endl;


    //wall normal circulation
    const scalar duplus = Ubar.value()[streamDir]*0.25/utau;
    //spanwise wavenumber: spacing z+ = 200
    //const scalar betaPlus = 2.0*constant::mathematical::pi*(1.0/200.0);
    const scalar sigma = 0.00055;
    //streamwise wave number: spacing x+ = 500
    //const scalar alphaPlus = 2.0*constant::mathematical::pi*(1.0/500.0);
    //const scalar epsilon = Ubar.value()[streamDir]/200.0;

    // Random number generator
    Random perturbation(1234567);

    const vectorField& centres = mesh.C();

    forAll(centres, celli)
    {
        // add a small (+/-20%) random component to enhance symetry breaking
        scalar deviation=1.0 + 0.2*perturbation.scalarNormal();

        const vector& cCentre = centres[celli];

        ////scalar zplus = cCentre[spanDir]*Retau/h;
        //scalar y = min(cCentre[heightDir], 2*h-cCentre[heightDir]); // when y = [0, 2h]
        //scalar z = min(cCentre[spanDir], 2*h-cCentre[spanDir]); // when z = [0, 2h]
        scalar y = min(cCentre[heightDir], 2*a-cCentre[heightDir]); // when y = [0, 2h]
        scalar z = min(cCentre[spanDir], 2*b-cCentre[spanDir]); // when z = [0, 2h]


        //scalar yplus = y*Retau/h;
        //scalar zplus = z*Retau/h;
        scalar yplus = y*Retau/a;
        scalar zplus = z*Retau/b;

        //scalar xplus = cCentre[streamDir]*Retau/h;

        if (setBulk)
        {
            // laminar parabolic profile
            U[celli] = vector::zero;

           // U[celli][streamDir] =
               // 3.0*Ubar.value()[streamDir] * (y/h - 0.5*sqr(y/h));
            for (double i = 1; i < 12; i += 2)
                {
                    U[celli][streamDir]  += 
                        48.0*Ubar.value()[streamDir] / (constant::mathematical::pi * constant::mathematical::pi * constant::mathematical::pi) * Foam::pow(-1.0, (i-1.0)/2.0)  
                        * (1.0 - Foam::cosh(i * constant::mathematical::pi * (z-b) / (2.0*b)) / Foam::cosh(i * constant::mathematical::pi * a / (2.0*b))) 
                        * Foam::cos(i * constant::mathematical::pi * (y-a) / (2.0*a)) / Foam::pow(i, 3.0);
                }
        }

        if (perturb)
        {
            // streak streamwise velocity
            U[celli][streamDir] +=
                (utau * duplus/2.0) * (yplus/40.0)
                * Foam::exp(-sigma * (Foam::sqr(yplus) + Foam::sqr(zplus))  + 0.5)
                * deviation;
        }
    }

    Info<< "Writing modified U field to " << runTime.timeName() << endl;
    U.write();

    Info<< endl;

    return(0);
}


// ************************************************************************* //
