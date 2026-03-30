/*---------------------------------------------------------------------------*\
Date: June 24, 2005
Author: Eugene de Villiers
Modified: December 2025
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
    Initialise duct velocity with optional streamwise streaks.
    To be used to force the flow to transition and reach a fully
    developed state sooner.

    Reads in perturbUDict.

    References:
        Schoppa, W. and Hussain, F.
        "Coherent structure dynamics in near wall turbulence",
        Fluid Dynamics Research, Vol 26, pp119-139, 2000.

Modifications:
    - Added yMid and zMid parameters for arbitrary duct center position
    - Improved coordinate transformations and wall distance calculations
    - Removed hardcoded offsets in Fourier series
    - Added validation checks for geometry

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

    // Read perturbation dictionary
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

    // Read parameters
    const scalar Retau(readScalar(perturbDict.lookup("Retau")));
    const scalar h(readScalar(perturbDict.lookup("h")));
    const bool setBulk(readBool(perturbDict.lookup("setBulk")));
    const bool perturb(readBool(perturbDict.lookup("perturb")));
    const direction streamDir(readLabel(perturbDict.lookup("streamwise")));
    const direction spanDir(readLabel(perturbDict.lookup("spanwise")));
    
    // Read duct center positions (with default values for backward compatibility)
    const scalar yMid
    (
        perturbDict.lookupOrDefault<scalar>("yMid", h)
    );
    const scalar zMid
    (
        perturbDict.lookupOrDefault<scalar>("zMid", h)
    );

    Info<< "Duct half height          = " << h << nl
        << "Duct center (y)           = " << yMid << nl
        << "Duct center (z)           = " << zMid << nl
        << "Re(tau)                   = " << Retau << nl
        << "Set bulk flow             = " << Switch(setBulk) << nl
        << "Perturb flow              = " << Switch(perturb) << nl
        << "Streamwise flow component = " << streamDir << nl
        << "Spanwise flow component   = " << spanDir << nl
        << endl;

    // Validation checks
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

    // Get component normal to streamDir and spanDir (height component)
    direction heightDir = 0;
    if (streamDir == heightDir)
    {
        heightDir++;
    }
    if (spanDir == heightDir)
    {
        heightDir++;
    }

    // Read velocity field
    IOobject Uheader
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );
    Info<< "Reading U" << endl;
    volVectorField U(Uheader, mesh);

    // Read transport properties
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
    const scalar utau = Retau*nu.value()/h;
    Info<< "u(tau)  = " << utau << endl;

    // Perturbation parameters
    const scalar duplus = Ubar.value()[streamDir]*0.25/utau;
    const scalar sigma = 0.00055;

    // Random number generator
    Random perturbation(1234567);

    const vectorField& centres = mesh.C();

    // Statistics for validation
    scalar minHeight = GREAT;
    scalar maxHeight = -GREAT;
    scalar minSpan = GREAT;
    scalar maxSpan = -GREAT;

    forAll(centres, celli)
    {
        // Add small (+/-20%) random component to enhance symmetry breaking
        scalar deviation = 1.0 + 0.2*perturbation.scalarNormal();

        const vector& cCentre = centres[celli];

        // Transform coordinates to duct-centered system
        // y_local and z_local are coordinates relative to duct center
        // They range from -h to +h in each direction
        scalar y_local = cCentre[heightDir] - yMid;
        scalar z_local = cCentre[spanDir] - zMid;

        // Distance from nearest wall in each direction
        // This accounts for the duct being centered at (yMid, zMid)
        scalar y_wall = h - Foam::mag(y_local);  // distance to top/bottom wall
        scalar z_wall = h - Foam::mag(z_local);  // distance to side walls

        // Wall units (based on distance from nearest wall)
        //scalar yplus = Foam::mag(y_local)*Retau/h;  // Distance from center line
        //scalar zplus = Foam::mag(z_local)*Retau/h;  // Distance from center line
        
        // For perturbation, we want distance from nearest wall
        scalar yplus_wall = y_wall*Retau/h;
        scalar zplus_wall = z_wall*Retau/h;

        // Track domain extent for validation
        minHeight = min(minHeight, y_wall);
        maxHeight = max(maxHeight, Foam::mag(y_local));
        minSpan = min(minSpan, z_wall);
        maxSpan = max(maxSpan, Foam::mag(z_local));

        if (setBulk)
        {
            // Initialize with fully developed laminar duct flow (Fourier series solution)
            // This is the analytical solution for Poiseuille flow in a rectangular duct
            U[celli] = vector::zero;

            // Fourier series: sum odd terms for symmetry
            // Note: coordinates are now properly centered at (0,0)
            for (double i = 1; i < 12; i += 2)
            {
                U[celli][streamDir] += 
                    48.0*Ubar.value()[streamDir] 
                    / (constant::mathematical::pi * constant::mathematical::pi * constant::mathematical::pi) 
                    * Foam::pow(-1.0, (i-1.0)/2.0)  
                    * (1.0 - Foam::cosh(i * constant::mathematical::pi * z_local / (2.0*h)) 
                            / Foam::cosh(i * constant::mathematical::pi * h / (2.0*h)))
                    * Foam::cos(i * constant::mathematical::pi * y_local / (2.0*h)) 
                    / Foam::pow(i, 3.0);
            }
        }

        if (perturb)
        {
            // Add streak perturbation based on wall distance
            // This perturbation is strongest near the wall and decays exponentially
            U[celli][streamDir] +=
                (utau * duplus/2.0) * (yplus_wall/40.0)
                * Foam::exp(-sigma * (Foam::sqr(yplus_wall) + Foam::sqr(zplus_wall)) + 0.5)
                * deviation;
        }
    }

    // Validation output
    Info<< nl << "Domain validation:" << nl
        << "  Max |y-yMid|        = " << maxHeight << " (should be <= " << h << ")" << nl
        << "  Min wall distance y = " << minHeight << " (should be >= 0)" << nl
        << "  Max |z-zMid|        = " << maxSpan << " (should be <= " << h << ")" << nl
        << "  Min wall distance z = " << minSpan << " (should be >= 0)" << nl;

    if (minHeight < -SMALL || minSpan < -SMALL)
    {
        WarningIn(args.executable())
            << "Some cells appear to be outside the duct boundaries!" << nl
            << "Check that yMid=" << yMid << " and zMid=" << zMid 
            << " are correct for your geometry." << endl;
    }

    if (maxHeight > h + SMALL || maxSpan > h + SMALL)
    {
        WarningIn(args.executable())
            << "Some cells extend beyond the duct half-height h=" << h << nl
            << "This might indicate a geometry mismatch." << endl;
    }

    Info<< nl << "Writing modified U field to " << runTime.timeName() << endl;
    U.write();

    Info<< nl << "End" << nl << endl;

    return(0);
}


// ************************************************************************* //
