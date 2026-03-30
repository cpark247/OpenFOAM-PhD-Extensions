/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::porousBafflePressureAMIFvPatchField

Description
    This boundary condition provides a porous baffle pressure jump on
    cyclicACMI (and cyclicAMI) patches, using the AMI interpolation
    framework to couple non-conformal mesh interfaces.

    This is an AMI-enabled re-implementation of porousBafflePressure.
    It inherits from fixedJumpAMIFvPatchField<scalar> instead of
    fixedJumpFvPatchField<scalar>, making it compatible with
    cyclicACMI patch types.

    The porous baffle pressure jump is defined by:

    \f[
        \Delta p = -\text{sign}(\phi) \left( D \mu |U_n| + \frac{1}{2} I \rho |U_n|^2 \right) L
    \f]

    \vartable
        \Delta p | pressure jump [Pa] or [m^2/s^2]
        \phi     | face flux, used to determine flow direction
        D        | Darcy coefficient [1/m^2]
        \mu      | dynamic viscosity [Pa.s] (compressible) or
                   kinematic viscosity [m^2/s] (incompressible)
        I        | inertial (Forchheimer) coefficient [1/m]
        \rho     | density [kg/m^3] (compressible case only)
        U_n      | face-normal velocity magnitude [m/s]
        L        | porous media length in the flow direction [m]
    \endvartable

    The flow direction is determined locally per-face from sign(phi),
    ensuring the pressure jump always opposes the flow.

    On cyclicACMI patches created from wall faces (e.g. via createBaffles
    on blade surfaces in rotating machinery), phi may have a persistently
    wrong sign because the face-normal velocity is a tiny residual of
    projecting a large tangential velocity onto a nearly perpendicular
    face normal.  The optional \c positive keyword enforces the expected
    jump sign: \c false forces jump <= 0 (pressure drop for owner-to-
    neighbour flow), \c true forces jump >= 0 (pressure drop for
    neighbour-to-owner flow).  Wrong-sign faces get their sign flipped
    while preserving the Darcy-Forchheimer magnitude.

    Additionally, phi on ACMI boundary faces inside an MRF zone can be
    contaminated by the rotational velocity correction, inflating Un
    (and hence the jump) far beyond the actual through-blade velocity.
    The optional \c jumpMax keyword bounds the magnitude to prevent the
    jump from exceeding the available pressure difference and reversing
    the pressure gradient.  Default is 1e3; for centrifugal pump MRF
    simulations, a value comparable to the expected pressure difference
    across the blade is recommended (e.g. 20-50 m^2/s^2).

    For multiphase solvers (e.g., multiphaseEulerFoam), use the optional
    \c phase keyword to specify which phase's turbulence model to use
    for the viscosity lookup. The \c phi and \c rho keywords can also
    be set to phase-specific fields (e.g., phi.water, rho.water).

Usage
    \table
        Property     | Description                          | Required | Default
        patchType    | underlying patch type (cyclicACMI)   | yes      |
        D            | Darcy coefficient [1/m^2]            | yes      |
        I            | inertial coefficient [1/m]           | yes      |
        length       | porous media length [m]              | yes      |
        phi          | flux field name                      | no       | phi
        rho          | density field name                   | no       | rho
        phase        | phase name for turbulence lookup     | no       | (none)
        positive     | enforce jump sign (true/false)       | no       | (none)
        jumpMax      | max allowed |jump| magnitude        | no       | 1e3
    \endtable

    Example for single-phase solver (e.g., simpleFoam):
    \verbatim
    hole_inlet
    {
        type            porousBafflePressureAMI;
        patchType       cyclicACMI;
        jump            uniform 0;
        D               0;
        I               100;
        length          0.015;
        value           uniform 0;
    }
    \endverbatim

    Example for multiphaseEulerFoam:
    \verbatim
    hole_inlet
    {
        type            porousBafflePressureAMI;
        patchType       cyclicACMI;
        jump            uniform 0;
        D               0;
        I               100;
        length          0.015;
        phi             phi.water;
        rho             rho.water;
        phase           water;
        value           uniform 0;
    }
    \endverbatim

    Example with enforced negative jump and magnitude bound:
    \verbatim
    hole_inlet
    {
        type            porousBafflePressureAMI;
        patchType       cyclicACMI;
        jump            uniform 0;
        D               0;
        I               250;
        length          0.005;
        positive        false;
        jumpMax         30;
        value           uniform 0;
    }
    \endverbatim

Note
    - The underlying patchType should be set to \c cyclicACMI or \c cyclicAMI
    - Both patches of the pair must use this boundary condition type
    - The jump is computed on the owner patch only; the neighbour
      automatically receives the negative jump via AMI interpolation
    - D = 0 is valid and recommended at high Reynolds numbers where
      the inertial term dominates
    - Based on the OpenFOAM Foundation v9 porousBafflePressure implementation

SourceFiles
    porousBafflePressureAMIFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef porousBafflePressureAMIFvPatchField_H
#define porousBafflePressureAMIFvPatchField_H

#include "fixedJumpAMIFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class porousBafflePressureAMIFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class porousBafflePressureAMIFvPatchField
:
    public fixedJumpAMIFvPatchField<scalar>
{
    // Private Data

        //- Name of flux field (default = phi)
        const word phiName_;

        //- Name of density field (default = rho)
        const word rhoName_;

        //- Name of phase for turbulence model lookup (default = "")
        const word phaseName_;

        //- Darcy pressure loss coefficient [1/m^2]
        scalar D_;

        //- Inertia pressure loss coefficient [1/m]
        scalar I_;

        //- Porous media length in the flow direction [m]
        scalar length_;

        //- Whether the 'positive' keyword was specified
        bool positiveSet_;

        //- If true, jump must be >= 0; if false, jump must be <= 0
        //  Only active when positiveSet_ is true.
        bool positive_;

        //- Maximum allowed |jump| magnitude [m^2/s^2 or Pa]
        scalar jumpMax_;


public:

    //- Runtime type information
    TypeName("porousBafflePressureAMI");


    // Constructors

        //- Construct from patch and internal field
        porousBafflePressureAMIFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        porousBafflePressureAMIFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given porousBafflePressureAMIFvPatchField
        //  onto a new patch
        porousBafflePressureAMIFvPatchField
        (
            const porousBafflePressureAMIFvPatchField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        porousBafflePressureAMIFvPatchField
        (
            const porousBafflePressureAMIFvPatchField&
        ) = delete;

        //- Copy constructor setting internal field reference
        porousBafflePressureAMIFvPatchField
        (
            const porousBafflePressureAMIFvPatchField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new porousBafflePressureAMIFvPatchField(*this, iF)
            );
        }


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
