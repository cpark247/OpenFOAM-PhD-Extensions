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

\*---------------------------------------------------------------------------*/

#include "porousBafflePressureAMIFvPatchField.H"
#include "surfaceFields.H"
#include "fvsPatchFields.H"
#include "momentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::porousBafflePressureAMIFvPatchField::
porousBafflePressureAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpAMIFvPatchField<scalar>(p, iF),
    phiName_("phi"),
    rhoName_("rho"),
    phaseName_(word::null),
    D_(0),
    I_(0),
    length_(0),
    positiveSet_(false),
    positive_(false),
    jumpMax_(1e3)
{}


Foam::porousBafflePressureAMIFvPatchField::
porousBafflePressureAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedJumpAMIFvPatchField<scalar>(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    phaseName_(dict.lookupOrDefault<word>("phase", word::null)),
    D_(dict.lookup<scalar>("D")),
    I_(dict.lookup<scalar>("I")),
    length_(dict.lookup<scalar>("length")),
    positiveSet_(dict.found("positive")),
    positive_(dict.lookupOrDefault<bool>("positive", false)),
    jumpMax_(dict.lookupOrDefault<scalar>("jumpMax", 1e3))
{
    // Sanity checks
    if (length_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Negative porous media length: " << length_
            << exit(FatalIOError);
    }

    if (D_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Negative Darcy coefficient D: " << D_
            << exit(FatalIOError);
    }

    if (I_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Negative inertial coefficient I: " << I_
            << exit(FatalIOError);
    }

    fvPatchField<scalar>::operator=
    (
        Field<scalar>("value", dict, p.size())
    );
}


Foam::porousBafflePressureAMIFvPatchField::
porousBafflePressureAMIFvPatchField
(
    const porousBafflePressureAMIFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedJumpAMIFvPatchField<scalar>(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    phaseName_(ptf.phaseName_),
    D_(ptf.D_),
    I_(ptf.I_),
    length_(ptf.length_),
    positiveSet_(ptf.positiveSet_),
    positive_(ptf.positive_),
    jumpMax_(ptf.jumpMax_)
{}


Foam::porousBafflePressureAMIFvPatchField::
porousBafflePressureAMIFvPatchField
(
    const porousBafflePressureAMIFvPatchField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpAMIFvPatchField<scalar>(ptf, iF),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    phaseName_(ptf.phaseName_),
    D_(ptf.D_),
    I_(ptf.I_),
    length_(ptf.length_),
    positiveSet_(ptf.positiveSet_),
    positive_(ptf.positive_),
    jumpMax_(ptf.jumpMax_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porousBafflePressureAMIFvPatchField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Only compute the jump on the owner patch.
    // The neighbour patch receives the (negated, AMI-interpolated) jump
    // automatically via the fixedJumpAMI base class.
    if (this->cyclicAMIPatch().owner())
    {
        // -----------------------------------------------------------------
        // 1. Face flux and face-normal velocity
        // -----------------------------------------------------------------
        const surfaceScalarField& phi =
            db().lookupObject<surfaceScalarField>(phiName_);

        const fvsPatchScalarField& phip =
            phi.boundaryField()[patch().index()];

        // Face-normal velocity magnitude.
        //
        // For incompressible solvers: phi is volumetric [m^3/s]
        //     Un = |phi| / Sf
        //
        // For multiphaseEulerFoam: phi.<phase> is the phase volumetric
        //     flux = alpha * U . Sf [m^3/s], so at alpha ~ 1 (hole filled
        //     with water), Un approximates the true phase velocity.
        //
        // NOTE: For compressible single-phase solvers where phi is mass
        //     flux [kg/s], this would need division by rho as well. This
        //     implementation assumes volumetric phi, which is correct for
        //     both incompressible solvers and multiphaseEulerFoam.

        const scalarField Un(mag(phip)/patch().magSf());


        // -----------------------------------------------------------------
        // 2. Turbulence model lookup (phase-aware)
        // -----------------------------------------------------------------
        // Single-phase:    looks up "momentumTransport"
        // Multiphase:      looks up "momentumTransport.<phaseName>"
        //   e.g. phase = water  ->  "momentumTransport.water"

        const word turbModelName =
        (
            phaseName_ != word::null
          ? IOobject::groupName
            (
                momentumTransportModel::typeName,
                phaseName_
            )
          : word(momentumTransportModel::typeName)
        );

        const momentumTransportModel& turbModel =
            db().lookupObject<momentumTransportModel>(turbModelName);


        // -----------------------------------------------------------------
        // 3. Compute Darcy-Forchheimer pressure jump
        // -----------------------------------------------------------------
        //
        // General form:
        //   dp = -sign(phi) * (D * mu * |Un| + 0.5 * I * rho * |Un|^2) * L
        //
        // sign(phi) ensures the jump always opposes the flow:
        //   - Flow owner -> neighbour (phi > 0): pressure drops
        //   - Flow neighbour -> owner (phi < 0): pressure rises
        //
        // Compressible (rho field found):
        //   mu = rho * nu,  jump in [Pa]
        //
        // Incompressible (no rho field):
        //   uses nu directly, jump in [m^2/s^2]

        if (db().foundObject<volScalarField>(rhoName_))
        {
            // --- Compressible / multiphase path ---
            const scalarField rhop
            (
                patch().lookupPatchField<volScalarField, scalar>(rhoName_)
            );

            const scalarField mup
            (
                rhop
              * turbModel.nu()().boundaryField()[patch().index()]
            );

            this->jump_ = -sign(phip)
                *(D_*mup*Un + 0.5*I_*rhop*sqr(Un))
                *length_;
        }
        else
        {
            // --- Incompressible path ---
            const scalarField nup
            (
                turbModel.nu()().boundaryField()[patch().index()]
            );

            this->jump_ = -sign(phip)
                *(D_*nup*Un + 0.5*I_*sqr(Un))
                *length_;
        }


        // -----------------------------------------------------------------
        // 4. Enforce jump sign (optional) and magnitude bound
        // -----------------------------------------------------------------
        //
        // On cyclicACMI patches created from wall faces (e.g. via
        // createBaffles on blade surfaces in rotating machinery), the
        // face flux phi can have a persistently wrong sign because the
        // face-normal velocity component is a tiny residual of
        // projecting a large tangential velocity onto a nearly
        // perpendicular face normal.  An incorrect phi sign produces
        // a wrong-sign jump that creates positive feedback and
        // destroys the pressure field.
        //
        // Additionally, phi on ACMI boundary faces inside an MRF zone
        // can be contaminated by the MRF rotational velocity correction
        //   phi_rel = phi_abs - (omega x r) . Sf
        // which inflates |phi| (and hence Un) far beyond the actual
        // through-blade velocity.  This produces jump magnitudes that
        // can exceed the available pressure difference, reversing the
        // local pressure gradient and destroying the field.
        //
        // The optional "positive" keyword enforces the expected sign:
        //   positive  false;   // jump must be <= 0
        //   positive  true;    // jump must be >= 0
        //
        // The optional "jumpMax" keyword bounds the magnitude:
        //   jumpMax   50;      // |jump| <= 50 on every face
        //
        // If not specified, defaults are: no sign enforcement,
        // jumpMax = 1e3.

        label nFlipped = 0;
        label nBounded = 0;
        const label nTotalFaces =
            returnReduce(patch().size(), sumOp<label>());

        if (positiveSet_)
        {
            if (positive_)
            {
                // Jump must be positive: flip any negative values
                forAll(this->jump_, facei)
                {
                    if (this->jump_[facei] < 0)
                    {
                        this->jump_[facei] = mag(this->jump_[facei]);
                        nFlipped++;
                    }
                }
            }
            else
            {
                // Jump must be negative: flip any positive values
                forAll(this->jump_, facei)
                {
                    if (this->jump_[facei] > 0)
                    {
                        this->jump_[facei] = -mag(this->jump_[facei]);
                        nFlipped++;
                    }
                }
            }

            reduce(nFlipped, sumOp<label>());
        }

        // Magnitude bounding
        forAll(this->jump_, facei)
        {
            if (mag(this->jump_[facei]) > jumpMax_)
            {
                this->jump_[facei] =
                    Foam::sign(this->jump_[facei]) * jumpMax_;
                nBounded++;
            }
        }
        reduce(nBounded, sumOp<label>());


        // -----------------------------------------------------------------
        // 5. Diagnostic output
        // -----------------------------------------------------------------
        Info<< "porousBafflePressureAMI " << patch().name()
            << " : jump min/avg/max = "
            << gMin(this->jump_) << " / "
            << gAverage(this->jump_) << " / "
            << gMax(this->jump_);

        if (nFlipped > 0)
        {
            Info<< "  [FLIPPED " << nFlipped << "/" << nTotalFaces << "]";
        }
        if (nBounded > 0)
        {
            Info<< "  [BOUNDED " << nBounded << "/" << nTotalFaces << "]";
        }

        Info<< endl;
    }

    // Base class handles AMI interpolation of jump to the neighbour patch
    fixedJumpAMIFvPatchField<scalar>::updateCoeffs();
}


void Foam::porousBafflePressureAMIFvPatchField::write(Ostream& os) const
{
    fixedJumpAMIFvPatchField<scalar>::write(os);

    // Explicitly write patchType from the underlying patch.
    // This is required because the dictionary constructor passes (p, iF)
    // to the base class (not the dict), so fvPatchField::patchType_ is
    // never set from the dictionary. Without this, decomposed fields
    // lack the "patchType cyclicACMI" keyword, causing the parallel
    // solver to fail the patch/patchField consistency check.
    writeEntry(os, "patchType", patch().type());

    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);

    if (phaseName_ != word::null)
    {
        writeEntry(os, "phase", phaseName_);
    }

    writeEntry(os, "D", D_);
    writeEntry(os, "I", I_);
    writeEntry(os, "length", length_);

    if (positiveSet_)
    {
        writeEntry(os, "positive", positive_);
    }

    if (jumpMax_ < 1e3 - SMALL)
    {
        writeEntry(os, "jumpMax", jumpMax_);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        porousBafflePressureAMIFvPatchField
    );
}

// ************************************************************************* //
