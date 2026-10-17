/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

#include "PrinceBlanchForSAS2.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "phaseDynamicMomentumTransportModel.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace coalescenceModels
{
    defineTypeNameAndDebug(PrinceBlanchForSAS2, 0);
    addToRunTimeSelectionTable
    (
        coalescenceModel,
        PrinceBlanchForSAS2,
        dictionary
    );
}
}
}

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::coalescenceModels::PrinceBlanchForSAS2::
PrinceBlanchForSAS2
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    coalescenceModel(popBal, dict),
    PMax_("PMax", dimless, dict.lookupOrDefault<scalar>("PMax", 0.8)),
    CPack_
    (
        IOobject
        (
            "CPack",
            popBal_.time().timeName(),
            popBal_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
        popBal_.mesh(),
        dimensionedScalar
        (
            "CPack",
            dimless,
            Zero
        )
    ),
    CPackMax_
    (
        "CPackMax", dimless, dict.lookupOrDefault<scalar>("CPackMax", 1e5)
    ),
    C1_("C1", dimless, dict.lookupOrDefault<scalar>("C1", 0.356)),
    h0_("h0", dimLength, dict.lookupOrDefault<scalar>("h0", 1e-4)),
    hf_("hf", dimLength, dict.lookupOrDefault<scalar>("hf", 1e-8)),
    turbulentCollisions_(dict.lookup("turbulentCollisions")),
    buoyantCollisions_(dict.lookup("buoyantCollisions")),
    laminarShearCollisions_(dict.lookup("laminarShearCollisions")),
    coalescenceRate_
    (
         IOobject
         (
            "coalescenceRateTotal",
             popBal_.time().timeName(),
             popBal_.mesh(),
             IOobject::NO_READ,
             IOobject::AUTO_WRITE  // Write instantaneous total
         ),
         popBal_.mesh(),
         dimensionedScalar
         (
            "coalescenceRateTotal",
            dimVolume/dimTime,
            Zero
         )
    ),
    coalescenceRates_(),
    ratesInitialized_(false),
    lastI_(-1),
    lastJ_(-1)
{
    Info << "PrinceBlanchForSAS2(): Initializing coalescence model" << endl;
    Info << "    PMax = " << PMax_.value() << endl;
    Info << "    CPackMax = " << CPackMax_.value() << endl;
    Info << "    C1 = " << C1_.value() << endl;
    Info << "    h0 = " << h0_.value() << endl;
    Info << "    hf = " << hf_.value() << endl;
    Info << "    turbulentCollisions = " << turbulentCollisions_ << endl;
    Info << "    buoyantCollisions = " << buoyantCollisions_ << endl;
    Info << "    laminarShearCollisions = " << laminarShearCollisions_ << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::coalescenceModels::PrinceBlanchForSAS2::precompute()
{
    // Initialize per-class fields on first call
    if (!ratesInitialized_)
    {
        const label nSizeGroups = popBal_.sizeGroups().size();
        
        Info << "PrinceBlanchForSAS2::precompute(): Initializing " 
             << nSizeGroups << " per-class coalescence rate fields" << endl;
        
        coalescenceRates_.setSize(nSizeGroups);
        
        forAll(popBal_.sizeGroups(), i)
        {
            const sizeGroup& fi = popBal_.sizeGroups()[i];
            
            // Field name: coalescenceRate_f0.air, coalescenceRate_f1.air, etc.
            word fieldName = "coalescenceRate_" + fi.name();
            
            Info << "    Creating field: " << fieldName 
                 << " (index " << i << ")" << endl;
            
            coalescenceRates_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        fieldName,
                        popBal_.time().timeName(),
                        popBal_.mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE  // Instantaneous per-class NOT written
                    ),
                    popBal_.mesh(),
                    dimensionedScalar(fieldName, dimVolume/dimTime, Zero)
                )
            );
        }
        
        ratesInitialized_ = true;
    }
}


void Foam::diameterModels::coalescenceModels::PrinceBlanchForSAS2::
addToCoalescenceRate
(
    volScalarField& coalescenceRate,
    const label i,
    const label j
)
{
    const phaseModel& continuousPhase = popBal_.continuousPhase();
    const sizeGroup& fi = popBal_.sizeGroups()[i];
    const sizeGroup& fj = popBal_.sizeGroups()[j];
    const uniformDimensionedVectorField& g =
        popBal_.mesh().lookupObject<uniformDimensionedVectorField>("g");

    const dimensionedScalar rij(1/(1/fi.dSph() + 1/fj.dSph()));

    // Detect the start of a new source-term sweep. populationBalanceModel
    // walks the pairs in lexicographic (i, j) order, so a pair that does not
    // advance that order marks a fresh sweep. Resetting on the time index
    // alone would not do: the pair loop runs once per population balance
    // correction AND once per PIMPLE outer corrector.
    if (ratesInitialized_ && (i < lastI_ || (i == lastI_ && j <= lastJ_)))
    {
        coalescenceRate_ = dimensionedScalar
        (
            "zero", coalescenceRate_.dimensions(), Zero
        );

        forAll(coalescenceRates_, k)
        {
            coalescenceRates_[k] = dimensionedScalar
            (
                "zero", coalescenceRates_[k].dimensions(), Zero
            );
        }
    }

    lastI_ = i;
    lastJ_ = j;

//********************************* HSM Addon **************************//
    CPack_ = pos0(popBal_.alphas() - PMax_) * min(sqr(1/max(1 - popBal_.alphas(), SMALL)), CPackMax_) + neg(popBal_.alphas() - PMax_);
   
    // Total dissipation rate seen by the kernel.
    // Preferred source: epsilonTotal.<phase>, registered by the resolved-LES
    // solver as
    //     epsilon = epsilonRes + epsilonMod = 2*(nu + nu_t)*s_ij*s_ij
    // with s_ij the INSTANTANEOUS resolved strain rate (Pope 2000). The
    // turbulence model's own epsilon() is an equilibrium estimate
    // (Ce*k_sgs^1.5/delta) and must not be mixed in when epsilonTotal is
    // available, otherwise the sub-grid part would be counted twice.
    //
    // Fallbacks, in order:
    //   1. turbulence-model epsilon + epsilonRes.<phase>  (legacy behaviour)
    //   2. turbulence-model epsilon alone                 (RANS)
    const objectRegistry& db = continuousPhase.db();

    volScalarField epsilonTot_(popBal_.continuousTurbulence().epsilon());

    if (db.foundObject<volScalarField>("epsilonTotal."+continuousPhase.name()))
    {
        epsilonTot_ =
            db.lookupObject<volScalarField>
            (
                "epsilonTotal."+continuousPhase.name()
            );
    }
    else if (db.foundObject<volScalarField>("epsilonRes."+continuousPhase.name()))
    {
        epsilonTot_ +=
            db.lookupObject<volScalarField>
            (
                "epsilonRes."+continuousPhase.name()
            );
    }
    else
    {
        Warning
            << "Neither epsilonTotal." << continuousPhase.name()
            << " nor epsilonRes." << continuousPhase.name()
            << " is registered." << endl
            << "Falling back to the turbulence model dissipation rate only."
            << endl;
    }

//********************************* End HSM Addon **************************//


 
    const volScalarField collisionEfficiency
    (
        exp
        (
          - sqrt
            (
                pow3(rij)*continuousPhase.rho()
               /(16.0*popBal_.sigmaWithContinuousPhase(fi.phase()))
            )
           *log(h0_/hf_)
           *cbrt(epsilonTot_)/pow(rij, 2.0/3.0)
        )
    );

    // Temporary field to store rate contribution for this (i,j) pair
    volScalarField rateContribution
    (
        IOobject
        (
            "rateContribution",
            popBal_.time().timeName(),
            popBal_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        popBal_.mesh(),
        dimensionedScalar("zero", dimVolume/dimTime, Zero)
    );

    if (turbulentCollisions_)
    {
        rateContribution +=
            (
                C1_*pi*sqr(fi.dSph() + fj.dSph())
               *cbrt(epsilonTot_)
               *sqrt(pow(fi.dSph(), 2.0/3.0) + pow(fj.dSph(), 2.0/3.0))
	       *CPack_
            )
           *collisionEfficiency;
    }

    if (buoyantCollisions_)
    {
        const dimensionedScalar Sij(pi/4*sqr(fi.dSph() + fj.dSph()));

        rateContribution +=
            (
                Sij
               *mag
                (
                    sqrt
                    (
                        2.14*popBal_.sigmaWithContinuousPhase(fi.phase())
                       /(continuousPhase.rho()*fi.dSph())
		      + 0.505*mag(g)*fi.dSph()
                    )
                  - sqrt
                    (
                        2.14*popBal_.sigmaWithContinuousPhase(fi.phase())
                       /(continuousPhase.rho()*fj.dSph()) 
		       + 0.505*mag(g)*fj.dSph()
                    )
                )
               *CPack_
	    )
           *collisionEfficiency;
    }

    if (laminarShearCollisions_)
    {
        FatalErrorInFunction
            << "Laminar shear collision contribution not implemented for "
            << this->type() << " coalescence model."
            << exit(FatalError);
    }

    // Add to the population balance coalescence rate (passed by reference)
    coalescenceRate += rateContribution;

    // ---- Diagnostics: size-fraction weighted mean coalescence kernel ----
    // rateContribution is the kernel a(v_i, v_j) [m^3/s]. Weighting it with
    // the size fractions gives the expected kernel for two independently drawn
    // bubbles, the class-method analogue of evaluating the kernel at a single
    // representative diameter. The solver converts it to a frequency with the
    // Buffo et al. (2016) relation, a*alpha_d/((pi/6)*d_32^3), whose
    // reciprocal is the coalescence time scale.
    //
    // addToCoalescenceRate() is called once per UNORDERED pair (j <= i), so
    // the off-diagonal pairs carry a factor 2 to complete the double sum
    //     sum_i sum_j f_i*f_j*a_ij.
    if (ratesInitialized_)
    {
        // Mean kernel seen by each participating class, sum_j f_j*a_ij
        coalescenceRates_[i] += rateContribution*fj;

        if (i != j)
        {
            coalescenceRates_[j] += rateContribution*fi;
        }

        if (i == j)
        {
            coalescenceRate_ += rateContribution*fi*fj;
        }
        else
        {
            coalescenceRate_ += 2.0*rateContribution*fi*fj;
        }
    }
}


// ************************************************************************* //
