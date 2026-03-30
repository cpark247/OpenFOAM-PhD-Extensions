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
    hf_("hf", dimLength, dict.lookupOrDefault<scalar>("h0", 1e-8)),
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
    lastResetTimeIndex_(-1)
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

    // Reset rate fields at the start of each new timestep
    const label currentTimeIndex = popBal_.mesh().time().timeIndex();
    
    if (currentTimeIndex != lastResetTimeIndex_)
    {
        // Reset total rate
        coalescenceRate_ = dimensionedScalar
        (
            "zero", coalescenceRate_.dimensions(), Zero
        );
        
        // Reset per-class rates
        forAll(coalescenceRates_, i)
        {
            coalescenceRates_[i] = dimensionedScalar
            (
                "zero", coalescenceRates_[i].dimensions(), Zero
            );
        }
        
        lastResetTimeIndex_ = currentTimeIndex;
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

//********************************* HSM Addon **************************//
    CPack_ = pos0(popBal_.alphas() - PMax_) * min(sqr(1/max(1 - popBal_.alphas(), SMALL)), CPackMax_) + neg(popBal_.alphas() - PMax_);
   
    volScalarField epsilonRes_ = scalar(0)*popBal_.continuousTurbulence().epsilon();
    volScalarField epsilonTot_ = popBal_.continuousTurbulence().epsilon();

    if (continuousPhase.db().foundObject<volScalarField>("epsilonRes."+continuousPhase.name()))
    {
	const objectRegistry& db = continuousPhase.db();
        epsilonRes_ = db.lookupObject<volScalarField>("epsilonRes."+continuousPhase.name());
        epsilonTot_ += epsilonRes_;
    }
    else
    {
	 Warning << "epsilonRes."+continuousPhase.name()+ " is not active in controlDict." << endl
	         << "Add field epsilonRes to the fieldAverage utility in controlDict." << endl;
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

    // Accumulate to total rate
    coalescenceRate_ += rateContribution;

    // Accumulate to per-class rates for BOTH participating classes
    // When bubbles from class i and j coalesce, both classes contribute
    if (ratesInitialized_)
    {
        if (i < coalescenceRates_.size())
        {
            coalescenceRates_[i] += rateContribution;
        }
        if (j < coalescenceRates_.size())
        {
            coalescenceRates_[j] += rateContribution;
        }
    }
}


// ************************************************************************* //
