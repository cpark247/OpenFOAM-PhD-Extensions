/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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

#include "LuoSvendsenForSAS2.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseDynamicMomentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace binaryBreakupModels
{
    defineTypeNameAndDebug(LuoSvendsenForSAS2, 0);
    addToRunTimeSelectionTable
    (
        binaryBreakupModel,
        LuoSvendsenForSAS2,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::binaryBreakupModels::LuoSvendsenForSAS2::LuoSvendsenForSAS2
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    binaryBreakupModel(popBal, dict),
    gammaUpperReg2by11_(),
    gammaUpperReg5by11_(),
    gammaUpperReg8by11_(),
    C4_(dimensionedScalar::lookupOrDefault("C4", dict, dimless, 0.923)),
    beta_(dimensionedScalar::lookupOrDefault("beta", dict, dimless, 2.05)),
    minEddyRatio_
    (
        dimensionedScalar::lookupOrDefault("minEddyRatio", dict, dimless, 11.4)
    ),
    kolmogorovLengthScale_
    (
        IOobject
        (
            "kolmogorovLengthScale",
            popBal_.time().timeName(),
            popBal_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        popBal_.mesh(),
        dimensionedScalar("kolmogorovLengthScale", dimLength, Zero)
    ),
    binaryBreakupRate_
    (
        IOobject
        (
            "binaryBreakupRateTotal",
            popBal_.time().timeName(),
            popBal_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE  // Write instantaneous total
        ),
        popBal_.mesh(),
        dimensionedScalar("binaryBreakupRateTotal", inv(dimVolume*dimTime), Zero)
    ),
    binaryBreakupRates_(),
    ratesInitialized_(false),
    lastResetTimeIndex_(-1)
{
    Info << "LuoSvendsenForSAS2(): Initializing breakup model" << endl;
    Info << "    C4 = " << C4_.value() << endl;
    Info << "    beta = " << beta_.value() << endl;
    Info << "    minEddyRatio = " << minEddyRatio_.value() << endl;

    // Build tabulated incomplete gamma function values
    List<Tuple2<scalar, scalar>> gammaUpperReg2by11Table;
    List<Tuple2<scalar, scalar>> gammaUpperReg5by11Table;
    List<Tuple2<scalar, scalar>> gammaUpperReg8by11Table;

    gammaUpperReg2by11Table.append(Tuple2<scalar, scalar>(0.0, 1.0));
    gammaUpperReg5by11Table.append(Tuple2<scalar, scalar>(0.0, 1.0));
    gammaUpperReg8by11Table.append(Tuple2<scalar, scalar>(0.0, 1.0));

    for (scalar z = 1e-2; z <= 10.0; z = z + 1e-2)
    {
        gammaUpperReg2by11Table.append
        (
            Tuple2<scalar, scalar>(z, incGammaRatio_Q(2.0/11.0, z))
        );
        gammaUpperReg5by11Table.append
        (
            Tuple2<scalar, scalar>(z, incGammaRatio_Q(5.0/11.0, z))
        );
        gammaUpperReg8by11Table.append
        (
            Tuple2<scalar, scalar>(z, incGammaRatio_Q(8.0/11.0, z))
        );
    }

    gammaUpperReg2by11_ =
        new Function1s::Table<scalar>
        (
            "gamma2by11",
            Function1s::tableBase::boundsHandling::clamp,
            linearInterpolationWeights::typeName,
            gammaUpperReg2by11Table
        );

    gammaUpperReg5by11_ =
        new Function1s::Table<scalar>
        (
            "gamma5by11",
            Function1s::tableBase::boundsHandling::clamp,
            linearInterpolationWeights::typeName,
            gammaUpperReg5by11Table
        );

    gammaUpperReg8by11_ =
        new Function1s::Table<scalar>
        (
            "gamma8by11",
            Function1s::tableBase::boundsHandling::clamp,
            linearInterpolationWeights::typeName,
            gammaUpperReg8by11Table
        );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::binaryBreakupModels::LuoSvendsenForSAS2::precompute()
{
    // Initialize per-class fields on first call
    if (!ratesInitialized_)
    {
        const label nSizeGroups = popBal_.sizeGroups().size();
        
        Info << "LuoSvendsenForSAS2::precompute(): Initializing " 
             << nSizeGroups << " per-class breakup rate fields" << endl;
        
        binaryBreakupRates_.setSize(nSizeGroups);
        
        forAll(popBal_.sizeGroups(), i)
        {
            const sizeGroup& fi = popBal_.sizeGroups()[i];
            
            // Field name: binaryBreakupRate_f0.air, binaryBreakupRate_f1.air, etc.
            word fieldName = "binaryBreakupRate_" + fi.name();
            
            Info << "    Creating field: " << fieldName 
                 << " (index " << i << ")" << endl;
            
            binaryBreakupRates_.set
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
                    dimensionedScalar(fieldName, inv(dimVolume*dimTime), Zero)
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
        binaryBreakupRate_ = dimensionedScalar
        (
            "zero", binaryBreakupRate_.dimensions(), Zero
        );
        
        // Reset per-class rates
        forAll(binaryBreakupRates_, i)
        {
            binaryBreakupRates_[i] = dimensionedScalar
            (
                "zero", binaryBreakupRates_[i].dimensions(), Zero
            );
        }
        
        lastResetTimeIndex_ = currentTimeIndex;
    }

    //******************** HSM Addon: Total dissipation rate ********************//
    const phaseModel& continuousPhase = popBal_.continuousPhase();
    volScalarField epsilonTot = popBal_.continuousTurbulence().epsilon();

    if (continuousPhase.db().foundObject<volScalarField>("epsilonRes."+continuousPhase.name()))
    {
        const volScalarField& epsilonRes = 
            continuousPhase.db().lookupObject<volScalarField>("epsilonRes."+continuousPhase.name());
        epsilonTot += epsilonRes;
    }
    //******************** End HSM Addon ********************//

    // Update Kolmogorov length scale
    kolmogorovLengthScale_ = pow025
    (
        pow3(popBal_.continuousPhase().thermo().nu())
       /max(epsilonTot, dimensionedScalar("minEps", epsilonTot.dimensions(), 1e-14))
    );
}


void
Foam::diameterModels::binaryBreakupModels::LuoSvendsenForSAS2::addToBinaryBreakupRate
(
    volScalarField& binaryBreakupRate,
    const label i,
    const label j
)
{
    // i = daughter class index, j = parent class index (j > i)
    const phaseModel& continuousPhase = popBal_.continuousPhase();
    const sizeGroup& fi = popBal_.sizeGroups()[i];  // Daughter
    const sizeGroup& fj = popBal_.sizeGroups()[j];  // Parent

    //******************** HSM Addon: Total dissipation rate ********************//
    volScalarField epsilonTot = popBal_.continuousTurbulence().epsilon();

    if (continuousPhase.db().foundObject<volScalarField>("epsilonRes."+continuousPhase.name()))
    {
        const volScalarField& epsilonRes = 
            continuousPhase.db().lookupObject<volScalarField>("epsilonRes."+continuousPhase.name());
        epsilonTot += epsilonRes;
    }
    //******************** End HSM Addon ********************//

    const volScalarField epsilonTotSafe = max
    (
        epsilonTot,
        dimensionedScalar("minEps", epsilonTot.dimensions(), 1e-14)
    );

    // Surface energy coefficient
    const dimensionedScalar cf
    (
        pow(fi.x()/fj.x(), 2.0/3.0) + pow((1 - fi.x()/fj.x()), 2.0/3.0) - 1
    );

    // Dimensionless surface energy parameter
    const volScalarField b
    (
        12*cf*popBal_.sigmaWithContinuousPhase(fi.phase())
       /(
            beta_*continuousPhase.rho()*pow(fj.dSph(), 5.0/3.0)
           *pow(epsilonTotSafe, 2.0/3.0)
        )
    );

    // Minimum eddy size ratio
    const volScalarField xiMin(minEddyRatio_*kolmogorovLengthScale_/fj.d());
    const volScalarField tMin(b/pow(xiMin, 11.0/3.0));

    // Compute integral using incomplete gamma functions
    volScalarField integral(3/(11*pow(b, 8.0/11.0)));

    forAll(integral, celli)
    {
        integral[celli] *=
            (
                tgamma(8.0/11.0)
               *(
                    gammaUpperReg8by11_->value(b[celli])
                  - gammaUpperReg8by11_->value(tMin[celli])
                )
              + 2*pow(b[celli], 3.0/11.0)*tgamma(5.0/11.0)
               *(
                    gammaUpperReg5by11_->value(b[celli])
                  - gammaUpperReg5by11_->value(tMin[celli])
                )
              + pow(b[celli], 6.0/11.0)*tgamma(2.0/11.0)
               *(
                    gammaUpperReg2by11_->value(b[celli])
                  - gammaUpperReg2by11_->value(tMin[celli])
                )
            );
    }

    // Compute breakup rate contribution for this (i,j) pair
    volScalarField rateContribution
    (
        C4_*(1 - popBal_.alphas())/fj.x()
       *cbrt(epsilonTotSafe/sqr(fj.dSph()))
       *integral
    );

    // Add to population balance rate (passed by reference)
    binaryBreakupRate += rateContribution;

    // Accumulate to total rate
    binaryBreakupRate_ += rateContribution;

    // Accumulate to per-daughter-class rate
    if (ratesInitialized_ && i < binaryBreakupRates_.size())
    {
        binaryBreakupRates_[i] += rateContribution;
    }
}


// ************************************************************************* //
