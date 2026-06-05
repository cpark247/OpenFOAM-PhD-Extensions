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

#include "LuoSvendsenForSST.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseDynamicMomentumTransportModel.H"
#include "linearInterpolationWeights.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace binaryBreakupModels
{
    defineTypeNameAndDebug(LuoSvendsenForSST, 0);
    addToRunTimeSelectionTable
    (
        binaryBreakupModel,
        LuoSvendsenForSST,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::binaryBreakupModels::LuoSvendsenForSST::LuoSvendsenForSST
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
            popBal_.mesh()
        ),
        popBal_.mesh(),
        dimensionedScalar("kolmogorovLengthScale", dimLength, Zero)
    ),
    binaryBreakupRate_
    (
        IOobject
        (
            "binaryBreakupRateTotal",   // FIX: was "binaryBreakupRate_"
            popBal_.time().timeName(),
            popBal_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        popBal_.mesh(),
        dimensionedScalar("binaryBreakupRateTotal", inv(dimVolume*dimTime), Zero)
    ),
    binaryBreakupRates_(),
    ratesInitialized_(false),
    lastResetTimeIndex_(-1)
{
    Info << "LuoSvendsenForSST(): Initializing breakup model" << endl;
    Info << "    C4 = " << C4_.value() << endl;
    Info << "    beta = " << beta_.value() << endl;
    Info << "    minEddyRatio = " << minEddyRatio_.value() << endl;

    List<Tuple2<scalar, scalar>> gammaUpperReg2by11Table;
    List<Tuple2<scalar, scalar>> gammaUpperReg5by11Table;
    List<Tuple2<scalar, scalar>> gammaUpperReg8by11Table;

    gammaUpperReg2by11Table.append(Tuple2<scalar, scalar>(0.0, 1.0));
    gammaUpperReg5by11Table.append(Tuple2<scalar, scalar>(0.0, 1.0));
    gammaUpperReg8by11Table.append(Tuple2<scalar, scalar>(0.0, 1.0));

    // Fine-resolution region (z = 0.01 to 10): step 0.01
    // Covers LES/high-turbulence cases where b < 10.
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

    // Coarser-resolution extension (z = 10.1 to 500): step 0.1
    // Required for RANS/SST flows where the dimensionless surface-energy
    // parameter b = 12*cf*sigma/(beta*rho*eps^(2/3)*d^(5/3)) >> 10.
    // Without this range, Q(a,b) and Q(a,t_min) are both clamped to
    // Q(a,10), their difference evaluates to zero, and the breakup
    // integral collapses to zero for all RANS cells.
    for (scalar z = 10.1; z <= 500.0; z = z + 0.1)
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

void Foam::diameterModels::binaryBreakupModels::LuoSvendsenForSST::precompute()
{
    // FIX: initialize per-class fields on first call
    if (!ratesInitialized_)
    {
        const label nSizeGroups = popBal_.sizeGroups().size();

        Info << "LuoSvendsenForSST::precompute(): Initializing "
             << nSizeGroups << " per-class breakup rate fields" << endl;

        binaryBreakupRates_.setSize(nSizeGroups);

        forAll(popBal_.sizeGroups(), i)
        {
            const sizeGroup& fi = popBal_.sizeGroups()[i];

            word fieldName = "binaryBreakupRate_" + fi.name();

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
                        IOobject::NO_WRITE
                    ),
                    popBal_.mesh(),
                    dimensionedScalar(fieldName, inv(dimVolume*dimTime), Zero)
                )
            );
        }

        ratesInitialized_ = true;
    }

    // FIX: reset rate fields once per timestep
    const label currentTimeIndex = popBal_.mesh().time().timeIndex();

    if (currentTimeIndex != lastResetTimeIndex_)
    {
        binaryBreakupRate_ = dimensionedScalar
        (
            "zero", binaryBreakupRate_.dimensions(), Zero
        );

        forAll(binaryBreakupRates_, i)
        {
            binaryBreakupRates_[i] = dimensionedScalar
            (
                "zero", binaryBreakupRates_[i].dimensions(), Zero
            );
        }

        lastResetTimeIndex_ = currentTimeIndex;
    }

    // Compute Kolmogorov length scale from RANS modelled epsilon
    volScalarField epsilonTot_ = popBal_.continuousTurbulence().epsilon();

    kolmogorovLengthScale_ =
        pow025
        (
            pow3(popBal_.continuousPhase().thermo().nu())
           /max
            (
                epsilonTot_,
                dimensionedScalar("minEps", epsilonTot_.dimensions(), 1e-14)
            )
        );
}


void
Foam::diameterModels::binaryBreakupModels::LuoSvendsenForSST::addToBinaryBreakupRate
(
    volScalarField& binaryBreakupRate,
    const label i,
    const label j
)
{
    const phaseModel& continuousPhase = popBal_.continuousPhase();
    const sizeGroup& fi = popBal_.sizeGroups()[i];   // daughter
    const sizeGroup& fj = popBal_.sizeGroups()[j];   // parent

    // RANS modelled epsilon only (no resolved component for SST)
    volScalarField epsilonTot_ = popBal_.continuousTurbulence().epsilon();

    const volScalarField epsilonTotSafe = max
    (
        epsilonTot_,
        dimensionedScalar("minEps", epsilonTot_.dimensions(), 1e-14)
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

    // FIX: compute integral with all three gamma-function terms
    // (previously only the 5/11 term was included)
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

    // Rate contribution for this (i,j) pair
    volScalarField rateContribution
    (
        C4_*(1 - popBal_.alphas())/fj.x()
       *cbrt(epsilonTotSafe/sqr(fj.dSph()))
       *integral
    );

    // Add to population balance rate (passed by reference)
    binaryBreakupRate += rateContribution;

    // FIX: accumulate (not overwrite) total and per-class diagnostic rates
    binaryBreakupRate_ += rateContribution;

    if (ratesInitialized_ && i < binaryBreakupRates_.size())
    {
        binaryBreakupRates_[i] += rateContribution;
    }
}


// ************************************************************************* //
