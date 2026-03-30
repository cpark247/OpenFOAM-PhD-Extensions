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

#include "PrinceBlanchMom.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "phaseDynamicMomentumTransportModel.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace coalescenceMomModels
{
    defineTypeNameAndDebug(PrinceBlanchMom, 0);
    addToRunTimeSelectionTable
    (
        coalescenceMomModel,
        PrinceBlanchMom,
        dictionary
    );
}
}
}

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::coalescenceMomModels::PrinceBlanchMom::
PrinceBlanchMom
(
    const populationBalanceMomModel& popBal,
    const dictionary& dict
)
:
    coalescenceMomModel(popBal, dict),
    C1_(dimensionedScalar::lookupOrDefault("C1", dict, dimless, 0.356)),
    h0_(dimensionedScalar::lookupOrDefault("h0", dict, dimLength, 1e-4)),
    hf_(dimensionedScalar::lookupOrDefault("hf", dict, dimLength, 1e-8)),
    turbulence_(dict.lookup("turbulence")),
    buoyancy_(dict.lookup("buoyancy")),
    laminarShear_(dict.lookup("laminarShear"))
{
    if (laminarShear_)
    {
        shearStrainRate_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "shearStrainRate",
                    popBal.time().timeName(),
                    popBal.mesh()
                ),
                popBal.mesh(),
                dimensionedScalar
                (
                    "shearStrainRate",
                    dimVelocity/dimLength,
                    Zero
                )
            )
        );
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::coalescenceMomModels::PrinceBlanchMom::precompute()
{
    if (laminarShear_)
    {
        shearStrainRate_() =
            sqrt(2.0)*mag(symm(fvc::grad(popBal_.continuousPhase().U())));
    }
}

void Foam::diameterModels::coalescenceMomModels::PrinceBlanchMom::
addToCoalescenceRate
(
    volScalarField& coalescenceRate,
    const volScalarField& di,
    const volScalarField& dj,
    const label i,
    const label j
)
{
    const phaseModel& continuousPhase = popBal_.continuousPhase();
    const uniformDimensionedVectorField& g =
        popBal_.mesh().lookupObject<uniformDimensionedVectorField>("g");

	volScalarField rij(1.0/(1.0/di + 1.0/dj));

    const volScalarField collisionEfficiency
    (
        exp
        (
          - sqrt
            (
                pow3(rij)*continuousPhase.rho()
               /(16.0*popBal_.sigmaWithContinuousPhase(popBal_.refPhase()))
            )
           *log(h0_/hf_)
           *cbrt(popBal_.continuousTurbulence().epsilon())/pow(rij, 2.0/3.0)
        )
    );

    if (turbulence_)
    {
        coalescenceRate +=
            (
                C1_*pi*sqr(di + dj)
               *cbrt(popBal_.continuousTurbulence().epsilon())
               *sqrt(pow(di, 2.0/3.0) + pow(dj, 2.0/3.0))
            )
           *collisionEfficiency;
    }

    if (buoyancy_)
    {
        volScalarField Sij(pi/4.0*sqr(di + dj));

        coalescenceRate +=
            (
                Sij
               *mag
                (
                    sqrt
                    (
                        2.14*popBal_.sigmaWithContinuousPhase(popBal_.refPhase())
                       /(continuousPhase.rho()*di) + 0.505*mag(g)*di
                    )
                  - sqrt
                    (
                        2.14*popBal_.sigmaWithContinuousPhase(popBal_.refPhase())
                       /(continuousPhase.rho()*dj) + 0.505*mag(g)*dj
                    )
                )
            )
           *collisionEfficiency;
    }

    if (laminarShear_)
    {
        coalescenceRate +=
            pow3(di + dj)/6
           *shearStrainRate_()*collisionEfficiency;
    }
}


// ************************************************************************* //
