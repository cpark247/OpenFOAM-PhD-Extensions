/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
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

#include "kOmegaSSTSAS_ext.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> kOmegaSSTSAS_ext<BasicMomentumTransportModel>::Qsas
(
    const volScalarField::Internal& S2,
    const volScalarField::Internal& gamma,
    const volScalarField::Internal& beta
) const
{
    volScalarField::Internal L
    (
        sqrt(this->k_())/(pow025(this->betaStar_)*this->omega_())
    );

    volScalarField::Internal Lvk
    (
        max
        (
            kappa_*sqrt(S2)
           /(
                mag(fvc::laplacian(this->U_))()()
              + dimensionedScalar
                (
                    "rootVSmall",
                    dimensionSet(0, -1, -1, 0, 0),
                    rootVSmall
                )
            ),
            Cs_*sqrt(kappa_*zeta2_/(beta/this->betaStar_ - gamma))*delta()()
        )
    );

    // Declaration of QsasSource here as tmp
    tmp<volScalarField::Internal> QsasSource = 
        this->alpha_()*this->rho_()
       *min
        (
            max
            (
                zeta2_*kappa_*S2*pow(L/Lvk, nExponent_)
              - (2*C_/sigmaPhi_)*this->k_()
               *max
                (
                    magSqr(fvc::grad(this->omega_)()())/sqr(this->omega_()),
                    magSqr(fvc::grad(this->k_)()())/sqr(this->k_())
                ),
                dimensionedScalar(dimensionSet(0, 0, -2, 0, 0), 0)
            ),
            this->omega_()/(0.1*this->omega_.time().deltaT())
        );
    
    // Use const_cast to modify output fields from const function
    const_cast<volScalarField&>(L_).primitiveFieldRef() = L;    
    const_cast<volScalarField&>(Lvk_).primitiveFieldRef() = Lvk;
	const_cast<volScalarField&>(L_by_Lvk_).primitiveFieldRef() = L/max(Lvk,dimensionedScalar(dimless, 1.0e-10));
    const_cast<volScalarField&>(Qsas_).primitiveFieldRef() = QsasSource() / this->alpha_() / this->rho_(); 

    // Return the matrix with QsasSource() dereferenced
    return fvm::Su(QsasSource(), this->omega_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
kOmegaSSTSAS_ext<BasicMomentumTransportModel>::kOmegaSSTSAS_ext
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& type
)
:
    kOmegaSST_ext<BasicMomentumTransportModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
    ),
    nExponent_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "nExponent",
            this->coeffDict_,
            2.0
        )
    ),
    Cs_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cs",
            this->coeffDict_,
            0.11
        )
    ),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    zeta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "zeta2",
            this->coeffDict_,
            3.51
        )
    ),
    sigmaPhi_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaPhi",
            this->coeffDict_,
            2.0/3.0
        )
    ),
    C_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C",
            this->coeffDict_,
            2
        )
    ),

    delta_
    (
        LESdelta::New
        (
            IOobject::groupName("delta", alphaRhoPhi.group()),
            *this,
            this->coeffDict_
        )
    ),
    L_
    (
        IOobject
        (
            IOobject::groupName("L_", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("L_", dimLength, 0.0)
    ),
    Lvk_
    (
        IOobject
        (
            IOobject::groupName("Lvk_", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("Lvk_", dimLength, 0.0)
    ),
	L_by_Lvk_
	(
        IOobject
        (
            IOobject::groupName("L_by_Lvk_", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("L_by_Lvk_", dimless, 0.0)
    ),
    Qsas_
    (
        IOobject
        (
            IOobject::groupName("Qsas_", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("Qsas", dimless/dimTime/dimTime, 0.0)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool kOmegaSSTSAS_ext<BasicMomentumTransportModel>::read()
{
    if (kOmegaSST_ext<BasicMomentumTransportModel>::read())
    {
        nExponent_.readIfPresent(this->coeffDict());
        Cs_.readIfPresent(this->coeffDict());
        kappa_.readIfPresent(this->coeffDict());
        sigmaPhi_.readIfPresent(this->coeffDict());
        zeta2_.readIfPresent(this->coeffDict());
        C_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
