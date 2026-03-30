/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::RASModels::BSL_EARSM

Group
    grpRASTurbulence

See also
    Foam::kOmegaSST

SourceFiles
    BSL_EARSM.C

\*---------------------------------------------------------------------------*/

#ifndef BSL_EARSM_H
#define BSL_EARSM_H

#include "kOmegaSSTBase.H"
#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                          Class BSL_EARSM Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class BSL_EARSM
:
    public Foam::kOmegaSST
    <
        eddyViscosity<RASModel<BasicMomentumTransportModel>>,
        BasicMomentumTransportModel
    >
{

protected:
            dimensionedScalar alphaK1_;
            dimensionedScalar alphaOmega2_;
            // dimensionedScalar gamma1_;
            // dimensionedScalar gamma2_;
            // dimensionedScalar Prt_;
            dimensionedScalar kappa_;
            dimensionedScalar C1_;
            dimensionedScalar A1_;

            volScalarField tau_;
            volTensorField gradU_;
            volSymmTensorField S_;
            volTensorField W_;
            volScalarField IIs_;
            volScalarField IIo_;
            volScalarField IV_;
            dimensionedScalar C1prime_;
            volScalarField P1_;
            volScalarField P2_;
            volScalarField N_;
            volScalarField Q_;
            volScalarField Q1_;
            volScalarField Beta1_;
            volScalarField Beta2_;
            volScalarField Beta3_;
            volScalarField Beta4_;
            volScalarField Beta6_;
            volScalarField Beta9_;

            volScalarField check0_;
            volScalarField check1_;
            volScalarField check2_;
            volScalarField check3_;
            volScalarField check4_;
            volScalarField check5_;

            volSymmTensorField aij_;

            volSymmTensorField tauij_;
            volSymmTensorField Nij_;
public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("BSL_EARSM");


    // Constructors

        //- Construct from components
        BSL_EARSM
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            // const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~BSL_EARSM()
    {}
    //- Return the source term for the momentum equation
    virtual tmp<volSymmTensorField> R() const;

    //- Return the source term for the momentum equation
    virtual tmp<volSymmTensorField> devRhoReff() const;

    //- Return the source term for the momentum equation
    virtual tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const;

    //- Solve the turbulence equations and correct the turbulence viscosity
    virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "BSL_EARSM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
