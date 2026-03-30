/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::LESModels::sigmaModel

Group
    grpLESTurbulence

Description
Sigma model,
according to Franck Nicoud et.al, 2011, Physics of Fluids, vol. 23, "Using singular values to build a subgrid-scale model for large eddy simulations"
and Baya Toda Hubert et.al, "A dynamic procedure for advanced subgrid-scale models and wall-bounded flows"

ATTENTION: D = grad(U)

Algebraic eddy viscosity SGS model founded on the assumption that
local equilibrium prevails.
Thus,
\verbatim
    B = 2/3*k*I - 2*nuSgs*dev(symm(D))
    Beff = 2/3*k*I - 2*nuEff*dev(symm(D))

where

    k = cI*delta^2*||dev(symm(D)||^2
  epsilon = nuEff*||dev(symm(D)||^2 (Diss. Brussies, 2014)
    nuSgs = cDSigma*delta^2*DSigma
    nuEff = nuSgs + nu

In the dynamic version of the choric  Smagorinsky model
the coefficients cI and cD are calculated during the simulation,

    cI=<K*m>/<m*m>

and

    cD=1/2*<L.M>/<M.M>,

where

    K = 0.5*(F(U.U) - F(U).F(U))
    m = delta^2*(4*||F(dev(symm(D)))||^2 - F(||dev(symm(D))||^2))
    L = dev(F(U*U) - F(U)*F(U))
    M = delta^2*(F(DSigma*dev(symm(D))) - 4*F(DSigma)*F(dev(symm(D))))

The averaging <...> is over the whole domain, i.e. homogeneous turbulence
is assumed, containing SVS weighting (cf. paper of Baya Toda Hubert et.al)
\endverbatim

SourceFiles
sigmaModel.C

\*---------------------------------------------------------------------------*/

#ifndef dynSigmaModel_H
#define dynSigmaModel_H

#include "LESModel.H"
#include "LESeddyViscosity.H"
#include "LESfilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class WALE Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class dynSigmaModel
:
    public LESeddyViscosity<BasicMomentumTransportModel>
{
    // Private Member Functions

        //- No copy construct
        dynSigmaModel(const dynSigmaModel&) = delete;

        //- No copy assignment
        void operator=(const dynSigmaModel&) = delete;


protected:

    // Protected data

        //dimensionedScalar cSigma_;
        dimensionedScalar epsilonMin_;
        autoPtr<LESfilter> filterPtr_;
        LESfilter& filter_;
        volScalarField DSigma_;
    // Protected Member Functions

        //- Return the deviatoric symmetric part of the square of the given
        //  velocity gradient field

        scalar cS(const volTensorField& D, volScalarField& DSigma) const;
  	    dimensionedScalar cI(const volTensorField& D) const;

        //- Calculate coefficients cD, cI from filtering velocity field
        volScalarField cD(const volTensorField& D) const;

        //- Return SGS kinetic energy
        //  calculated from the given velocity gradient
        tmp<volScalarField> k(const volTensorField& D) const;

        //- Update the SGS eddy-viscosity
        virtual void correctNut();


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("dynSigmaModel");


    // Constructors

        //- Construct from components
        dynSigmaModel
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& type = typeName
        );


    //- Destructor
    virtual ~dynSigmaModel() = default;


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k(fvc::grad(this->U_));
        }

        //- Return sub-grid dissipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Return sub-grid specific dissipation rate
        virtual tmp<volScalarField> omega() const;

        //- Correct Eddy-Viscosity and related properties
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dynSigmaModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
