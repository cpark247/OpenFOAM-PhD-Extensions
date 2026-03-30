/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "sigmaModel.H"
#include "fvConstraints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
dimensionedScalar sigmaModel<BasicMomentumTransportModel>::cI
(
   const volTensorField& D
) const
{
  #include "cI.H"
}

template<class BasicMomentumTransportModel>
volScalarField sigmaModel<BasicMomentumTransportModel>::cD
(
     const volTensorField& D
) const
{
  #include "cD.H"
}




template<class BasicMomentumTransportModel>
tmp<volScalarField> sigmaModel<BasicMomentumTransportModel>::k
(
    const volTensorField& D
) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("k", this->alphaRhoPhi_.group()),
                this->runTime_.timeName(),
                this->mesh_
            ),
            cI(D)*sqr(this->delta())*magSqr(dev(symm(D)))
        )
    );
}


template<class BasicMomentumTransportModel>
void sigmaModel<BasicMomentumTransportModel>::correctNut()
{
    DSigma_ = cD(fvc::grad(this->U_));
    Info << "cSigma = " << cSigma_ << endl;
    this->nut_ = DSigma_*sqr(cSigma_)*sqr(this->delta());	//eq. 7 in TSFP_BayaToda
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
sigmaModel<BasicMomentumTransportModel>::sigmaModel
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
    LESeddyViscosity<BasicMomentumTransportModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
    ),
    cSigma_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cSigma",
            this->coeffDict_,
            1.5
        )
    ),
    epsilonMin_("epsilonMin",sqr(dimVelocity)/dimTime, SMALL),
    filterPtr_(LESfilter::New(this->mesh_, this->coeffDict())),
    filter_(filterPtr_()),
    DSigma_
    (
    	IOobject
    	(
        IOobject::groupName("DSigma_", this->alphaRhoPhi_.group()),
        this->runTime_.timeName(),
        this->mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    	),
    	this->mesh_,
      dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0, 0, 0),0.0)
    )
{
    // In old version bound(k_,k0_); bound(epsilon_, epsilonMin_); here. Might be in basic correct
    correctNut();
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool sigmaModel<BasicMomentumTransportModel>::read()
{
    if (LESeddyViscosity<BasicMomentumTransportModel>::read())
    {
        filter_.read(this->coeffDict());
        cSigma_.readIfPresent(this->coeffDict());
        return true;
    }

    return false;
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> sigmaModel<BasicMomentumTransportModel>::epsilon() const
{
    volScalarField k(this->k(fvc::grad(this->U_)));

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            2.0*this->nuEff()*magSqr(dev(symm(fvc::grad(this->U_))))
        )
    );
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> sigmaModel<BasicMomentumTransportModel>::omega() const
{
    volScalarField k(this->k(fvc::grad(this->U_)));
    volScalarField epsilon(this->epsilon());

    return tmp<volScalarField>
    (
	 	new volScalarField
		(
        	IOobject
        	(
        	    IOobject::groupName("omega", this->alphaRhoPhi_.group()),
        	    this->runTime_.timeName(),
        	    this->mesh_
        	),
        	epsilon/(0.09*k)
		)
    );
}



template<class BasicMomentumTransportModel>
void sigmaModel<BasicMomentumTransportModel>::correct()
{
    LESeddyViscosity<BasicMomentumTransportModel>::correct();
    // In old version bound(k_,k0_); bound(epsilon_, epsilonMin_); here. Might be in basic correct
    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
