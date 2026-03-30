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

#include "dynSigmaModel.H"
#include "fvConstraints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
dimensionedScalar dynSigmaModel<BasicMomentumTransportModel>::cI
(
   const volTensorField& D
) const
{
  #include "cI.H"
}

template<class BasicMomentumTransportModel>
volScalarField dynSigmaModel<BasicMomentumTransportModel>::cD
(
     const volTensorField& D
) const
{
  #include "cD.H"
}


template<class BasicMomentumTransportModel>
scalar dynSigmaModel<BasicMomentumTransportModel>::cS
(
     const volTensorField& D, volScalarField& DSigma
) const
{
    tensorField g(D.size());
    scalar cSigma_ = 0;

    scalarField SVS(D.size());
    const volSymmTensorField s
    (
    	dev(symm(D))
    );

    const volSymmTensorField MM // according to TSFP_BayaToda, bottom left of page 2
    (
	sqr(this->delta())*(filter_(DSigma*(s)) - 4*filter_(DSigma)*filter_(s))
    );

    dimensionedScalar MMMM
    (
       	"MMMM",
       	dimensionSet(0, 4, -4, 0, 0, 0, 0),
       	scalar(0)
    );

    // SVS-averaging, according to eq. 9 in TSFP_BayaToda
    scalar sumSVS = 0;
    forAll(D,celli)
    {
	    g[celli] = dev(symm(D[celli]&D[celli]));
    	//SVS[celli]  = pow((g[celli]&&g[celli]),(3./2.))/(pow((g[celli]&&g[celli]),(3./2.)) + pow((s[celli]&&s[celli]),(3.)));
    	SVS[celli]  = pow((g[celli]&&g[celli]),(3./2.))/(pow((g[celli]&&g[celli]),(3./2.)) + pow((s[celli]&&s[celli]),(3.))+SMALL);
        MMMM.value() += SVS[celli]*(MM[celli]&&MM[celli]);
        sumSVS += SVS[celli];
    }

    reduce(sumSVS,sumOp<scalar>());
    reduce(MMMM.value(),sumOp<scalar>());

    MMMM.value() = MMMM.value()/(sumSVS+SMALL);
    Info << "MMMM = " << MMMM << endl;

    dimensionedScalar LLMM
    (
       	"LLMM",
       	dimensionSet(0, 4, -4, 0, 0, 0, 0),
       	scalar(0)
    );
    if (MMMM.value() > VSMALL)
    {
        const volTensorField LL
	(
		filter_(this->U_*this->U_)-sqr(filter_(this->U_))
		//dev(filter_(sqr(U()))-sqr(filter_(U())))
	);

	forAll(D,celli)
    	{
		LLMM.value() += SVS[celli]*(LL[celli]&&MM[celli]);
    	}
    // summation over all processor domains
    reduce(LLMM.value(),sumOp<scalar>());

	// eq. 10 in TSFP_BayaToda
	    if ((0.5*LLMM/(MMMM*sumSVS)).value() >= 0.0)
	    {
		cSigma_= (0.5*LLMM/(MMMM*sumSVS)).value();
	    }
	return (cSigma_);
    }
    else
    {
        return 0.0;
    }
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> dynSigmaModel<BasicMomentumTransportModel>::k
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
void dynSigmaModel<BasicMomentumTransportModel>::correctNut()
{
    DSigma_ = cD(fvc::grad(this->U_));
    scalar cSigma = cS(fvc::grad(this->U_),DSigma_);
    Info << "avg(DSigma) = " << average(DSigma_) << endl;
    Info << "cSigma = " << sqrt(cSigma) << endl;
    this->nut_ = DSigma_*cSigma*sqr(this->delta());	//eq. 7 in TSFP_BayaToda // cSigma sqrt or not??? Old without sqrt sigmaModel -> sqrt(cSigma_)
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
dynSigmaModel<BasicMomentumTransportModel>::dynSigmaModel
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
    // In old version bound(k_,k0_); bound(epsilon_, epsilonMin_); here.
    correctNut();
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool dynSigmaModel<BasicMomentumTransportModel>::read()
{
    if (LESeddyViscosity<BasicMomentumTransportModel>::read())
    {
        filter_.read(this->coeffDict());
        //cSigma_.readIfPresent(this->coeffDict());
        return true;
    }

    return false;
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> dynSigmaModel<BasicMomentumTransportModel>::epsilon() const
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
tmp<volScalarField> dynSigmaModel<BasicMomentumTransportModel>::omega() const
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
void dynSigmaModel<BasicMomentumTransportModel>::correct()
{
    LESeddyViscosity<BasicMomentumTransportModel>::correct();
    // In old version bound(k_,k0_); bound(epsilon_, epsilonMin_); here. Might be in basic correct
    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
