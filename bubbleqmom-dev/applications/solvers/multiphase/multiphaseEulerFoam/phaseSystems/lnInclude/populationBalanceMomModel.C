/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "populationBalanceMomModel.H"
#include "coalescenceMomModel.H"
#include "breakupMomModel.H"
#include "phaseSystem.H"
#include "surfaceTensionModel.H"
#include "fvm.H"
#include "fvcDdt.H"
#include "phaseDynamicMomentumTransportModel.H"
#include "bound.H"
// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //
#define SMALLVOLQMOM 1e-60 // the SMALL from OpenFOAM is too large 

void
Foam::diameterModels::populationBalanceMomModel::registerPhasesAndQmom()
{
    forAll(fluid_.phases(), phasei)
    {
	Info << "Foam::diameterModels::populationBalanceMomModel::registerPhasesAndQmom:  phase " << phasei << " -- " << fluid_.phases()[phasei].name() << "\n";
	Info << "Foam::diameterModels::populationBalanceMomModel::registerPhasesAndQmom:  the continuousPhase is  " << continuousPhase_.name() << "\n";
        if (isA<qmomVelocity>(fluid_.phases()[phasei].dPtr()()))
        {
            const qmomVelocity& qmomVelGroup =
            refCast<const qmomVelocity>(fluid_.phases()[phasei].dPtr()());

            if (qmomVelGroup.popBalName() == this->name())
            {
		Info << "Foam::diameterModels::populationBalanceMomModel::registerPhasesAndQmom:  added to phase " << phasei << " -- " << fluid_.phases()[phasei].name() << "  popBalName of " << qmomVelGroup.popBalName() << "\n";
                qmomVel_.reset(&const_cast<qmomVelocity&>(qmomVelGroup));

        	quadrature_.reset
        	(
        	    new monoKineticQuadratureApproximation
        	    (
				 	fluid_.phases()[phasei].name(),
					fluid_.mesh(),
					"RPlus"
        	    )
        	);
		//this->calcAlphas();
		//qmomVel_().updateSauterMeanDiam(dsm_); TODO
		// qmom considers only one refrence phase
		refPhaseIndex_ = phasei;
            }
        }
    }
}


////not required here since all bubbles are in one qmom object, but it might be helful later for phase transfer, see sectional method
//void Foam::diameterModels::populationBalanceMomModel::createPhasePairs()
//{
//}


void Foam::diameterModels::populationBalanceMomModel::precompute()
{

    forAll(coalescence_, model)
    {
        coalescence_[model].precompute();
    }

    forAll(breakup_, model)
    {
        breakup_[model].precompute();

        breakup_[model].dsdPtr()().precompute();
    }

}


void Foam::diameterModels::populationBalanceMomModel::
evalCoalescenceST()
{
    const label nMoments = quadrature_().nMoments();
    const PtrList<volScalarNode>& nodes = quadrature_().nodes();
	
    if (coalescence_.size() != 0)
    {
        for(label k = 0; k < nMoments; k++)
        {
	    forAll(nodes, nodei)
	    {
	        const volScalarNode& node1 = nodes[nodei];
	    	const dimensionedScalar verySmallNodeVal("ApproxZero", node1.primaryAbscissae()[0].dimensions(), SMALLVOLQMOM);
	        const volScalarField abscissa1(Foam::max(node1.primaryAbscissae()[0], verySmallNodeVal));
	        const volScalarField n1(node1.n(node1.primaryWeight(), abscissa1));
	        const volScalarField d1(node1.d(abscissa1));

	        forAll(nodes, nodej)
	        {
		    const volScalarNode& node2 = nodes[nodej];
		    const volScalarField abscissa2(Foam::max(node2.primaryAbscissae()[0], verySmallNodeVal));
	    
		    coalescenceRate_() = Zero;

		    forAll(coalescence_, model)
		    {
		        coalescence_[model].addToCoalescenceRate
		        (
		            coalescenceRate_(),
		    			d1,
		    			node2.d(abscissa2),
		            nodei,
		            nodej
		        );
		    }

		    // update moment source term, mass based formulation
		    Sm_[k] += 0.5*n1
		       *(
		            node2.n(node2.primaryWeight(), abscissa2)
		           *(
		                pow
		                (
		                    abscissa1 + abscissa2,
		                    k
		                )
		              - pow(abscissa1, k)
		              - pow(abscissa2, k)
		            )
		        )*coalescenceRate_();
	    	}
	    }
    	}
    }
}


void Foam::diameterModels::populationBalanceMomModel::
evalBreakupST()
{
    const label nMoments = quadrature_().nMoments();
    const PtrList<volScalarNode>& nodes = quadrature_().nodes();
	
    if (breakup_.size() != 0)
    {
        forAll(breakup_, model)
        {
        	for(label k = 0; k < nMoments; k++)
        	{
    			forAll(nodes, nodei)
    			{
        	    	const volScalarNode& node = nodes[nodei];
    				const dimensionedScalar verySmallNodeVal("ApproxZero", node.primaryAbscissae()[0].dimensions(), SMALLVOLQMOM);
        	    	const volScalarField abscissa(Foam::max(node.primaryAbscissae()[0], verySmallNodeVal));
    		    	const volScalarField d(node.d(abscissa));

					breakupRate_() *= 0.;
        	   		    breakup_[model].setBreakupRate
        	   		    (
        	   		        breakupRate_(),
							d,
							nodei
        	   		    );
						breakup_[model].dsdPtr()().calcIntegralm
        	   		    (
        	   		        breakDsdIntegr_[k],
							abscissa,
							nodei,
							k
        	   		    );
					// update moment source term
					Sm_[k] += node.n(node.primaryWeight(), abscissa)
							   *breakupRate_()
        	    			   *(
									breakDsdIntegr_[k]             // birth
        	    			        -pow
        	    			        (
        	    			            abscissa,
        	    			            k
        	    			        )                              // death
        	    			    );
       			}
    		}
        }
	}
}


void Foam::diameterModels::populationBalanceMomModel::calcAlphas()
{
    alphas_() = Zero;

    const phaseModel& phase = fluid_.phases()[refPhaseIndex_];
    const volScalarField& rho = phase.thermo().rho();
    const volScalarField& alpha = phase;
	
    // total alpha
    alphas_() == quadrature_().moments()[1]/rho;
    // node - alphas
    if (quadrature_().nodes().size() == 1)
    {
        nodeAlphas_[0] = alpha;
        ds_[0] =
            Foam::min
            (
                Foam::max
                (
                    Foam::pow
                    (
                        quadrature_().nodes()[0].primaryAbscissae()[0]*6.0
                       /(rho*Foam::constant::mathematical::pi)
                      + dimensionedScalar("smallVolume", dimVolume, SMALLVOLQMOM),
                        1.0/3.0
                    ),
                    minD_
                ),
                maxD_
            );
        dsm_() = ds_[0];
    }
    else
    {
        dsm_() = dimensionedScalar("zero", dimLength, 0.0);
        volScalarField scale
        (
            alpha
           /Foam::max
            (
                quadrature_().moments()[1]/rho,
                residualAlpha_
            )
        );

        forAll(quadrature_().nodes(), nodei)
        {
            const volScalarNode& node = quadrature_().nodes()[nodei];

            // Set alpha values such that the moment.1 is equal to the bounded alpha
            nodeAlphas_[nodei] = node.primaryWeight()*node.primaryAbscissae()[0]/rho*scale;
            nodeAlphas_[nodei].max(0);
            nodeAlphas_[nodei].min(1);

            //  Calculate bubble diameter based on bubble mass (abscissa)
            ds_[nodei] =
                Foam::min
                (
                    Foam::max
                    (
                        Foam::pow
                        (
                            node.primaryAbscissae()[0]*6.0
                           /(rho*Foam::constant::mathematical::pi)
                          + dimensionedScalar("smallVolume", dimVolume, SMALLVOLQMOM),
                            1.0/3.0
                        ),
                        minD_
                    ),
                    maxD_
                );
            dsm_() += nodeAlphas_[nodei]*ds_[nodei];
        }

        dsm_() /= Foam::max(alpha, residualAlpha_);
        dsm_().max(minD_);
    }
    qmomVel_().updateSauterMeanDiam(dsm_());

    // write out nodes and weights
    forAll(quadrature_().nodes(), nodei)
    {
        const volScalarNode& node = quadrature_().nodes()[nodei];
	QMOMNode_[nodei]   = node.primaryAbscissae()[0];
	QMOMWeight_[nodei] = node.primaryWeight();
    }
}


void Foam::diameterModels::populationBalanceMomModel::scaleMomentsToAlpha()
{
    const phaseModel& phase = fluid_.phases()[refPhaseIndex_];
    const volScalarField& rho = phase.thermo().rho();
    const volScalarField& alphaRef = phase;
    volScalarField corrFactor  = (alphaRef * rho / Foam::max(quadrature_().moments()[1], residualAlpha_*rho));
    forAll(quadrature_().moments(), mi)
    {
    	quadrature_().moments()[mi] *= corrFactor;
    }
}


//void Foam::diameterModels::populationBalanceMomModel::calcVelocity()
//{
//    U_() = Zero;
//    const phaseModel& phase = fluid_.phases()[refPhaseIndex_];
//    U_() = phase.U();
//}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::populationBalanceMomModel::populationBalanceMomModel
(
    const phaseSystem& fluid,
    const word& name,
    HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash>& pDmdt
)
:
    regIOobject
    (
        IOobject
        (
            name,
            fluid.mesh().time().constant(),
            fluid.mesh()
        )
    ),
    fluid_(fluid),
    refPhaseIndex_(-1),
    pDmdt_(pDmdt),
    mesh_(fluid.mesh()),
    name_(name),
    dict_
    (
        fluid.subDict("populationBalanceCoeffs").subDict(name_)
    ),
    pimple_(mesh_.lookupObject<pimpleControl>("solutionControl")),
    continuousPhase_
    (
        mesh_.lookupObject<phaseModel>
        (
            IOobject::groupName("alpha", dict_.lookup("continuousPhase"))
        )
    ),
    quadrature_(),
    qmomVel_(),
    Sm_(),
    coalescence_
    (
        dict_.lookup("coalescenceModels"),
        coalescenceMomModel::iNew(*this)
    ),
    coalescenceRate_(),
    breakup_
    (
        dict_.lookup("breakupModels"),
        breakupMomModel::iNew(*this)
    ),
    breakupRate_(),
    breakDsdIntegr_(),
    alphas_(),
    nodeAlphas_(),
    ds_(),
    QMOMNode_(),
    QMOMWeight_(),
    dsm_(),
    maxD_("maxD", dimLength, dict_.lookupOrDefault<scalar>("maxD", 1e-2)),
    minD_("minD", dimLength, dict_.lookupOrDefault<scalar>("minD", 1e-5)),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        dict_.lookupOrDefault<scalar>("residualAlpha", 1e-6)
    )
    //U_()
{
	Info << "Make Population Balance Moment Model name_ = " << name_ << "\n";
    this->registerPhasesAndQmom();

    //this->createPhasePairs();

    for(int nodei=0; nodei<quadrature_().nodes().size(); nodei++)
    {
        nodeAlphas_.append
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "alphaNode",
                        IOobject::groupName
                        (
                            name_,
                            Foam::name(nodei)
                        )
                    ),
                    fluid.mesh().time().timeName(),
                    fluid.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluid.mesh(),
                dimensionedScalar("alpha", dimless, 0.0)
            )
        );
        ds_.append
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "diQmom",
                        IOobject::groupName
                        (
                            name_,
                            Foam::name(nodei)
                        )
                    ),
                    fluid.mesh().time().timeName(),
                    fluid.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluid.mesh(),
                dimensionedScalar("diqmom", dimLength, 0.0)
            )
        );
        QMOMNode_.append
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "xiQmom",
                        IOobject::groupName
                        (
                            name_,
                            Foam::name(nodei)
                        )
                    ),
                    fluid.mesh().time().timeName(),
                    fluid.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluid.mesh(),
                dimensionedScalar("xiqmom", dimMass, 0.0)
            )
        );
        QMOMWeight_.append
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "wQmom",
                        IOobject::groupName
                        (
                            name_,
                            Foam::name(nodei)
                        )
                    ),
                    fluid.mesh().time().timeName(),
                    fluid.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluid.mesh(),
                dimensionedScalar("wqmom", quadrature_().moments()[0].dimensions(), 0.0)
            )
        );

    const label nMoments = quadrature_().nMoments();
    for(int momi=0; momi<nMoments; momi++)
    {
 	   Sm_.append
	    (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "mST",
                        IOobject::groupName
                        (
                            name_,
                            Foam::name(momi)
                        )
                    ),
                    fluid.mesh().time().timeName(),
                    fluid.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluid.mesh(),
                dimensionedScalar("momST", quadrature_().moments()[momi].dimensions()/dimTime, 0.0)
            )
		);
 	   breakDsdIntegr_.append
	    (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "breakDsdIntegral",
                        IOobject::groupName
                        (
                            name_,
                            Foam::name(momi)
                        )
                    ),
                    fluid.mesh().time().timeName(),
                    fluid.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluid.mesh(),
                dimensionedScalar("breakDsdIntegral", pow(dimMass,momi), 0.0)
            )
		);
	}
    }

    if (coalescence_.size() != 0)
    {
        coalescenceRate_.reset
        (
            new volScalarField
            (
                IOobject
                (
                     "coalescenceRate",
                     fluid.time().timeName(),
                     fluid.mesh()
                ),
                fluid.mesh(),
                dimensionedScalar("coalescenceRate", dimVolume/dimTime, Zero)
            )
        );
    }

    if (breakup_.size() != 0)
    {
        breakupRate_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "breakupRate",
                    fluid.time().timeName(),
                    fluid.mesh()
                ),
                fluid.mesh(),
                dimensionedScalar("breakupRate", inv(dimTime), Zero)
            )
        );
    }

    alphas_.set
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("alpha", this->name()),
                fluid_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, Zero)
        )
    );
    
    dsm_.set
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("d", this->name()),
                fluid_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimLength, Zero)
        )
    );
    
    //U_.set
    //(
    //    new volVectorField
    //    (
    //        IOobject
    //        (
    //            IOobject::groupName("U", this->name()),
    //            fluid_.time().timeName(),
    //            mesh_,
    //            IOobject::NO_READ,
    //            IOobject::NO_WRITE
    //        ),
    //        mesh_,
    //        dimensionedVector(dimVelocity, Zero)
    //    )
    //);
	
    //const dictionary& pimpleDict =
    //    fluid_.mesh().solutionDict().subDict("PIMPLE");
    //label nCorrectors = pimpleDict.lookupOrDefault<label>("nFluxCorrectors", 0);
    //if (nCorrectors > 0)
    //{
    //    word patchName
    //    (
    //        pimpleDict.lookupOrDefault
    //        (
    //            "corrPatch",
    //            U_.boundaryField()[0].patch().name()
    //        )
    //    );
    //    wordList boundaries(U_.boundaryField().size(), "zeroGradient");
    //    forAll(boundaries, patchi)
    //    {
    //        if (U_.boundaryField()[patchi].patch().name() == patchName)
    //        {
    //            boundaries[patchi] = "fixedFace";
    //        }
    //    }

    //    corr_ = tmp<volScalarField>
    //    (
    //        new volScalarField
    //        (
    //            IOobject
    //            (
    //                "corr",
    //                fluid_.mesh().time().timeName(),
    //                fluid_.mesh(),
    //                IOobject::NO_READ,
    //                IOobject::NO_WRITE
    //            ),
    //            fluid_.mesh(),
    //            dimensionedScalar("0", sqr(dimLength)/dimTime, 0.0),
    //            boundaries
    //        )
    //    );
    //}
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::populationBalanceMomModel::~populationBalanceMomModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::diameterModels::populationBalanceMomModel>
Foam::diameterModels::populationBalanceMomModel::clone() const
{
    notImplemented("populationBalance::clone() const");
    return autoPtr<populationBalanceMomModel>(nullptr);
}


bool Foam::diameterModels::populationBalanceMomModel::writeData(Ostream& os) const
{
    return os.good();
}



const Foam::tmp<Foam::volScalarField>
Foam::diameterModels::populationBalanceMomModel::sigmaWithContinuousPhase
(
    const phaseModel& dispersedPhase
) const
{
    return
        fluid_.lookupSubModel<surfaceTensionModel>
        (
            phasePair(dispersedPhase, continuousPhase_)
        ).sigma();
}


const Foam::phaseCompressible::momentumTransportModel&
Foam::diameterModels::populationBalanceMomModel::continuousTurbulence() const
{
    return
        mesh_.lookupObject<phaseCompressible::momentumTransportModel>
        (
            IOobject::groupName
            (
                momentumTransportModel::typeName,
                continuousPhase_.name()
            )
        );
}


const phaseModel&
Foam::diameterModels::populationBalanceMomModel::refPhase() const
{
    return (fluid_.phases()[refPhaseIndex_]);
}


void Foam::diameterModels::populationBalanceMomModel::averageTransport()
{
    const label nNodes = quadrature_().nodes().size();
    const phaseModel& phase = fluid_.phases()[refPhaseIndex_];
    // Correct mean flux
    const PtrList<surfaceScalarNode>& nodesOwn = quadrature_().nodesOwn();
    const PtrList<surfaceScalarNode>& nodesNei = quadrature_().nodesNei();
    dimensionedScalar zeroPhi("zero", phase.phi()().dimensions(), 0.0);
    dimensionedScalar globalDt("dtToSolve",quadrature_().moments()[0].mesh().time().deltaT());
    surfaceScalarField phi(phase.phi());

	// may include if (corr_.valid()) in +1010 OQbMM-polydispersePhaseModel.C
	
    quadrature_().interpolateNodes();

    // Mean moment advection
    Info<< "Transporting moments with average velocity" << endl;
    forAll(quadrature_().moments(), mEqni)
    {
        volScalarField& m = quadrature_().moments()[mEqni];

        volScalarField meanDivUbMp
        (
            IOobject
            (
                "meanDivUbMp",
                fluid_.mesh().time().timeName(),
                fluid_.mesh()
            ),
            fluid_.mesh(),
            dimensionedScalar("zero", m.dimensions()/dimTime, Zero)
        );
        for (label nodei = 0; nodei < nNodes; nodei++)
        {
            // Update average size moment flux
            surfaceScalarField aFluxMp
            (
                "aFluxMp",
                nodesNei[nodei].primaryWeight()
               *(
                    pow
                    (
                        nodesNei[nodei].primaryAbscissae()[0],
                        mEqni
                    )
                )*Foam::min(phi, zeroPhi)
              + nodesOwn[nodei].primaryWeight()
               *pow
                (
                    nodesOwn[nodei].primaryAbscissae()[0],
                    mEqni
                )*Foam::max(phi, zeroPhi)
            );

            meanDivUbMp += fvc::surfaceIntegrate(aFluxMp);

        }

        // Solve average size moment transport
        fvScalarMatrix mEqn
        (
            fvm::ddt(m)
          - fvc::ddt(m)
          + meanDivUbMp
        );
        mEqn.relax();
        mEqn.solve();
    }

    // update now the velocity moments
    // reset all velocity moments to the mean velocity
    //forAll(quadrature_().velocityMoments(), mi)
    //{
    //    quadrature_().velocityMoments()[mi] = U_()*quadrature_().moments()[mi];
    //    quadrature_().velocityMoments()[mi].correctBoundaryConditions();
    //}

    //quadrature_().updateQuadrature();

    quadrature_().updateAllQuadrature();
	if(areThereNodesSmallerThanZero()){
		correctNodesAndMoments();
	}

    quadrature_().updateAllMoments();

    //- Update mean velocity
    //this->calcVelocity(); 

}
// realizableCo() wird nicht mehr verwendet.
/*
scalar Foam::diameterModels::populationBalanceMomModel::realizableCo()
{
    surfaceVectorField Sf(mesh_.Sf());
    //scalarField maxCoNum(mesh_.nCells(), 1.0);
    volScalarField maxCoNum
    (
//	new volScalarField
//	(
	    IOobject
	    (
		"maxCoNum",
	        mesh_.time().timeName(),
	        mesh_,
	        IOobject::NO_READ,
	        IOobject::NO_WRITE
	    ),
	    mesh_,
	    dimensionedScalar("maxCoNum", dimless, 1.0)
//	)
    );
    word weightInterpolation = mesh_.interpolationScheme("reconstruct(weight)");

    forAll(quadrature_().nodes(), nodei)
    {
    
        surfaceScalarField phiOwn
        (
            mag(quadrature_().velocitiesOwn()[nodei] & Sf)
        );
    
        surfaceScalarField phiNei
        (
            mag(quadrature_().velocitiesNei()[nodei] & Sf)
        );
    
	if(weightInterpolation == "upwind")
	{ 
	    forAll(quadrature_().moments()[0], celli)
	    {
	        const labelList& cell = mesh_.cells()[celli];
	
	        scalar den = 0;
	
	        forAll(cell, facei)
	        {
	            if (cell[facei] < mesh_.nInternalFaces())
	            {
	                den +=
	                    max
	                    (
	                        phiOwn[cell[facei]],
	                        phiNei[cell[facei]]
	                    );
	            }
	
	            den = max(den, small);
	
	            maxCoNum[celli] =
	                min
	                (
	                    maxCoNum[celli],
	                    mesh_.V()[celli]
	                   /(den*mesh_.time().deltaTValue())
	                );
	        }
	    }
	}

	else
	{
	    forAll(quadrature_().moments()[0], celli)
	    {
	        const labelList& cell = mesh_.cells()[celli];
	    
	        scalar num = quadrature_().nodes()[nodei].primaryWeight()[celli];
	        scalar den = 0;
	        forAll(cell, facei)
	        {
	            if (cell[facei] < mesh_.nInternalFaces())
	            {
	                if (mesh_.owner()[cell[facei]] == celli)
	                {
	                    den +=
	                        quadrature_().nodesOwn()[nodei].primaryWeight()[cell[facei]]
	                       *max(phiOwn[cell[facei]], 0.0);
	                }
	                else if (mesh_.neighbour()[cell[facei]] == celli)
	                {
	                    den -=
	                        quadrature_().nodesNei()[nodei].primaryWeight()[cell[facei]]
	                       *min(phiNei[cell[facei]], 0.0);
	                }
	            }
	            if (num > 1e-6)
	            {
	                den = max(den, small);
	                maxCoNum[celli] =
	                    min
	                    (
	                        maxCoNum[celli],
	                        num*mesh_.V()[celli]
	                       /(den*mesh_.time().deltaTValue())
	                    );
	            }
	        }
	    }
	}
    }

    Info<<"realizableCo "<<name_<<": "<<gMin(maxCoNum)<<"\n";
    if(mesh_.time().outputTime())
    {
	maxCoNum.write();
    }
    return gMin(maxCoNum);
}
*/

void Foam::diameterModels::populationBalanceMomModel::resetSTs()
{
    forAll(Sm_,k)
    {
    	Sm_[k] *= 0.;
    }
}

bool Foam::diameterModels::populationBalanceMomModel::areThereNodesSmallerThanZero()
{
    //check
    forAll(quadrature_().nodes(), nodei)
    {
        const volScalarNode& node = quadrature_().nodes()[nodei];
	dimensionedScalar minQuad = min(node.primaryAbscissae()[0]);
	if(minQuad.value() < 0.0)
	{
	     return true;
	}
    }
    return false;
}


void Foam::diameterModels::populationBalanceMomModel::correctNodesAndMoments()
{
    // correct nodes
    forAll(quadrature_().nodes(), nodei)
    {
        volScalarNode& node = quadrature_().nodes()[nodei];
	const dimensionedScalar zeroNodeVal("zeroValNode", node.primaryAbscissae()[0].dimensions(), 0.);
	bound( node.primaryAbscissae()[0], zeroNodeVal );
    }

    // correct moments using the nodes
    const PtrList<volScalarNode>& nodes = quadrature_().nodes();
    forAll(quadrature_().moments(), k)
    {
    	quadrature_().moments()[k] *= 0.;
	forAll(nodes, nodei)
	{
	    const volScalarNode& node = nodes[nodei];
	    const volScalarField abscissa(node.primaryAbscissae()[0]);
	    const volScalarField n(node.n(node.primaryWeight(), abscissa));
	    quadrature_().moments()[k] += n * pow(abscissa, k);
	}
    }
}


void Foam::diameterModels::populationBalanceMomModel::solve()
{
    Info<<"solve in populationBalanceMomModel"<<endl;
    quadrature_().updateAllQuadrature();

    if(areThereNodesSmallerThanZero())
    {
	correctNodesAndMoments();
    }
    //correct();
    calcAlphas();
    //calcVelocity();
    resetSTs();

    const dictionary& solutionControls = mesh_.solverDict(name_);
    bool solveOnFinalIterOnly = solutionControls.lookupOrDefault<bool>("solveOnFinalIterOnly", false);
    const label nCorr = mesh_.solverDict(name_).lookup<label>("nCorr");

    if (!solveOnFinalIterOnly || pimple_.finalPimpleIter())
    {
	Info << "solve now qmom Eqs" << endl;
	if(nCorr > 0)
	{
	    precompute();
	}
	
	// transport
	averageTransport();
	
	evalCoalescenceST();
	evalBreakupST();
	// update the moments by the break-up and coalescence source terms
	volScalarMomentFieldSet& moments = quadrature_().moments();
	dimensionedScalar globalDt("dtToSolve",moments[0].mesh().time().deltaT());
	forAll(moments, k)
	{
	    moments[k] += globalDt * Sm_[k];
	}
	// keep consistency to moments and alpha
	scaleMomentsToAlpha();
	// compute Co for next time step
	//realizableCo();
	    
	Info << "solving qmom Eqs finished" << endl;

    }
}

void Foam::diameterModels::populationBalanceMomModel::correct()
{
}

// ************************************************************************* //
