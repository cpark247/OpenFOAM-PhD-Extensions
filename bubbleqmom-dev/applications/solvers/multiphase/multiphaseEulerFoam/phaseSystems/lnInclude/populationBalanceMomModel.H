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

Class
    Foam::diameterModels::populationBalanceMomModel

Description
    Class that solves the univariate population balance equation by means of
    a qmom-method. The internal coordinate is set to the particle volume, so the equation is based on
    a transport equation of the volume-based number density function.
	The source terms are written in a way that particle
    number and mass are preserved. Coalescence (aggregation), breakup, drift
    (growth and surface loss) as well as nucleation are supported.
    For the discrete breakup term two recipes are available, depending on the
    model choice. For models which state a total breakup rate and a separate
    daughter size distribution function, the formulation of Kumar and Ramkrishna
    (1996) is applied which is applicable for binary and multiple breakup
    events. The second formulation is given by Liao et al. (2017). It is useful
    for binary breakup models which give the breakup rate between a sizeGroup
    pair directly, without an explicit expression for the daughter size
    distribution. The drift term is implemented using a finite difference upwind
    scheme. Although it is diffusive, it ensures a stable and
    number-conservative solution.


Usage
    Example excerpt from a phaseProperties dictionary.
    \verbatim
    type populationBalanceTwoPhaseSystem;

    phases (air water);

    populationBalances (bubbles);

    air
    {
        type            purePhaseModel;
        diameterModel   qmomVelocity;
        velocityGroupCoeffs
        {
            populationBalance    bubbles;

            shapeModel		spherical;

        }

        residualAlpha   1e-6;
    }

    populationBalanceCoeffs
    {
        bubbles
        {
            continuousPhase water;

            coalescenceModels
            (
                hydrodynamic
                {
                    C 0.25;
                }
            );

            binaryBreakupModels
            ();

            breakupModels
            (
                exponential
                {
                    C 0.5;
                    exponent 0.01;
                    daughterSizeDistributionModel uniformBinary;
                }
            );

            driftModels
            (
                densityChange{}
            );

            nucleationModels
            ();
        }
    }
    \endverbatim

See also
    Foam::diameterModels::sizeGroup
    Foam::diameterModels::velocityGroup

SourceFiles
    populationBalanceMomModel.C

\*---------------------------------------------------------------------------*/

#ifndef populationBalanceMomModel_H
#define populationBalanceMomModel_H

#include "qmomVelocity.H"
#include "phasePair.H"
#include "pimpleControl.H"
#include "phaseDynamicMomentumTransportModel.H"
#include "HashPtrTable.H"
#include "monoKineticQuadratureApproximation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phaseSystem;

namespace diameterModels
{

class coalescenceMomModel;
class breakupMomModel;
class binaryBreakupMomModel;
class driftMomModel;
class nucleationMomModel;

/*---------------------------------------------------------------------------*\
                   Class populationBalanceMomModel Declaration
\*---------------------------------------------------------------------------*/

class populationBalanceMomModel
:
    public regIOobject
{
    // Private typedefs

        typedef
            HashTable<autoPtr<phasePair>, phasePairKey, phasePairKey::hash>
            phasePairTable;

	typedef
	    HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash>
	    pDmdtTable;


    // Private data

        //- Reference to the phaseSystem
        const phaseSystem& fluid_;
		//
        //- Reference phase index
        label refPhaseIndex_;

        //- Interfacial Mass transfer rate
        pDmdtTable& pDmdt_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Name of the populationBalance
        word name_;

        //- Dictionary
        dictionary dict_;

        //- Reference to pimpleControl
        const pimpleControl& pimple_;

        //- Continuous phase
        const phaseModel& continuousPhase_;

	//- QuadratureApproximation
	//monoKineticQuadratureApproximation quadrature_;
        autoPtr<monoKineticQuadratureApproximation> quadrature_;

        //- velocityGroup belonging to this populationBalance
        autoPtr<qmomVelocity> qmomVel_;

        //- List of unordered phasePairs in this populationBalance
        phasePairTable phasePairs_;

        //- sources for size moments
        PtrList<volScalarField> Sm_;

        //- Coalescence models
        PtrList<coalescenceMomModel> coalescence_;

        //- Coalescence rate
        autoPtr<volScalarField> coalescenceRate_;

        //- BreakupModels
        PtrList<breakupMomModel> breakup_;

        //- Breakup rate
        autoPtr<volScalarField> breakupRate_;

        //- Breakup daughter distribution integral
        PtrList<volScalarField> breakDsdIntegr_;

        //- Total void fraction of phases belonging to this populationBalance
        autoPtr<volScalarField> alphas_;
		
        //- List of quadrature node alphas
        PtrList<volScalarField> nodeAlphas_;

        //- List of diameters
        PtrList<volScalarField> ds_;

        //- List of qmom mass nodes
        PtrList<volScalarField> QMOMNode_;

        //- List of qmom mass Weights
        PtrList<volScalarField> QMOMWeight_;

        //- mean Diameter (alpha weighted diameter for polydisperse)
        autoPtr<volScalarField> dsm_;

        //- Maximum diameter
        dimensionedScalar maxD_;

        //- Minimum diameter
        dimensionedScalar minD_;

        //- Return the residual phase-fraction for given phase
        //  Used to stabilize the phase momentum as the phase-fraction -> 0
        dimensionedScalar residualAlpha_;

        //- Average velocity
        //autoPtr<volVectorField> U_;

    // Private member functions

        void registerPhasesAndQmom();

	void precompute();

	//-  add coalescence source term on moment source term fields
        void evalCoalescenceST();

	//-  add coalescence source term on moment source term fields
        void evalBreakupST();

        void calcAlphas();

	//- checks wether a single node is less than zero
        bool areThereNodesSmallerThanZero();
	//- bound nodes, such they are greater than zero. Update then the moments
        void correctNodesAndMoments();

        //void calcVelocity();

	// scale moments so that m1/rho and alpha are consistent
        void scaleMomentsToAlpha();

public:

    friend class sizeGroup;
    friend class velocityGroup;

    // Constructor

        populationBalanceMomModel
        (
            const phaseSystem& fluid,
            const word& name,
            pDmdtTable& pDmdt
        );

        //- Return clone
        autoPtr<populationBalanceMomModel> clone() const;

        //- Return a pointer to a new populationBalanceMomModel object created on
        //  freestore from Istream
        class iNew
        {
            const phaseSystem& fluid_;

            pDmdtTable& pDmdt_;

        public:

            iNew
            (
                const phaseSystem& fluid,
                pDmdtTable& pDmdt
            )
            :
                fluid_(fluid),
                pDmdt_(pDmdt)
            {}

            autoPtr<populationBalanceMomModel> operator()(Istream& is) const
            {
                return autoPtr<populationBalanceMomModel>
                (
                    new populationBalanceMomModel(fluid_, word(is), pDmdt_)
                );
            }
        };


    //- Destructor
    virtual ~populationBalanceMomModel();

    // Member Functions

        //- Dummy write for regIOobject
        bool writeData(Ostream&) const;

        //- Return reference to the phaseSystem
        inline const phaseSystem& fluid() const;

        //- Return reference to the mesh
        inline const fvMesh& mesh() const;

        //- Return reference to the mesh
        const phaseModel& refPhase() const;

        //- Return populationBalanceCoeffs dictionary
        inline const dictionary& dict() const;

        //- Return continuous phase
        inline const phaseModel& continuousPhase() const;

        //- Return list of unordered phasePairs in this populationBalance
        inline const phasePairTable& phasePairs() const;

        //- Return total void of phases belonging to this populationBalance
        inline const volScalarField& alphas() const;

        //- Return average velocity
        //inline const volVectorField& U() const;

	//- returns the qudrature object
	inline monoKineticQuadratureApproximation& getQuadratureObj();

	//- returns the qudrature object
	inline const monoKineticQuadratureApproximation& getQuadratureObj() const;
		
        //- Return the number of nodes
        inline label nNodes() const;

        //- Return alpha field for nodei
        inline const volScalarField& nodeAlphas(const label nodei) const;

        //- Return non-constant access to alpha field for nodei
        inline volScalarField& nodeAlphas(const label nodei);

        //- Return the diameter for nodei
        inline const volScalarField& ds(const label nodei) const;

        //- Return the sauter mean diameter
        const volScalarField& dsm() const;

        //- Return the surface tension coefficient between a given dispersed
        //  and the continuous phase
        const tmp<volScalarField> sigmaWithContinuousPhase
        (
            const phaseModel& dispersedPhase
        ) const;

        //- Return reference to turbulence model of the continuous phase
        const phaseCompressible::momentumTransportModel& continuousTurbulence() const;

        //- Average transport of moments (i.e. convection with phase mean velocity)
        virtual void averageTransport();

        //- Reset the source term vectors 
        void resetSTs();

        //- Solve the population balance equation
        void solve();

	//- Correct derived quantities
        void correct();

//realizableCo() wird nicht mehr verwendet
	//- calculate the max Co number that satisfies the time step criteria by Vikas
//	scalar realizableCo();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "populationBalanceMomModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
