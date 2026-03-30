/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

Application
    pimpleFoam_calcResLES

Description
    Transient solver for incompressible, turbulent flow of Newtonian fluids,
    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.
    
    Extended with resolved field calculations and time averaging.
    Hence recommended for scale resolving methods.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "kinematicMomentumTransportModel.H"
#include "pimpleControl.H"
#include "pressureReference.H"
#include "CorrectPhi.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class TimeAverager
{
private:
    scalar startTime_;
    scalar actualAveragingTime_;
    bool initialized_;
    const scalar SMALL_TIME = 1e-14;

public:
    // Default constructor
    TimeAverager()
    :
        startTime_(0.0),
        actualAveragingTime_(0.0),
        initialized_(false)
    {}

    // Constructor with start time
    TimeAverager(scalar startTime, scalar currentTime)
    :
        startTime_(startTime),
        actualAveragingTime_(max(0.0, currentTime - startTime)),
        initialized_(true)
    {
        Info<< "TimeAverager initialized with start time: " << startTime_ 
        << " s, accumulated time: " << actualAveragingTime_ << " s" << endl;
    }

    inline void update(scalar currentTime, scalar deltaT)
    {
        // Initialize on first call if not yet initialized (only needed for default constructor)
        if (!initialized_)
        {
            startTime_ = currentTime;
            initialized_ = true;
            Info<< "TimeAverager starting at t = " << startTime_ << " s (default)" << endl;
        }

        // Normal accumulation
        scalar timeInAveraging = max(0.0, currentTime - startTime_);
        actualAveragingTime_ = min(actualAveragingTime_ + deltaT, timeInAveraging);
    }

    // Separate method to handle user-specified start time
    void setStartTime(scalar newStartTime, scalar currentTime)
    {
        if (mag(newStartTime - startTime_) > SMALL_TIME)
        {
            startTime_ = newStartTime;
            initialized_ = true;
            // Recalculate accumulated time
            actualAveragingTime_ = max(0.0, min(actualAveragingTime_, currentTime - startTime_));
            Info<< "Averaging start time updated to: " << startTime_ << " s" << endl;
        }
    }

    inline scalar getActualAveragingTime() const { return actualAveragingTime_; }
    inline scalar getStartTime() const { return startTime_; }
    inline bool isAveraging() const { return actualAveragingTime_ > SMALL_TIME; }

    template<typename T>
    inline void updateAverage(T& avgField, const T& instField, scalar deltaT)
    {
        if (actualAveragingTime_ <= deltaT + SMALL_TIME)
        {
            avgField = instField;
        }
        else
        {
            scalar oldWeight = (actualAveragingTime_ - deltaT) / actualAveragingTime_;
            avgField = oldWeight * avgField + (1.0 - oldWeight) * instField;
        }
    }
};

scalar readAveragingStartTimeOnce(const Time& runTime)
{
    const dictionary& controlDict = runTime.controlDict();

    if (controlDict.found("averagingStartTime"))
    {
        scalar startTime = controlDict.lookup<scalar>("averagingStartTime");
        Info<< "Averaging will start at t = " << startTime << " s (user-specified)" << endl;
        return startTime;
    }
    else
    {
        // Default to current simulation time
        scalar startTime = runTime.value();
        Info<< "No averagingStartTime specified - averaging will start at t = "
            << startTime << " s" << endl;
        return startTime;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createResolvedFields.H"
    #include "createUfIfPresent.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // LES constants
    //const scalar C_k = 0.094;

    // Calculate LES delta field once (mesh-dependent, doesn't change during simulation)
    volScalarField delta_LES
    (
        IOobject
        (
            "delta_LES",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("delta", dimLength, 0.0)
    );

    // Create delta field from cube root of cell volumes
    forAll(delta_LES, celli)
    {
        delta_LES[celli] = Foam::cbrt(mesh.V()[celli]);
    }

    // Add safety bounds for delta to prevent division by zero
    volScalarField delta_safe = max
    (
        delta_LES,
        dimensionedScalar("minDelta", dimLength, 1e-12)
    );

    // Initialize time averager with start time from controlDict
    scalar averagingStartTime = readAveragingStartTimeOnce(runTime);
    TimeAverager averager(averagingStartTime, runTime.value());

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (pimple.run(runTime))
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        const scalar currentDeltaT = runTime.deltaTValue();
        const scalar currentTime = runTime.value();

        // Always update the averager
        averager.update(currentTime, currentDeltaT);

        // Only check for start time changes occasionally (expensive file I/O)
        if (runTime.timeIndex() % 10 == 0)
        {
            const dictionary& controlDict = runTime.controlDict();
            if (controlDict.found("averagingStartTime"))
            {
                scalar newStartTime = controlDict.lookup<scalar>("averagingStartTime");
                if (mag(newStartTime - averagingStartTime) > 1e-10)
                {
                    averager.setStartTime(newStartTime, currentTime);
                    averagingStartTime = newStartTime;
                    Info<< "Averaging start time updated to: " << newStartTime << endl;
                }
            }
        }

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
            {
                fvModels.preUpdateMesh();

                mesh.update();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        #include "correctPhi.H"
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            fvModels.correct();

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        // ========== SAS RELEVANT FIELDS AND TIME AVERAGING =========
        // Other terms such as L, L_vk, QSAS need to be averaged via function object

        tmp<volTensorField> tgradU = fvc::grad(U);
        volScalarField S2_(2*magSqr(symm(tgradU())));
        volScalarField W2_(2*magSqr(skew(tgradU())));
        volScalarField S2byW2_(S2_/max(W2_, dimensionedScalar("small", W2_.dimensions(), 1e-9)));

        averager.updateAverage(S2Mean_, S2_, currentDeltaT);
        averager.updateAverage(W2Mean_, W2_, currentDeltaT);
        averager.updateAverage(S2byW2Mean_, S2byW2_, currentDeltaT);

        // ========== RESOLVED FIELD CALCULATIONS AND TIME AVERAGING ==========

        // Update time-averaged pressure
        averager.updateAverage(p_Mean_, p, currentDeltaT);

        // Update time-averaged velocity
        averager.updateAverage(UMean_, U, currentDeltaT);

        // Calculate fluctuating velocity (U' = U - <U>)
        UPrime_ = U - UMean_;

        // Get kinematic viscosity
        volScalarField nu = laminarTransport.nu();

        // ---- RESOLVED TURBULENT KINETIC ENERGY AND DISSIPATION ----
        volTensorField gradUPrime = fvc::grad(UPrime_);
        volSymmTensorField Sres(symm(gradUPrime));

        kRes_ = 0.5*Foam::magSqr(UPrime_);
        epsilonRes_ = 2.0*nu*Foam::magSqr(Sres);

        averager.updateAverage(kRes_Mean_, kRes_, currentDeltaT);
        averager.updateAverage(epsilonRes_Mean_, epsilonRes_, currentDeltaT);

        // ---- MODELLED/SGS FIELDS ----
        // Access turbulence fields directly from turbulence model
        kMod_ = turbulence->k();
        nuMod_ = turbulence->nut();
        epsilonMod_ = turbulence->epsilon();

        // Time-averaged modelled fields
        averager.updateAverage(kMod_Mean_, kMod_, currentDeltaT);
        averager.updateAverage(epsilonMod_Mean_, epsilonMod_, currentDeltaT);
        averager.updateAverage(nut_Mean_, nuMod_, currentDeltaT);

        // ---- TOTAL TURBULENT KINETIC ENERGY AND DISSIPATION ----
        kTotal_ = kRes_ + kMod_;
        epsilonTotal_ = epsilonRes_ + epsilonMod_;

        // Time-averaged totals
        averager.updateAverage(kTotal_Mean_, kTotal_, currentDeltaT);
        averager.updateAverage(epsilonTotal_Mean_, epsilonTotal_, currentDeltaT);

        // ---- LES GRID RESOLUTION ASSESSMENT ----
        volScalarField epsTotMean_safe = max
        (
            epsilonTotal_Mean_,
            dimensionedScalar("minEps", dimVelocity*dimVelocity/dimTime, 1e-12)
        );
        volScalarField kolmogorovScale = pow(pow(nu, 3)/epsTotMean_safe, 0.25);
        deltaByKolmog_ = delta_safe / kolmogorovScale;
        averager.updateAverage(deltaByKolmog_Mean_, deltaByKolmog_, currentDeltaT);

        // ---- REYNOLDS STRESS COMPONENTS ----
        uu_ = Foam::sqr(UPrime_.component(vector::X));
        vv_ = Foam::sqr(UPrime_.component(vector::Y));
        ww_ = Foam::sqr(UPrime_.component(vector::Z));
        uv_ = UPrime_.component(vector::X)*UPrime_.component(vector::Y);
        uw_ = UPrime_.component(vector::X)*UPrime_.component(vector::Z);
        vw_ = UPrime_.component(vector::Y)*UPrime_.component(vector::Z);

        // Time-averaged stresses
        averager.updateAverage(uu_Mean_, uu_, currentDeltaT);
        averager.updateAverage(vv_Mean_, vv_, currentDeltaT);
        averager.updateAverage(ww_Mean_, ww_, currentDeltaT);
        averager.updateAverage(uv_Mean_, uv_, currentDeltaT);
        averager.updateAverage(uw_Mean_, uw_, currentDeltaT);
        averager.updateAverage(vw_Mean_, vw_, currentDeltaT);

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
