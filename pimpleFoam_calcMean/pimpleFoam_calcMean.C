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

        // ========== RESOLVED FIELD CALCULATIONS AND TIME AVERAGING ==========

        // Update time-averaged pressure
        averager.updateAverage(p_Mean_, p, currentDeltaT);

        // Update time-averaged velocity
        averager.updateAverage(UMean_, U, currentDeltaT);

        // Calculate fluctuating velocity (U' = U - <U>)
        UPrime_ = U - UMean_;

        // Get kinematic viscosity
        volScalarField nu = laminarTransport.nu();

        // ---- MODELLED/SGS FIELDS ----
        // Access turbulence fields directly from turbulence model
        kMod_ = turbulence->k();
        nuMod_ = turbulence->nut();
        epsilonMod_ = turbulence->epsilon();

        // Time-averaged modelled fields
        averager.updateAverage(kMod_Mean_, kMod_, currentDeltaT);
        averager.updateAverage(epsilonMod_Mean_, epsilonMod_, currentDeltaT);
        averager.updateAverage(nut_Mean_, nuMod_, currentDeltaT);

        // ---- UPrime2Mean ----
        UPrime2_ = Foam::sqr(UPrime_);
        averager.updateAverage(UPrime2_Mean_, UPrime2_, currentDeltaT);

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
