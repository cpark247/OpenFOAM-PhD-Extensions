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
    multiphaseEulerFoamHSM_RANS

Description
    Solver for a system of any number of compressible fluid phases with a
    common pressure, but otherwise separate properties. The type of phase model
    is run time selectable and can optionally represent multiple species and
    in-phase reactions. The phase system is also run time selectable and can
    optionally represent different types of momentum, heat and mass transfer.

    Extended with RANS field calculations and time averaging for multiphase flows.
    Includes population balance rate averaging (coalescence and breakup per size class).
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "phaseSystem.H"
#include "phaseDynamicMomentumTransportModel.H"
#include "pimpleControl.H"
#include "pressureReference.H"
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

void readBlendingParameters(const Time& runTime, scalar& alphaMax, scalar& alphaMin)
{
    // Read phase properties
    IOdictionary phaseProperties
    (
        IOobject
        (
            "phaseProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Default values
    alphaMax = 0.8;
    alphaMin = 1e-4;

    // Try to read from blending subdictionary
    if (phaseProperties.found("blending"))
    {
        const dictionary& blending = phaseProperties.subDict("blending");

        if (blending.found("drag"))
        {
            const dictionary& drag = blending.subDict("drag");

            if (drag.found("alphaMax"))
            {
                alphaMax = drag.lookup<scalar>("alphaMax");
                Info<< "Read alphaMax = " << alphaMax << " from phaseProperties" << endl;
            }

            if (drag.found("alphaMin"))
            {
                alphaMin = drag.lookup<scalar>("alphaMin");
                Info<< "Read alphaMin = " << alphaMin << " from phaseProperties" << endl;
            }
        }
    }

    if (alphaMax == 0.8 && alphaMin == 1e-4)
    {
        Info<< "Using default blending parameters: alphaMax = " << alphaMax
            << ", alphaMin = " << alphaMin << endl;
    }
}

int main(int argc, char *argv[])
{
    Info<< "hsm RANS solver with time averaging" << endl;

    #include "postProcess.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createMultiphaseFieldsRANS.H"
    #include "createPopBalRateFields.H"              // NEW: Population balance rate fields
    #include "createFieldRefs.H"

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // Switch declarations - OUTSIDE the time loop
    Switch faceMomentum
    (
        pimple.dict().lookupOrDefault<Switch>("faceMomentum", false)
    );
    Switch partialElimination
    (
        pimple.dict().lookupOrDefault<Switch>("partialElimination", false)
    );

    #include "createRDeltaTf.H"

    // Read blending parameters for CBeta calculation
    scalar alphaMax, alphaMin;
    readBlendingParameters(runTime, alphaMax, alphaMin);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Initialize time averager with start time from controlDict
    scalar averagingStartTime = readAveragingStartTimeOnce(runTime);
    TimeAverager averager(averagingStartTime, runTime.value());

    Info<< nl << "Starting time loop" << nl << endl;

    while (pimple.run(runTime))
    {
        #include "readDyMControls.H"

        // if (LTS) { #include "setRDeltaT.H"; if (faceMomentum) { #include "setRDeltaTf.H" } }
        // else
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl;
        Info<< "deltaT = " << runTime.deltaTValue() << nl << endl;

        const scalar currentDeltaT = runTime.deltaTValue();
        const scalar currentTime = runTime.value();

        // Always update the averager (fast, minimal branching)
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
                // Store divU from the previous mesh so that it can be
                // mapped and used in correctPhi to ensure the corrected phi
                // has the same divergence
                tmp<volScalarField> divU;

                if (correctPhi)
                {
                    // Construct and register divU for mapping
                    divU = new volScalarField
                    (
                        "divU0",
                        fvc::div
                        (
                            fvc::absolute
                            (
                                phi,
                                fluid.movingPhases()[0].U()
                            )
                        )
                    );
                }

                fvModels.preUpdateMesh();
                mesh.update();
            }

            if (pimple.models())
            {
                fvModels.correct();
            }

            fluid.solve(rAUs, rAUfs);
            fluid.correct();
            fluid.correctContinuityError();

            #include "pU/UEqns.H"

            /******************* pEqn inclusive Matrix Relaxation ***************/
            #include "pU/pEqnRelax.H"

            fluid.correctKinematics();

            if (pimple.turbCorr())
            {
                fluid.correctTurbulence();
            }
        } // end pimple.loop()

        // Update time-averaged pressure
        averager.updateAverage(p_rgh_Mean_, p_rgh, currentDeltaT);
        averager.updateAverage(p_Mean_, p, currentDeltaT);

        // ========== MULTIPHASE FIELD CALCULATIONS ==========

        // Loop through all phases for averaging
        forAll(phases, phasei)
        {
            const phaseModel& phase = phases[phasei];
            const word phaseName = phase.name();

            // Get phase velocity
            const volVectorField& U_phase = phase.U();

            // Update time-averaged phase velocity
            averager.updateAverage(UMean_[phasei], U_phase, currentDeltaT);

            // Try to access turbulence fields for this phase
            const word kFieldName     = "k."     + phaseName;
            const word nutFieldName   = "nut."   + phaseName;
            const word omegaFieldName = "omega." + phaseName;

            if (phases[phasei].name() == "air")
            {
                // Get phase fraction as volScalarField
                const volScalarField& alpha_phase = phases[phasei];
                averager.updateAverage(alphaMean_[phasei], alpha_phase, currentDeltaT);

                // Calculate isovolume fields (1 where alpha.air > threshold, 0 otherwise)
                isoVol_80 = pos(alpha_phase - 0.8);
                isoVol_90 = pos(alpha_phase - 0.9);
                isoVol_95 = pos(alpha_phase - 0.95);

                // Calculate CBeta field
                volScalarField alpha_phase_safe =
                    max(alpha_phase, dimensionedScalar("small", dimless, 1e-12));
                volScalarField one_minus_alpha_phase = 1.0 - alpha_phase;

                volScalarField term1 =
                    (one_minus_alpha_phase*alphaMax)
                  / ((1.0 - alphaMax)*alpha_phase_safe);

                volScalarField term2 = alphaMin/alpha_phase_safe;

                CBeta = min
                (
                    max(term1, term2),
                    dimensionedScalar("one", dimless, 1.0)
                );

                // Update time-averaged CBeta
                averager.updateAverage(CBeta_Mean_, CBeta, currentDeltaT);
            }
            else // for water
            {
                // Get modelled turbulent fields from object registry
                kMod_[phasei]     = mesh.lookupObject<volScalarField>(kFieldName);
                omegaMod_[phasei] = mesh.lookupObject<volScalarField>(omegaFieldName);
                nuMod_[phasei]    = mesh.lookupObject<volScalarField>(nutFieldName);

                // Update time-averaged modelled fields
                averager.updateAverage(kMod_Mean_[phasei],     kMod_[phasei],     currentDeltaT);
                averager.updateAverage(omegaMod_Mean_[phasei], omegaMod_[phasei], currentDeltaT);
                averager.updateAverage(nut_Mean_[phasei],      nuMod_[phasei],    currentDeltaT);
            }
        } // end phases loop

        // ========== POPULATION BALANCE RATE AVERAGING ========== (NEW)
        #include "updatePopBalRateAverages.H"

        // ========== D32 CALCULATION ========== (UPDATED)
        #include "calculateD32.H"

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime()  << " s"
            << "  ClockTime = "   << runTime.elapsedClockTime() << " s"
            << nl << endl;

    } // end time loop

    Info<< "End\n" << endl;
    return 0;
}
// ************************************************************************* //
