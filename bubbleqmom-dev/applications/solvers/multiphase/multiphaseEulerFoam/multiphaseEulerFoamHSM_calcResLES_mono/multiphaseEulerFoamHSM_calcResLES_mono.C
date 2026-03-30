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
    multiphaseEulerFoamHSM

Description
    Solver for a system of any number of compressible fluid phases with a
    common pressure, but otherwise separate properties. The type of phase model
    is run time selectable and can optionally represent multiple species and
    in-phase reactions. The phase system is also run time selectable and can
    optionally represent different types of momentum, heat and mass transfer.

    Extended with resolved field calculations and time averaging for multiphase flows.
    
    Supports both LES models and RANS SRS models (kOmegaSSTSAS):
    - LES mode: kMod computed from nut and delta
    - RANS SRS mode: kMod directly from turbulence model's k field,
                     epsilonMod from k and omega
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
    Info<< "hsm solver with time averaging (LES + RANS SRS support)" << endl;

    #include "postProcess.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createMultiphaseResolvedLESFields.H"
    #include "createFieldRefs.H"

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

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

    // LES constants (only used if not in RANS SRS mode)
    const scalar C_k = 0.094;
    
    // k-omega model constant for epsilon calculation: epsilon = beta_star * k * omega
    const scalar beta_star = 0.09;

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

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Initialize time averager with start time from controlDict - re-read every timestep for run-time update
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
            if (controlDict.found("averagingStartTime"))  // Only check if it exists
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
            // if (!pimple.flow()) { ... } else {
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
                            fvc::absolute(phi, fluid.movingPhases()[0].U())
                        )
                    );
                }

                fvModels.preUpdateMesh();
                mesh.update();

                // if (mesh.changing()) { ... } // (kept commented)
            }

            if (pimple.models())
            {
                fvModels.correct();
            }

            fluid.solve(rAUs, rAUfs);
            fluid.correct();
            fluid.correctContinuityError();

            // if (pimple.thermophysics()) { #include "YEqns.H" }

            // if (faceMomentum) { #include "pUf/UEqns.H"; ... #include "pUf/pEqn.H" }
            // else
            #include "pU/UEqns.H"

            // if (pimple.thermophysics()) { #include "EEqns.H" }

            /******************* pEqn inclusive Matrix Relaxation ***************/
            #include "pU/pEqnRelax.H"
            // } // end faceMomentum branch

            fluid.correctKinematics();

            if (pimple.turbCorr())
            {
                fluid.correctTurbulence();
            }
        } // end pimple.loop()

        // Update time-averaged pressure
        averager.updateAverage(p_rgh_Mean_, p_rgh, currentDeltaT);
        averager.updateAverage(p_Mean_,     p,     currentDeltaT);

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

            // Calculate phase fluctuating velocity (ASSIGN, do not declare)
            UPrime_[phasei] = U_phase - UMean_[phasei];

            // Field names for turbulence variables
            const word kFieldName       = "k."       + phaseName;
            const word omegaFieldName   = "omega."   + phaseName;
            const word nutFieldName     = "nut."     + phaseName;
            const word epsilonFieldName = "epsilon." + phaseName;

            if (phases[phasei].name() == "air")
            {
                // Get phase fraction as volScalarField
                const volScalarField& alpha_phase = phases[phasei];
                averager.updateAverage(alphaMean_[phasei], alpha_phase, currentDeltaT);

                // Isovolume fields
                isoVol_80 = pos(alpha_phase - 0.8);
                isoVol_90 = pos(alpha_phase - 0.9);
                isoVol_95 = pos(alpha_phase - 0.95);

            }
            else // water phase
            {
                // Common calculations for both LES and RANS SRS modes
                volScalarField nu_phase = phase.thermo().nu();
                volTensorField gradUPrime_phase = fvc::grad(UPrime_[phasei]);
                volSymmTensorField Sres_phase(symm(gradUPrime_phase));
                
                // Resolved turbulent kinetic energy: k_res = 0.5 * |U'|^2
                kRes_[phasei] = 0.5*Foam::magSqr(UPrime_[phasei]);
                
                // Resolved dissipation rate: epsilon_res = 2 * nu * |S'|^2
                epsilonRes_[phasei] = 2.0*nu_phase*Foam::magSqr(Sres_phase);

                averager.updateAverage(kRes_Mean_[phasei],       kRes_[phasei],       currentDeltaT);
                averager.updateAverage(epsilonRes_Mean_[phasei], epsilonRes_[phasei], currentDeltaT);

                // Get nut (turbulent viscosity) - available for both LES and RANS
                nuMod_[phasei] = mesh.lookupObject<volScalarField>(nutFieldName);

                // ============================================================
                // MODELLED k AND epsilon CALCULATION
                // Branch based on turbulence model type
                // ============================================================
                if (useRANSSRS_)
                {
                    // ---------------------------------------------------------
                    // RANS SRS MODE (kOmegaSSTSAS)
                    // k and omega fields exist and are directly solved
                    // ---------------------------------------------------------
                    
                    // Get k directly from turbulence model
                    const volScalarField& k_turb = 
                        mesh.lookupObject<volScalarField>(kFieldName);
                    kMod_[phasei] = k_turb;
                    
                    // Get omega directly from turbulence model
                    const volScalarField& omega_turb = 
                        mesh.lookupObject<volScalarField>(omegaFieldName);
                    omegaMod_[phasei] = omega_turb;
                    
                    // Calculate epsilon from k and omega using:
                    // epsilon = beta_star * k * omega (where beta_star = 0.09)
                    epsilonMod_[phasei] = beta_star * k_turb * omega_turb;
                    
                    // Time-averaged omega (only for RANS SRS)
                    averager.updateAverage(omegaMod_Mean_[phasei], omegaMod_[phasei], currentDeltaT);
                }
                else
                {
                    // ---------------------------------------------------------
                    // LES MODE
                    // Compute kMod from nut and delta (Smagorinsky-like relation)
                    // ---------------------------------------------------------
                    
                    // kMod = (nut / (C_k * delta))^2
                    volScalarField denominator_k = C_k * delta_safe;
                    kMod_[phasei] = sqr(nuMod_[phasei] / denominator_k);
                    
                    // epsilonMod = 2 * nut * |S'|^2  (SGS dissipation estimate)
                    epsilonMod_[phasei] = 2.0 * nuMod_[phasei] * Foam::magSqr(Sres_phase);
                    
                    // Alternative: epsilonMod from kMod (commented out)
                    // epsilonMod_[phasei] = C_eps * pow(max(kMod_[phasei], 
                    //     dimensionedScalar("minK", dimVelocity*dimVelocity, SMALL_K)), 1.5) / delta_safe;
                }

                // Time-averaged modelled fields (common for both modes)
                averager.updateAverage(kMod_Mean_[phasei],       kMod_[phasei],       currentDeltaT);
                averager.updateAverage(epsilonMod_Mean_[phasei], epsilonMod_[phasei], currentDeltaT);
                averager.updateAverage(nut_Mean_[phasei],        nuMod_[phasei],      currentDeltaT);

                // Total fields (resolved + modelled)
                kTotal_[phasei]       = kRes_[phasei] + kMod_[phasei];
                epsilonTotal_[phasei] = epsilonRes_[phasei] + epsilonMod_[phasei];

                // Time-averaged totals
                averager.updateAverage(kTotal_Mean_[phasei],       kTotal_[phasei],       currentDeltaT);
                averager.updateAverage(epsilonTotal_Mean_[phasei], epsilonTotal_[phasei], currentDeltaT);

                // LES/SRS grid resolution assessment using Kolmogorov scale
                volScalarField epsTotMean_safe = max(
                    epsilonTotal_Mean_[phasei], 
                    dimensionedScalar("minEps", dimVelocity*dimVelocity/dimTime, 1e-12)
                );
                volScalarField kolmogorovScale = pow(pow(nu_phase, 3)/epsTotMean_safe, 0.25);
                deltaByKolmog = delta_safe / kolmogorovScale;
                averager.updateAverage(deltaByKolmog_Mean_, deltaByKolmog, currentDeltaT);
      
                //Info<<"deltaByKolmog: max = " << max(deltaByKolmog).value() 
                //    << ", average = " << average(deltaByKolmog).value() << endl;

                // Phase Reynolds stress components
                uu_[phasei] = Foam::sqr(UPrime_[phasei].component(vector::X));
                vv_[phasei] = Foam::sqr(UPrime_[phasei].component(vector::Y));
                ww_[phasei] = Foam::sqr(UPrime_[phasei].component(vector::Z));
                uv_[phasei] = UPrime_[phasei].component(vector::X)*UPrime_[phasei].component(vector::Y);
                uw_[phasei] = UPrime_[phasei].component(vector::X)*UPrime_[phasei].component(vector::Z);
                vw_[phasei] = UPrime_[phasei].component(vector::Y)*UPrime_[phasei].component(vector::Z);

                // Time-averaged stresses
                averager.updateAverage(uu_Mean_[phasei], uu_[phasei], currentDeltaT);
                averager.updateAverage(vv_Mean_[phasei], vv_[phasei], currentDeltaT);
                averager.updateAverage(ww_Mean_[phasei], ww_[phasei], currentDeltaT);
                averager.updateAverage(uv_Mean_[phasei], uv_[phasei], currentDeltaT);
                averager.updateAverage(uw_Mean_[phasei], uw_[phasei], currentDeltaT);
                averager.updateAverage(vw_Mean_[phasei], vw_[phasei], currentDeltaT);
            }
        } // end phases loop

        
    
        runTime.write();
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime()
            << " s  ClockTime = "   << runTime.elapsedClockTime()
            << " s" << nl << endl;

    } // end time loop

    Info<< "End\n" << endl;
    return 0;
}
// ************************************************************************* //
