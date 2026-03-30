/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
\*---------------------------------------------------------------------------*/

Info<< "Creating multiphase resolved turbulence fields\n" << endl;

// Get number of phases
const label nPhases = phases.size();

// ============================================================================
// TURBULENCE MODEL TYPE DETECTION
// ============================================================================
// Check if we're using a RANS SRS model (kOmegaSSTSAS) by checking for k and
// omega fields. For LES models, these fields don't exist and we compute kMod
// from nut and delta.

bool useRANSSRS_ = false;
word continuousPhaseName = "water";  // Adjust if your continuous phase has a different name

// Check for k and omega fields for the continuous phase
const word kFieldName = "k." + continuousPhaseName;
const word omegaFieldName = "omega." + continuousPhaseName;

if (mesh.foundObject<volScalarField>(kFieldName) && 
    mesh.foundObject<volScalarField>(omegaFieldName))
{
    useRANSSRS_ = true;
    Info<< "Detected RANS SRS model (kOmegaSSTSAS): k and omega fields found" << nl
        << "    -> kMod will use turbulence model's k field directly" << nl
        << "    -> epsilonMod will be computed from k and omega" << endl;
}
else
{
    Info<< "Using LES model: k and omega fields not found" << nl
        << "    -> kMod will be computed from nut and delta" << nl
        << "    -> epsilonMod will be computed from resolved strain rate" << endl;
}

// ============================================================================
// TIME-AVERAGED PRESSURE
// ============================================================================
volScalarField p_rgh_Mean_
(
    IOobject
    (
        "p_rgh_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    p_rgh
);

volScalarField p_Mean_
(
    IOobject
    (
        "p_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    p
);

// ============================================================================
// PHASE VELOCITIES (MEAN) + FLUCTUATIONS (TRANSIENT ONLY)
// ============================================================================
PtrList<volVectorField> UMean_(nPhases);

// Container for fluctuating velocity (U') — transient only, not written
PtrList<volVectorField> UPrime_(phases.size());

// ============================================================================
// MIXTURE PROPERTIES (MEAN)
// ============================================================================
PtrList<volScalarField> alphaMean_(nPhases);

// ============================================================================
// MODELLED/SGS FIELDS (INSTANTANEOUS) PER PHASE
// ============================================================================
PtrList<volScalarField> kMod_(nPhases);
PtrList<volScalarField> nuMod_(nPhases);
PtrList<volScalarField> epsilonMod_(nPhases);
PtrList<volScalarField> omegaMod_(nPhases);  // For RANS SRS models

// ============================================================================
// MODELLED/SGS FIELDS (MEAN)
// ============================================================================
PtrList<volScalarField> kMod_Mean_(nPhases);
PtrList<volScalarField> epsilonMod_Mean_(nPhases);
PtrList<volScalarField> nut_Mean_(nPhases);
PtrList<volScalarField> omegaMod_Mean_(nPhases);  // For RANS SRS models

// ============================================================================
// RESOLVED FIELDS (INSTANTANEOUS) PER PHASE
// ============================================================================
PtrList<volScalarField> kRes_(nPhases);
PtrList<volScalarField> epsilonRes_(nPhases);

// ============================================================================
// RESOLVED FIELDS (MEAN)
// ============================================================================
PtrList<volScalarField> kRes_Mean_(nPhases);
PtrList<volScalarField> epsilonRes_Mean_(nPhases);

// ============================================================================
// TOTAL FIELDS (INSTANTANEOUS) PER PHASE
// ============================================================================
PtrList<volScalarField> kTotal_(nPhases);
PtrList<volScalarField> epsilonTotal_(nPhases);

// ============================================================================
// TOTAL FIELDS (MEAN)
// ============================================================================
PtrList<volScalarField> kTotal_Mean_(nPhases);
PtrList<volScalarField> epsilonTotal_Mean_(nPhases);

// ============================================================================
// REYNOLDS STRESS COMPONENTS (INSTANTANEOUS) PER PHASE
// ============================================================================
PtrList<volScalarField> uu_(nPhases);
PtrList<volScalarField> vv_(nPhases);
PtrList<volScalarField> ww_(nPhases);
PtrList<volScalarField> uv_(nPhases);
PtrList<volScalarField> uw_(nPhases);
PtrList<volScalarField> vw_(nPhases);

// Time-averaged Reynolds stress components
PtrList<volScalarField> uu_Mean_(nPhases);
PtrList<volScalarField> vv_Mean_(nPhases);
PtrList<volScalarField> ww_Mean_(nPhases);
PtrList<volScalarField> uv_Mean_(nPhases);
PtrList<volScalarField> uw_Mean_(nPhases);
PtrList<volScalarField> vw_Mean_(nPhases);

// ============================================================================
// BUILD PER-PHASE FIELDS
// ============================================================================
forAll(phases, phasei)
{
    const word phaseName = phases[phasei].name();

    // ---- UMean_[phasei]
    UMean_.set
    (
        phasei,
        new volVectorField
        (
            IOobject
            (
                "UMean_." + phaseName,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            phases[phasei].U()
        )
    );

    // ---- UPrime_[phasei] (transient only)
    UPrime_.set
    (
        phasei,
        new volVectorField
        (
            IOobject
            (
                "UPrime_tmp",              // dummy name, not written
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE         // ensure not written to disk
            ),
            mesh,
            dimensionedVector("zero", dimVelocity, vector::zero)
        )
    );

    if (phases[phasei].name() == "water")
    {
        // ---------------- Modelled (instantaneous)
        kMod_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "kMod." + phaseName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("kMod", dimVelocity*dimVelocity, 0.0)
            )
        );

        nuMod_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "nuMod." + phaseName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("nuMod", dimViscosity, 0.0)
            )
        );

        epsilonMod_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "epsilonMod." + phaseName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    "epsilonMod",
                    dimVelocity*dimVelocity/dimTime,
                    0.0
                )
            )
        );

        // omega field for RANS SRS models
        omegaMod_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "omegaMod." + phaseName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("omegaMod", dimless/dimTime, 0.0)
            )
        );

        // ---------------- Modelled (mean)
        kMod_Mean_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "kMod_Mean_." + phaseName,
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                kMod_[phasei]
            )
        );

        epsilonMod_Mean_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "epsilonMod_Mean_." + phaseName,
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                epsilonMod_[phasei]
            )
        );

        nut_Mean_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "nut_Mean_." + phaseName,
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                nuMod_[phasei]
            )
        );

        // omega mean for RANS SRS models
        omegaMod_Mean_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "omegaMod_Mean_." + phaseName,
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    useRANSSRS_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
                ),
                omegaMod_[phasei]
            )
        );

        // ---------------- Resolved (instantaneous)
        kRes_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "kRes." + phaseName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("kRes", dimVelocity*dimVelocity, 0.0)
            )
        );

        epsilonRes_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "epsilonRes." + phaseName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    "epsilonRes",
                    dimVelocity*dimVelocity/dimTime,
                    0.0
                )
            )
        );

        // ---------------- Resolved (mean)
        kRes_Mean_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "kRes_Mean_." + phaseName,
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                kRes_[phasei]
            )
        );

        epsilonRes_Mean_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "epsilonRes_Mean_." + phaseName,
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                epsilonRes_[phasei]
            )
        );

        // ---------------- Total (instantaneous)
        kTotal_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "kTotal." + phaseName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("kTotal", dimVelocity*dimVelocity, 0.0)
            )
        );

        epsilonTotal_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "epsilonTotal." + phaseName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    "epsilonTotal",
                    dimVelocity*dimVelocity/dimTime,
                    0.0
                )
            )
        );

        // ---------------- Total (mean)
        kTotal_Mean_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "kTotal_Mean_." + phaseName,
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                kTotal_[phasei]
            )
        );

        epsilonTotal_Mean_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "epsilonTotal_Mean_." + phaseName,
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                epsilonTotal_[phasei]
            )
        );

        // ---------------- Reynolds stresses (instantaneous)
        uu_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "uu." + phaseName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("uu", dimVelocity*dimVelocity, 0.0)
            )
        );

        vv_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "vv." + phaseName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("vv", dimVelocity*dimVelocity, 0.0)
            )
        );

        ww_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "ww." + phaseName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("ww", dimVelocity*dimVelocity, 0.0)
            )
        );

        uv_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "uv." + phaseName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("uv", dimVelocity*dimVelocity, 0.0)
            )
        );

        uw_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "uw." + phaseName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("uw", dimVelocity*dimVelocity, 0.0)
            )
        );

        vw_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "vw." + phaseName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("vw", dimVelocity*dimVelocity, 0.0)
            )
        );

        // ---------------- Reynolds stresses (mean)
        uu_Mean_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "uu_Mean_." + phaseName,
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                uu_[phasei]
            )
        );

        vv_Mean_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "vv_Mean_." + phaseName,
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                vv_[phasei]
            )
        );

        ww_Mean_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "ww_Mean_." + phaseName,
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                ww_[phasei]
            )
        );

        uv_Mean_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "uv_Mean_." + phaseName,
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                uv_[phasei]
            )
        );

        uw_Mean_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "uw_Mean_." + phaseName,
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                uw_[phasei]
            )
        );

        vw_Mean_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "vw_Mean_." + phaseName,
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                vw_[phasei]
            )
        );
    }
    else // air phase
    {
        alphaMean_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "alpha." + phaseName + "_Mean",
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                phases[phasei]
            )
        );
    }
}

// ============================================================================
// POST-PROCESS FIELDS (AIR) — ISOVOLUMES, CBeta, d32, deltaByKolmog
// ============================================================================
volScalarField isoVol_80
(
    IOobject
    (
        "isoVol_air_0.8",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("isoVol_80", dimless, 0.0)
);

volScalarField isoVol_90
(
    IOobject
    (
        "isoVol_air_0.9",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("isoVol_90", dimless, 0.0)
);

volScalarField isoVol_95
(
    IOobject
    (
        "isoVol_air_0.95",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("isoVol_95", dimless, 0.0)
);

volScalarField deltaByKolmog
(
    IOobject
    (
        "deltaByKolmog",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("deltaByKolmog", dimless, 0.0)
);

volScalarField deltaByKolmog_Mean_
(
    IOobject
    (
        "deltaByKolmog_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    deltaByKolmog
);


Info<< "Multiphase resolved turbulence fields created for "
    << nPhases << " phases" << nl;

if (useRANSSRS_)
{
    Info<< "Using RANS SRS mode (kOmegaSSTSAS detected)" << endl;
}
else
{
    Info<< "Using LES mode" << endl;
}

Info<< endl;
