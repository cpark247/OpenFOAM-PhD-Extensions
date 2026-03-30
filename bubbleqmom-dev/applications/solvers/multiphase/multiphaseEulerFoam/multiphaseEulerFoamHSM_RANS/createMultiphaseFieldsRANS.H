/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
-------------------------------------------------------------------------------
Description
    Creates multiphase fields for RANS simulations with time averaging.
    Extended with population balance rate averaging support.
\*---------------------------------------------------------------------------*/

Info<< "Creating multiphase RANS fields\n" << endl;

// Get number of phases
const label nPhases = phases.size();

// Storage for means of phase velocities
PtrList<volVectorField> UMean_(nPhases);

forAll(phases, phasei)
{
    const word phaseName = phases[phasei].name();

    // ========== TIME-AVERAGED PHASE VELOCITIES ==========
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
}

// ========== MODELLED/SGS FIELDS PER PHASE ==========
PtrList<volScalarField> kMod_(nPhases);
PtrList<volScalarField> nuMod_(nPhases);
PtrList<volScalarField> omegaMod_(nPhases);

forAll(phases, phasei)
{
    if (phases[phasei].name() == "water")
    {
        const word phaseName = phases[phasei].name();

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
                    IOobject::AUTO_WRITE
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
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("nuMod", dimViscosity, 0.0)
            )
        );

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
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("omegaMod", dimless/dimTime, 0.0)
            )
        );
    }
}

// ========== TIME-AVERAGED MODELLED FIELDS ==========
PtrList<volScalarField> kMod_Mean_(nPhases);
PtrList<volScalarField> omegaMod_Mean_(nPhases);
PtrList<volScalarField> nut_Mean_(nPhases);

forAll(phases, phasei)
{
    if (phases[phasei].name() == "water")
    {
        const word phaseName = phases[phasei].name();

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
                    IOobject::AUTO_WRITE
                ),
                omegaMod_[phasei]
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
    }
}

// ========== MIXTURE PROPERTIES ==========
PtrList<volScalarField> alphaMean_(nPhases);
forAll(phases, phasei)
{
    if (phases[phasei].name() == "air")
    {
        const word phaseName = phases[phasei].name();

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

// ========== ISOVOLUME CALCULATIONS ==========
// Create isovolume fields for air phase above certain thresholds
volScalarField isoVol_80
(
    IOobject
    (
        "isoVol_air_0.8",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
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
        IOobject::AUTO_WRITE
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
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("isoVol_95", dimless, 0.0)
);

// ========== CBETA FIELD ==========
volScalarField CBeta
(
    IOobject
    (
        "CBeta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("CBeta", dimless, 0.0)
);

// Time-averaged CBeta
volScalarField CBeta_Mean_
(
    IOobject
    (
        "CBeta_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    CBeta
);

// ========== D32 FIELD ==========
// IMPORTANT: d_32 has dimension [length], not dimensionless!
volScalarField d_32
(
    IOobject
    (
        "d_32",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE  // Only time-averaged version is written
    ),
    mesh,
    dimensionedScalar("d_32", dimLength, 0.0)  // <-- CORRECT: dimLength
);

// Time-averaged d_32
volScalarField d_32_Mean_
(
    IOobject
    (
        "d_32_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("d_32_Mean_", dimLength, 0.0)  // <-- Explicit dimensions
);

Info<< "Multiphase RANS fields created for "
    << nPhases << " phases\n" << endl;
