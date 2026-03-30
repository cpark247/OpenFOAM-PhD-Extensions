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

Info<< "Creating resolved turbulence fields\n" << endl;

// ============================================================================
// TIME-AVERAGED PRESSURE
// ============================================================================
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
// VELOCITY (MEAN) + FLUCTUATIONS (TRANSIENT ONLY)
// ============================================================================
volVectorField UMean_
(
    IOobject
    (
        "UMean_",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    U
);

// Container for fluctuating velocity (U') – transient only, not written
volVectorField UPrime_
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
);

// ============================================================================
// MODELLED/SGS FIELDS (INSTANTANEOUS)
// ============================================================================
volScalarField kMod_
(
    IOobject
    (
        "kMod",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("kMod", dimVelocity*dimVelocity, 0.0)
);

volScalarField nuMod_
(
    IOobject
    (
        "nuMod",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("nuMod", dimViscosity, 0.0)
);

volScalarField epsilonMod_
(
    IOobject
    (
        "epsilonMod",
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
);

// ============================================================================
// MODELLED/SGS FIELDS (MEAN)
// ============================================================================
volScalarField kMod_Mean_
(
    IOobject
    (
        "kMod_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    kMod_
);

volScalarField epsilonMod_Mean_
(
    IOobject
    (
        "epsilonMod_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    epsilonMod_
);

volScalarField nut_Mean_
(
    IOobject
    (
        "nut_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    nuMod_
);

// ============================================================================
// RESOLVED FIELDS (INSTANTANEOUS)
// ============================================================================
volScalarField kRes_
(
    IOobject
    (
        "kRes",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("kRes", dimVelocity*dimVelocity, 0.0)
);

volScalarField epsilonRes_
(
    IOobject
    (
        "epsilonRes",
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
);

// ============================================================================
// RESOLVED FIELDS (MEAN)
// ============================================================================
volScalarField kRes_Mean_
(
    IOobject
    (
        "kRes_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    kRes_
);

volScalarField epsilonRes_Mean_
(
    IOobject
    (
        "epsilonRes_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    epsilonRes_
);

// ============================================================================
// TOTAL FIELDS (INSTANTANEOUS)
// ============================================================================
volScalarField kTotal_
(
    IOobject
    (
        "kTotal",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("kTotal", dimVelocity*dimVelocity, 0.0)
);

volScalarField epsilonTotal_
(
    IOobject
    (
        "epsilonTotal",
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
);

// ============================================================================
// TOTAL FIELDS (MEAN)
// ============================================================================
volScalarField kTotal_Mean_
(
    IOobject
    (
        "kTotal_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    kTotal_
);

volScalarField epsilonTotal_Mean_
(
    IOobject
    (
        "epsilonTotal_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    epsilonTotal_
);

// ============================================================================
// REYNOLDS STRESS COMPONENTS (INSTANTANEOUS)
// ============================================================================
volScalarField uu_
(
    IOobject
    (
        "uu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("uu", dimVelocity*dimVelocity, 0.0)
);

volScalarField vv_
(
    IOobject
    (
        "vv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("vv", dimVelocity*dimVelocity, 0.0)
);

volScalarField ww_
(
    IOobject
    (
        "ww",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("ww", dimVelocity*dimVelocity, 0.0)
);

volScalarField uv_
(
    IOobject
    (
        "uv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("uv", dimVelocity*dimVelocity, 0.0)
);

volScalarField uw_
(
    IOobject
    (
        "uw",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("uw", dimVelocity*dimVelocity, 0.0)
);

volScalarField vw_
(
    IOobject
    (
        "vw",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("vw", dimVelocity*dimVelocity, 0.0)
);

// Time-averaged Reynolds stress components
volScalarField uu_Mean_
(
    IOobject
    (
        "uu_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    uu_
);

volScalarField vv_Mean_
(
    IOobject
    (
        "vv_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    vv_
);

volScalarField ww_Mean_
(
    IOobject
    (
        "ww_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    ww_
);

volScalarField uv_Mean_
(
    IOobject
    (
        "uv_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    uv_
);

volScalarField uw_Mean_
(
    IOobject
    (
        "uw_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    uw_
);

volScalarField vw_Mean_
(
    IOobject
    (
        "vw_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    vw_
);

// ============================================================================
// LES GRID RESOLUTION ASSESSMENT
// ============================================================================
volScalarField deltaByKolmog_
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
    deltaByKolmog_
);

Info<< "Resolved turbulence fields created\n" << endl;
