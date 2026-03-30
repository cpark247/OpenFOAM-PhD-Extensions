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
        dimVelocity*dimVelocity/dimTime, 0.0)
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
// U squared for rms
// ============================================================================
volSymmTensorField UPrime2_
(
    IOobject
    (
        "UPrime2_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedSymmTensor("zero", dimVelocity*dimVelocity, symmTensor::zero)
);


// UPrime squared Mean - standard deviation
volSymmTensorField UPrime2_Mean_
(
    IOobject
    (
        "UPrime2_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    UPrime2_
);

