
volScalarField kMod_
(
    IOobject
    (
        "kMod_",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField uPrime_uPrime_Mean_
(
    IOobject
    (
        "uPrime_uPrime_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField uPrime_vPrime_Mean_
(
    IOobject
    (
        "uPrime_vPrime_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField uPrime_wPrime_Mean_
(
    IOobject
    (
        "uPrime_wPrime_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField vPrime_vPrime_Mean_
(
    IOobject
    (
        "vPrime_vPrime_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField vPrime_wPrime_Mean_
(
    IOobject
    (
        "vPrime_wPrime_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField wPrime_wPrime_Mean_
(
    IOobject
    (
        "wPrime_wPrime_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volSymmTensorField R_
(
    IOobject
    (
        "R_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedSymmTensor("zero", dimensionSet(0,2,-2,0,0,0,0), symmTensor::zero)
);
