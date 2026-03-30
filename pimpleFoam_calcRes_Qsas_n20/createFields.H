#include "createRDeltaT.H"

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField UPrime_
(
    IOobject
    (
        "UPrime_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("UPrime_", dimVelocity, vector(0.0,0.0,0.0))
);

volVectorField UMean_
(
    IOobject
    (
        "UMean_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("UMean_", dimVelocity, vector(0.0,0.0,0.0))
);

IOdictionary transportProperties
(
    IOobject
    (
       "transportProperties",
   	runTime.constant(),
   	mesh,
   	IOobject::MUST_READ,
   	IOobject::NO_WRITE
    )
);

dimensionedScalar nu
(
    transportProperties.lookup("nu")
);


///// epsilon /////
volScalarField epsilonRes_
(
    IOobject
    (
        "epsilonRes_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("epsilonRes_", dimensionSet(0,2,-3,0,0,0,0), 0.0)
);

volScalarField epsilonRes_Mean_
(
    IOobject
    (
        "epsilonRes_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("epsilonRes_Mean_", dimensionSet(0,2,-3,0,0,0,0), 0.0)
);

volScalarField epsilonResControl_
(
    IOobject
    (
        "epsilonResControl_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("epsilonResControl_", dimensionSet(0,2,-3,0,0,0,0), 0.0)
);

volScalarField epsilonMod_
(
    IOobject
    (
        "epsilonMod_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("epsilonMod_", dimVelocity*dimVelocity/dimTime, 0.0)
);

volScalarField epsilonMod_Mean_
(
    IOobject
    (
        "epsilonMod_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("epsilonMod_Mean_", dimVelocity*dimVelocity/dimTime, 0.0)
);

volScalarField epsilonTot_
(
    IOobject
    (
        "epsilonTot_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("epsilonTot_", dimensionSet(0,2,-3,0,0,0,0), 0.0)
);

volScalarField epsilonTot_Mean_
(
    IOobject
    (
        "epsilonTot_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("epsilonTot_Mean_", dimensionSet(0,2,-3,0,0,0,0), 0.0)
);

////// TKE //////
volScalarField kRes_
(
    IOobject
    (
        "kRes_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("kRes_", dimVelocity*dimVelocity, 0.0)
);

volScalarField kRes_Mean_
(
    IOobject
    (
        "kRes_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("kRes_Mean_", dimVelocity*dimVelocity, 0.0)
);


volScalarField kMod_
(
    IOobject
    (
        "kMod_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE //same as k printed from kOmegaSSTSAS model
    ),
    mesh,
    dimensionedScalar("kMod_", dimVelocity*dimVelocity, 0.0)
);

volScalarField kMod_Mean_
(
    IOobject
    (
        "kMod_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("kMod_Mean_", dimVelocity*dimVelocity, 0.0)
);


volScalarField kTot_
(
    IOobject
    (
        "kTot_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("kTot_", dimVelocity*dimVelocity, 0.0)
);

volScalarField kTot_Mean_
(
    IOobject
    (
        "kTot_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("kTot_Mean_", dimVelocity*dimVelocity, 0.0)
);

///// omega /////
volScalarField omegaMod_
(
    IOobject
    (
        "omegaMod_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE // same as omega of kOmegaSSTSAS model
    ),
    mesh,
    dimensionedScalar("omegaMod_", dimTime/dimTime/dimTime, 0.0)
);

volScalarField omegaMod_Mean_
(
    IOobject
    (
        "omegaMod_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("omegaMod_Mean_", dimTime/dimTime/dimTime, 0.0)
);

volScalarField omegaRes_
(
    IOobject
    (
        "omegaRes_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("omegaRes_", dimTime/dimTime/dimTime, 0.0)
);

volScalarField omegaRes_Mean_
(
    IOobject
    (
        "omegaRes_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("omegaRes_Mean_", dimTime/dimTime/dimTime, 0.0)
);

volScalarField omegaTot_
(
    IOobject
    (
        "omegaTot_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("omegaTot_", dimTime/dimTime/dimTime, 0.0)
);

volScalarField omegaTot_Mean_
(
    IOobject
    (
        "omegaTot_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("omegaTot_Mean_", dimTime/dimTime/dimTime, 0.0)
);

/*
Info<< "Reading field sigmaMod_\n" << endl;
volSymmTensorField sigmaMod_
(
     IOobject
     (
        "sigmaMod_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
     ),
    mesh,
    dimensionedSymmTensor("sigmaMod_", dimVelocity*dimVelocity, symmTensor::zero)
);
*/

///// eddy visc /////
volScalarField nuMod_
(
    IOobject
    (
        "nuMod_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE // same as nut from kOmegaSSTSAS model 
    ),
    mesh,
    dimensionedScalar("nuMod_", dimVelocity*dimVelocity*dimTime, 0.0)
);

volScalarField nut_Mean_
(
    IOobject
    (
        "nut_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE 
    ),
    mesh,
    dimensionedScalar("nut_Mean_", dimVelocity*dimVelocity*dimTime, 0.0)
);



volScalarField u2_
(
    IOobject
    (
        "u2_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("u2_", dimVelocity*dimVelocity, 0.0)
);

volScalarField v2_
(
    IOobject
    (
        "v2_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("v2_", dimVelocity*dimVelocity, 0.0)
);

volScalarField w2_
(
    IOobject
    (
        "w2_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("w2_", dimVelocity*dimVelocity, 0.0)
);

volScalarField uPrime2_
(
    IOobject
    (
        "uPrime2_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("uPrime2_", dimVelocity*dimVelocity, 0.0)
);

volScalarField vPrime2_
(
    IOobject
    (
        "vPrime2_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("vPrime2_", dimVelocity*dimVelocity, 0.0)
);

volScalarField wPrime2_
(
    IOobject
    (
        "wPrime2_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("wPrime2_", dimVelocity*dimVelocity, 0.0)
);

volScalarField Qsas_
(
    IOobject
    (
        "Qsas_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Qsas_", dimTime/dimTime/dimTime/dimTime, 0.0)
);

volScalarField Qsas_Mean_
(
    IOobject
    (
        "Qsas_Mean_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Qsas_Mean_", dimTime/dimTime/dimTime/dimTime, 0.0)
);

#include "createPhi.H"

pressureReference pressureReference(p, pimple.dict());

mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"
