volScalarField epsilonMod
(
    IOobject
    (
        "epsilonMod",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("epsilonMod", dimVelocity*dimVelocity/dimTime, 0.0)
);
