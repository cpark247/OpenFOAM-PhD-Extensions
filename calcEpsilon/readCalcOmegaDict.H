    Info<< "Reading calcOmegaDict\n" << endl;

    IOdictionary calcOmegaDict
    (
        IOobject
        (
            "calcOmegabDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    //scalar Ea(boxTurbDict.lookup<scalar>("Ea"));

    //scalar k0(boxTurbDict.lookup<scalar>("k0"));
