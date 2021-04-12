  Info<< "Reading physicalProperties\n" << endl;

    IOdictionary physicalProperties
    (
        IOobject
        (
            "physicalProperties",
            runTime.constant(),
            meshAir,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar epsilon0
    (
        physicalProperties.lookup("epsilon0")
    );

    dimensionedScalar k_loads
    (
        physicalProperties.lookup("k_loads")
    );

// Variables for check of the convergence
int convphiEl=1;
int convrhoq=1;

Info<< "Reading field phiEl\n" << endl;
volScalarField phiEl
(
    IOobject
    (
        "phiEl",
        runTime.timeName(),
        meshAir,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshAir
);

volScalarField k_Air
(
    IOobject
    (
        "k_Air",
        runTime.constant(),
        meshAir,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    meshAir
);

    Info<< "Reading field rhoq\n" << endl;
    volScalarField rhoq
    (
        IOobject
        (
            "rhoq",
            runTime.timeName(),
            meshAir,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        meshAir
    );


    Info<< "Calculating field rhoFlux\n" << endl;
    surfaceScalarField rhoFlux
    (
        IOobject
        (
            "rhoFlux",
            runTime.timeName(),
            meshAir,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        -k_loads*meshAir.magSf()*fvc::snGrad(phiEl)
    );

    Info<< "Calculating field EA\n" << endl;
    volVectorField EA
    (
        IOobject
        (
            "EA",
            runTime.timeName(),
            meshAir,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        -fvc::grad(phiEl)
    );

    Info<< "Calculating electrostatic force\n" << endl;
    volVectorField f_el
    (
        IOobject
        (
            "f_el",
            runTime.timeName(),
            meshAir,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        meshAir
    );



 

