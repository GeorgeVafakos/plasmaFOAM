// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Create fields in the air region
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Info<< "Reading physicalProperties\n" << endl;

IOdictionary physicalProperties
(
    IOobject
    (
        "physicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar e0
(
    physicalProperties.lookup("e0")
);

// Read streamer time and frequency
scalar streamerTime =
    runTime.controlDict().lookupOrDefault<scalar>("streamerTime", 0.5e-05);

scalar freq =
    runTime.controlDict().lookupOrDefault<scalar>("freq", 1.0e04);

// Convergence variables
int convVoltAextMag=1;
int convVoltArho=1;
int convVoltDextMag=1;
int convVoltDrho=1;
int convVoltIextMag=1;
int convVoltIrho=1;
int counter=0;

Info<< "Reading field voltAextMag\n" << endl;
volScalarField voltAextMag
(
    IOobject
    (
        "voltAextMag",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field voltArho\n" << endl;
volScalarField voltArho
(
    IOobject
    (
        "voltArho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field voltA\n" << endl;
volScalarField voltA
(
    IOobject
    (
        "voltA",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("temp",dimensionSet(1,2,-3,0,0,-1,0),scalar(0))
);

volScalarField epsilonA
(
    IOobject
    (
        "epsilonA",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("temp",dimless,scalar(1))
);

Info<< "Reading field rhoq\n" << endl;
volScalarField rhoq
(
    IOobject
    (
        "rhoq",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Calculating field EAextMag\n" << endl;
volVectorField EAextMag
(
    IOobject
    (
        "EAextMag",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(voltAextMag)
);

Info<< "Calculating field EArho\n" << endl;
volVectorField EArho
(
    IOobject
    (
        "EArho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(voltArho)
);

Info<< "Calculating field EA\n" << endl;
volVectorField EA
(
    IOobject
    (
        "EA",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(voltA)
);

Info<< "Calculating field Fc\n" << endl;
volVectorField Fc
(
    IOobject
    (
        "Fc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rhoq*EA
);
