// ************************************************************************* //
//
// Create fields in the air region
//
// ************************************************************************* //

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

dimensionedScalar erD
(
    physicalProperties.lookup("erD")
);

// Convergence variables
scalar iterCount=0.0;
int convVolt=1;
int convVoltD=1;
int convVoltA=1;

Info<< "Reading field volt\n" << endl;
volScalarField volt
(
    IOobject
    (
        "volt",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field voltD\n" << endl;
volScalarField voltD
(
    IOobject
    (
        "voltD",
        runTime.timeName(),
        meshD,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshD
);

Info<< "Reading field voltA\n" << endl;
volScalarField voltA
(
    IOobject
    (
        "voltA",
        runTime.timeName(),
        meshA,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshA
);

volScalarField epsilon
(
    IOobject
    (
        "epsilon",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("temp",dimless,scalar(1))
);

volScalarField epsilonD
(
    IOobject
    (
        "epsilonD",
        runTime.timeName(),
        meshD,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    meshD,
    erD
);

volScalarField epsilonA
(
    IOobject
    (
        "epsilonA",
        runTime.timeName(),
        meshA,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    meshA,
    dimensionedScalar("temp",dimless,scalar(1))
);

Info<< "Calculating field E\n" << endl;
volVectorField E
(
    IOobject
    (
        "E",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(volt)
);

Info<< "Calculating field ED\n" << endl;
volVectorField ED
(
    IOobject
    (
        "ED",
        runTime.timeName(),
        meshD,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(voltD)
);

Info<< "Calculating field EA\n" << endl;
volVectorField EA
(
    IOobject
    (
        "EA",
        runTime.timeName(),
        meshA,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(voltA)
);
