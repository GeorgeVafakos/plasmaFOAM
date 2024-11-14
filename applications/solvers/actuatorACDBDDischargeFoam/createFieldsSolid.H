// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Create fields in the solid regions
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensionedScalar erD
(
    "erD",
    dimless,
    physicalProperties
);

dimensionedScalar erE
(
    "erE",
    dimless,
    physicalProperties
);

Info<< "Reading field voltDextMag\n" << endl;
volScalarField voltDextMag
(
    IOobject
    (
        "voltDextMag",
        runTime.timeName(),
        meshD,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshD
);

Info<< "Reading field voltDrho\n" << endl;
volScalarField voltDrho
(
    IOobject
    (
        "voltDrho",
        runTime.timeName(),
        meshD,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshD
);

Info<< "Reading field voltD\n" << endl;
volScalarField voltD
(
    IOobject
    (
        "voltD",
        runTime.timeName(),
        meshD,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    meshD,
    dimensionedScalar("temp",dimensionSet(1,2,-3,0,0,-1,0),scalar(0))
);

Info<< "Reading field voltEextMag\n" << endl;
volScalarField voltEextMag
(
    IOobject
    (
        "voltEextMag",
        runTime.timeName(),
        meshE,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshE
);

Info<< "Reading field voltErho\n" << endl;
volScalarField voltErho
(
    IOobject
    (
        "voltErho",
        runTime.timeName(),
        meshE,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshE
);

Info<< "Reading field voltE\n" << endl;
volScalarField voltE
(
    IOobject
    (
        "voltE",
        runTime.timeName(),
        meshE,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    meshE,
    dimensionedScalar("temp",dimensionSet(1,2,-3,0,0,-1,0),scalar(0))
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

volScalarField epsilonE
(
    IOobject
    (
        "epsilonE",
        runTime.timeName(),
        meshE,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    meshE,
    erE
);

Info<< "Calculating field EDextMag\n" << endl;
volVectorField EDextMag
(
    IOobject
    (
        "EDextMag",
        runTime.timeName(),
        meshD,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(voltDextMag)
);

Info<< "Calculating field EDrho\n" << endl;
volVectorField EDrho
(
    IOobject
    (
        "EDrho",
        runTime.timeName(),
        meshD,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(voltDrho)
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

Info<< "Calculating field EEextMag\n" << endl;
volVectorField EEextMag
(
    IOobject
    (
        "EEextMag",
        runTime.timeName(),
        meshE,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(voltEextMag)
);

Info<< "Calculating field EErho\n" << endl;
volVectorField EErho
(
    IOobject
    (
        "EErho",
        runTime.timeName(),
        meshE,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(voltErho)
);

Info<< "Calculating field EE\n" << endl;
volVectorField EE
(
    IOobject
    (
        "EE",
        runTime.timeName(),
        meshE,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(voltE)
);
