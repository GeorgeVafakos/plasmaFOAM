// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Create fields in the solid regions
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


dimensionedScalar erD
(
    physicalProperties.lookup("erD")
);

dimensionedScalar erI
(
    physicalProperties.lookup("erI")
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

Info<< "Reading field voltIextMag\n" << endl;
volScalarField voltIextMag
(
    IOobject
    (
        "voltIextMag",
        runTime.timeName(),
        meshI,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshI
);

Info<< "Reading field voltIrho\n" << endl;
volScalarField voltIrho
(
    IOobject
    (
        "voltIrho",
        runTime.timeName(),
        meshI,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshI
);

Info<< "Reading field voltI\n" << endl;
volScalarField voltI
(
    IOobject
    (
        "voltI",
        runTime.timeName(),
        meshI,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    meshI,
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

volScalarField epsilonI
(
    IOobject
    (
        "epsilonI",
        runTime.timeName(),
        meshI,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    meshI,
    erI
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

Info<< "Calculating field EIextMag\n" << endl;
volVectorField EIextMag
(
    IOobject
    (
        "EIextMag",
        runTime.timeName(),
        meshI,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(voltIextMag)
);

Info<< "Calculating field EIrho\n" << endl;
volVectorField EIrho
(
    IOobject
    (
        "EIrho",
        runTime.timeName(),
        meshI,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(voltIrho)
);

Info<< "Calculating field EI\n" << endl;
volVectorField EI
(
    IOobject
    (
        "EI",
        runTime.timeName(),
        meshI,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(voltI)
);
