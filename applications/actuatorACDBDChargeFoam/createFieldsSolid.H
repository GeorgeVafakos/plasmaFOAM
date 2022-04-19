// ************************************************************************* //
//
// Create fields in the solid regions
//
// ************************************************************************* //


dimensionedScalar erD
(
    physicalProperties.lookup("erD")
);

dimensionedScalar erI
(
    physicalProperties.lookup("erI")
);

Info<< "Reading field voltDext\n" << endl;
volScalarField voltDext
(
    IOobject
    (
        "voltDext",
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

Info<< "Reading field voltDrho_strm\n" << endl;
volScalarField voltDrho_strm
(
    IOobject
    (
        "voltDrho_strm",
        runTime.timeName(),
        meshD,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    voltDrho
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

Info<< "Reading field voltIext\n" << endl;
volScalarField voltIext
(
    IOobject
    (
        "voltIext",
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

Info<< "Reading field voltIrho_strm\n" << endl;
volScalarField voltIrho_strm
(
    IOobject
    (
        "voltIrho_strm",
        runTime.timeName(),
        meshI,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    voltIrho
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

Info<< "Calculating field EDext\n" << endl;
volVectorField EDext
(
    IOobject
    (
        "EDext",
        runTime.timeName(),
        meshD,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(voltDext)
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

Info<< "Calculating field EIext\n" << endl;
volVectorField EIext
(
    IOobject
    (
        "EIext",
        runTime.timeName(),
        meshI,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(voltIext)
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
