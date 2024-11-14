// ************************************************************************* //
//
// Create fields in the solid regions
//
// ************************************************************************* //
// Dielectric material relative permittivity
dimensionedScalar erD
(
    "erD",
    dimless,
    physicalProperties
);

// Encapsulator material relative permittivity
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
        IOobject::NO_WRITE
    ),
    meshD
);

Info<< "Reading field voltDext\n" << endl;
volScalarField voltDext
(
    IOobject
    (
        "voltDext",
        runTime.timeName(),
        meshD,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    voltDextMag*Foam::sin(2*M_PI*(1.0/endTime)*runTime.value())
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

Info<< "Reading field voltEextMag\n" << endl;
volScalarField voltEextMag
(
    IOobject
    (
        "voltEextMag",
        runTime.timeName(),
        meshE,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    meshE
);

Info<< "Reading field voltEext\n" << endl;
volScalarField voltEext
(
    IOobject
    (
        "voltEext",
        runTime.timeName(),
        meshE,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    voltEextMag*Foam::sin(2*M_PI*(1.0/endTime)*runTime.value())
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

Info<< "Reading field voltErho_strm\n" << endl;
volScalarField voltErho_strm
(
    IOobject
    (
        "voltErho_strm",
        runTime.timeName(),
        meshE,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    voltErho
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

Info<< "Calculating field EEext\n" << endl;
volVectorField EEext
(
    IOobject
    (
        "EEext",
        runTime.timeName(),
        meshE,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(voltEext)
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
