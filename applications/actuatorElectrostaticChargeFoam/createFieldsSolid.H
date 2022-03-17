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

Info<< "Reading field voltI\n" << endl;
volScalarField voltI
(
    IOobject
    (
        "voltI",
        runTime.timeName(),
        meshI,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshI
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
