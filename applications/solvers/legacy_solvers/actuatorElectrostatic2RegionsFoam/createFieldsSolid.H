// ************************************************************************* //
//
// Create fields in the solid regions
//
// ************************************************************************* //


dimensionedScalar erD
(
    physicalProperties.lookup("erD")
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
