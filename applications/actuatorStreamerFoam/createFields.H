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

dimensionedScalar k
(
    physicalProperties.lookup("k")
);

// Convergence variables
scalar iterCount=0.0;
int convVoltA=1;
int convVoltD=1;
int convVoltI=1;
int convVoltR=1;

Info<< "Reading field voltA\n" << endl;
volScalarField voltA
(
    IOobject
    (
        "voltA",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
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

Info<< "Reading field voltR\n" << endl;
volScalarField voltR
(
    IOobject
    (
        "voltR",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Calculating field ER\n" << endl;
volVectorField ER
(
    IOobject
    (
        "ER",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(voltR)
);
