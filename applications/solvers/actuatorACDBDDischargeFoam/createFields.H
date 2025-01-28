#include "createDielectricFields.H"

// Vacuum permittivity
dimensionedScalar e0
(
    "e0",
    dimensionSet(-1, -3, 4, 0, 0, 2, 0),
    scalar(8.85419e-12)
);

// Region convergence variables
int regionLoopCounter=0;
std::vector<int> voltEqnIter = {1};

// Read controlDict for maximum region loops
int maxRegionLoops =
    runTime.controlDict().lookupOrDefault("maxRegionLoops", 30000);

// Read controlDict for writing current solution
int writeRegionLoop =
    runTime.controlDict().lookupOrDefault("writeRegionLoop", 5000);

volScalarField epsilon
(
    IOobject
    (
        "epsilon",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("temp", dimless, scalar(1))
);

Info<< "Reading field voltInd\n" << endl;
volScalarField voltInd
(
    IOobject
    (
        "voltInd",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
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

Info<< "Calculating field EInd\n" << endl;
volVectorField EInd
(
    IOobject
    (
        "EInd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(voltInd)
);
