#include "createDielectricFields.H"

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

// Vacuum permittivity
dimensionedScalar e0
(
    "e0",
    dimensionSet(-1, -3, 4, 0, 0, 2, 0),
    physicalProperties
);

// Ion mobility
dimensionedScalar muc
(
    "muc",
    dimensionSet(-1, 0, 2, 0, 0, 1, 0),
    physicalProperties
);

// Ion diffusivity
dimensionedScalar Dc
(
    "Dc",
    dimensionSet(0, 2, -1, 0, 0, 0, 0),
    physicalProperties
);

// Convergence variables
int regionLoopCounter=0;
std::vector<int> voltIndIter = {1};

Info<< "Reading field voltExtAmp\n" << endl;
volScalarField voltExtAmp
(
    IOobject
    (
        "voltExtAmp",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field voltExt\n" << endl;
volScalarField voltExt
(
    IOobject
    (
        "voltExt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    voltExtAmp*Foam::sin(2*M_PI*(1.0/endTime)*runTime.value())
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

Info<< "Reading field voltIndDischarge\n" << endl;
volScalarField voltIndDischarge
(
    IOobject
    (
        "voltIndDischarge",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    voltInd
);

Info<< "Reading field volt\n" << endl;
volScalarField volt
(
    IOobject
    (
        "volt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("temp", dimensionSet(1,2,-3,0,0,-1,0), scalar(0))
);

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

Info<< "Reading field rhoqDischarge\n" << endl;
volScalarField rhoqDischarge
(
    IOobject
    (
        "rhoqDischarge",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rhoq
);

Info<< "Calculating field rhoqFlux\n" << endl;
surfaceScalarField rhoqFlux
(
    IOobject
    (
        "rhoqFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    -sign(linearInterpolate(rhoq))*muc*mesh.magSf()*fvc::snGrad(volt)
);

Info<< "Calculating field EExt\n" << endl;
volVectorField EExt
(
    IOobject
    (
        "EExt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(voltExt)
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
    rhoq*E
);


