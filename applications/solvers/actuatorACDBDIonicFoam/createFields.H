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
int convVoltArho=1;
int convVoltDrho=1;
int convVoltErho=1;
int convRhoq=1;
int counter=0;

Info<< "Reading field voltAextMag\n" << endl;
volScalarField voltAextMag
(
    IOobject
    (
        "voltAextMag",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Reading field voltAext\n" << endl;
volScalarField voltAext
(
    IOobject
    (
        "voltAext",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    //mesh
    voltAextMag*Foam::sin(2*M_PI*(1.0/endTime)*runTime.value())
);

Info<< "Reading field voltArho\n" << endl;
volScalarField voltArho
(
    IOobject
    (
        "voltArho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field voltArho_strm\n" << endl;
volScalarField voltArho_strm
(
    IOobject
    (
        "voltArho_strm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    voltArho
);

Info<< "Reading field voltA\n" << endl;
volScalarField voltA
(
    IOobject
    (
        "voltA",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("temp",dimensionSet(1,2,-3,0,0,-1,0),scalar(0))
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

Info<< "Reading field rhoq_strm\n" << endl;
volScalarField rhoq_strm
(
    IOobject
    (
        "rhoq_strm",
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
    -sign(linearInterpolate(rhoq))*muc*mesh.magSf()*fvc::snGrad(voltA)
);

Info<< "Calculating field EAext\n" << endl;
volVectorField EAext
(
    IOobject
    (
        "EAext",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(voltAext)
);

Info<< "Calculating field EArho\n" << endl;
volVectorField EArho
(
    IOobject
    (
        "EArho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(voltArho)
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
    rhoq*EA
);


