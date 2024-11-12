// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Create fields in air and solid regions
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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
dimensionedScalar epsilon0
(
    physicalProperties.lookup("epsilon0")
);

// Electron charge
dimensionedScalar e
(
    physicalProperties.lookup("e")
);

// Electron Diffusion coefficient
dimensionedScalar De
(
    physicalProperties.lookup("De")
);

// Ion Diffusion coefficient
dimensionedScalar Di
(
    physicalProperties.lookup("Di")
);

// Electron mobility
dimensionedScalar mue
(
    physicalProperties.lookup("mue")
);

// Ion mobility
dimensionedScalar mui
(
    physicalProperties.lookup("mui")
);

// Recombination coefficient
dimensionedScalar r
(
    physicalProperties.lookup("r")
);

// Ionization coefficient
dimensionedScalar alpha
(
    physicalProperties.lookup("alpha")
);

// Dielectric elative permittivity
dimensionedScalar erD
(
    physicalProperties.lookup("erD")
);

// Convergence variables
int convVoltAext=1;
int convVoltDext=1;
int convVoltArho=1;
int convVoltDrho=1;
int counter=0;

// Electron number density
Info<< "Reading field ne\n" << endl;
volScalarField ne
(
    IOobject
    (
        "ne",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Ion number density
Info<< "Reading field ni\n" << endl;
volScalarField ni
(
    IOobject
    (
        "ni",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Air external electric potential
Info<< "Reading field voltAext\n" << endl;
volScalarField voltAext
(
    IOobject
    (
        "voltAext",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Dielectric external electric potential
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

// Air induced electric potential
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

// Dielectric electric potential
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

// Air total electric potential
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

// Dielectric total potential
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

// Electron flux
Info<< "Calculating field eleFlux\n" << endl;
surfaceScalarField eleFlux
(
    IOobject
    (
        "eleFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    -mue*mesh.magSf()*fvc::snGrad(voltA)
);

// Ion flux
Info<< "Calculating field ionFlux\n" << endl;
surfaceScalarField ionFlux
(
    IOobject
    (
        "ionFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    -mui*mesh.magSf()*fvc::snGrad(voltA)
);

// Air external electric field
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

// Dielectric external electric field
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

// Air induced electric field
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

// Dielectric induced electric field
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

// Air total electric field
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

// Dielectric total electric field
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

// Air relative permittivity field
volScalarField epsilonA
(
    IOobject
    (
        "epsilonA",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("temp",dimless,scalar(1))
);

// Dielectric relative permittivity field
volScalarField epsilonD
(
    IOobject
    (
        "epsilonD",
        runTime.timeName(),
        meshD,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    meshD,
    erD
);

// // Ionization frequency
// volScalarField ionizFreq
// (
//     IOobject
//     (
//         "ionizFreq",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("temp",dimensionSet(0,0,-1,0,0,0,0),scalar(0))
// );
