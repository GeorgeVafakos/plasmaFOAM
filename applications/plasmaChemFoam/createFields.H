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

// Boltzmann constant
dimensionedScalar kB
(
    physicalProperties.lookup("kB")
);

// Electron charge
dimensionedScalar e
(
    physicalProperties.lookup("e")
);

// Ambient Pressure
dimensionedScalar press
(
    physicalProperties.lookup("press")
);

// Ambient temperature
dimensionedScalar Tamb
(
    physicalProperties.lookup("Tamb")
);

// Electron temperature
dimensionedScalar Te
(
    physicalProperties.lookup("Te")
);

// Recombination coefficient
dimensionedScalar rnp
(
    physicalProperties.lookup("rnp")
);

// Recombination coefficient
dimensionedScalar rep
(
    physicalProperties.lookup("rep")
);

// Dielectric elative permittivity
dimensionedScalar erD
(
    physicalProperties.lookup("erD")
);

// Insulator elative permittivity
dimensionedScalar erI
(
    physicalProperties.lookup("erI")
);

// Convergence variables
int convVoltAext=1;
int convVoltDext=1;
int convVoltIext=1;
int convVoltAind=1;
int convVoltDind=1;
int convVoltIind=1;
int counter=0;

// scalar cond=1.0e-3;

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

// Insulator external electric potential
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

// Air induced electric potential
Info<< "Reading field voltAind\n" << endl;
volScalarField voltAind
(
    IOobject
    (
        "voltAind",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Dielectric electric potential
Info<< "Reading field voltDind\n" << endl;
volScalarField voltDind
(
    IOobject
    (
        "voltDind",
        runTime.timeName(),
        meshD,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshD
);

// Insulator electric potential
Info<< "Reading field voltIind\n" << endl;
volScalarField voltIind
(
    IOobject
    (
        "voltIind",
        runTime.timeName(),
        meshI,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshI
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

// Insulator total potential
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

// Insulator external electric field
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

// Air induced electric field
Info<< "Calculating field EAind\n" << endl;
volVectorField EAind
(
    IOobject
    (
        "EAind",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(voltAind)
);

// Dielectric induced electric field
Info<< "Calculating field EDind\n" << endl;
volVectorField EDind
(
    IOobject
    (
        "EDind",
        runTime.timeName(),
        meshD,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(voltDind)
);

// Insulator induced electric field
Info<< "Calculating field EIind\n" << endl;
volVectorField EIind
(
    IOobject
    (
        "EIind",
        runTime.timeName(),
        meshI,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(voltIind)
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

// Insulator total electric field
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

// Insulator relative permittivity field
volScalarField epsilonI
(
    IOobject
    (
        "epsilonI",
        runTime.timeName(),
        meshI,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    meshI,
    erI
);







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

// Positive ion number density
Info<< "Reading field np\n" << endl;
volScalarField np
(
    IOobject
    (
        "np",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Negative ion number density
Info<< "Reading field nn\n" << endl;
volScalarField nn
(
    IOobject
    (
        "nn",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Electron mobility
Info<< "Calculating field mue\n" << endl;
volScalarField mue
(
    IOobject
    (
        "mue",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("temp",dimensionSet(-1, 0, 2, 0, 0, 1, 0),scalar(0))
);

// Positive ion mobility
Info<< "Calculating field mup\n" << endl;
volScalarField mup
(
    IOobject
    (
        "mup",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("temp",dimensionSet(-1, 0, 2, 0, 0, 1, 0),scalar(0))
);

// Negative ion mobility
Info<< "Calculating field mun\n" << endl;
volScalarField mun
(
    IOobject
    (
        "mun",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("temp",dimensionSet(-1, 0, 2, 0, 0, 1, 0),scalar(0))
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
    -linearInterpolate(mue)*mesh.magSf()*fvc::snGrad(voltA)
);

// Positive ion flux
Info<< "Calculating field posFlux\n" << endl;
surfaceScalarField posFlux
(
    IOobject
    (
        "posFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    -linearInterpolate(mup)*mesh.magSf()*fvc::snGrad(voltA)
);

// Negative ion flux
Info<< "Calculating field negFlux\n" << endl;
surfaceScalarField negFlux
(
    IOobject
    (
        "negFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    -linearInterpolate(mun)*mesh.magSf()*fvc::snGrad(voltA)
);

// Ionisation coefficient
Info<< "Calculating field alpha\n" << endl;
volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("temp",dimensionSet(0,-1,0,0,0,0,0),scalar(0))
);

// Recombination coefficient
Info<< "Calculating field heta\n" << endl;
volScalarField heta
(
    IOobject
    (
        "heta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("temp",dimensionSet(0,-1,0,0,0,0,0),scalar(0))
);

// Electron diffusion coefficient
Info<< "Calculating field De\n" << endl;
volScalarField De
(
    IOobject
    (
        "De",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("temp",dimensionSet(0,2,-1,0,0,0,0),scalar(0))
);

// Positive Ion diffusion coefficient
Info<< "Calculating field Dp\n" << endl;
volScalarField Dp
(
    IOobject
    (
        "Dp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("temp",dimensionSet(0,2,-1,0,0,0,0),scalar(0))
);

// Negative Ion diffusion coefficient
Info<< "Calculating field Dn\n" << endl;
volScalarField Dn
(
    IOobject
    (
        "Dn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("temp",dimensionSet(0,2,-1,0,0,0,0),scalar(0))
);

// Electron flux vector
Info<< "Calculating field GammaEle\n" << endl;
volVectorField GammaEle
(
    IOobject
    (
        "GammaEle",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -mue*ne*EA - De*fvc::grad(ne)
);

// Positive ion flux vector
Info<< "Calculating field GammaPos\n" << endl;
volVectorField GammaPos
(
    IOobject
    (
        "GammaPos",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mup*np*EA - Dp*fvc::grad(np)
);

// Negative ion flux vector
Info<< "Calculating field GammaNeg\n" << endl;
volVectorField GammaNeg
(
    IOobject
    (
        "GammaNeg",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -mun*nn*EA - Dn*fvc::grad(nn)
);

// Electric field dimentions
Info<< "Calculating field dimE\n" << endl;
volScalarField dimE
(
    IOobject
    (
        "dimE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    // dimensionedScalar("temp",dimless,scalar(0))
    dimensionedScalar("temp",dimensionSet(1,1,-3,0,0,-1,0),scalar(1))
);

// Voltage dimentions
Info<< "Calculating field dimVolt\n" << endl;
volScalarField dimVolt
(
    IOobject
    (
        "dimVolt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    // dimensionedScalar("temp",dimless,scalar(0))
    dimensionedScalar("temp",dimensionSet(1,2,-3,0,0,-1,0),scalar(1))
);

// Ionisation coefficient dimentions
Info<< "Calculating field dimAlpha\n" << endl;
volScalarField dimAlpha
(
    IOobject
    (
        "dimAlpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    // dimensionedScalar("temp",dimless,scalar(0))
    dimensionedScalar("temp",dimensionSet(0,-1,0,0,0,0,0),scalar(1))
);

// Mobility dimentions
Info<< "Calculating field dimMu\n" << endl;
volScalarField dimMu
(
    IOobject
    (
        "dimMu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    // dimensionedScalar("temp",dimless,scalar(0))
    dimensionedScalar("temp",dimensionSet(-1,0,2,0,0,1,0),scalar(1))
);

