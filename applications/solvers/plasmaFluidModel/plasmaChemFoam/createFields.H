// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Create fields in air and solid regions
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


volScalarField ones
(
    IOobject
    ("ones",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("temp",dimensionSet(0,0,0,0,0,0,0),scalar(1))
);

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
    "epsilon0",
    dimensionSet(-1, -3, 4, 0, 0, 2, 0),
    physicalProperties
);

// Boltzmann constant
dimensionedScalar kB
(
    "kB",
    dimensionSet(1, 2, -2, -1, 0, 0, 0),
    physicalProperties
);

// Electron charge
dimensionedScalar e
(
    "e",
    dimensionSet(0, 0, 1, 0, 0, 1, 0),
    physicalProperties
);

// Ambient Pressure
dimensionedScalar press
(
    "press",
    dimensionSet(0, 0, 0, 0, 0, 0, 0),
    physicalProperties
);

// Ambient temperature
dimensionedScalar Tamb
(
    "Tamb",
    dimensionSet(0, 0, 0, 1, 0, 0, 0),
    physicalProperties
);

// Electron temperature
dimensionedScalar Te
(
    "Te",
    dimensionSet(0, 0, 0, 1, 0, 0, 0),
    physicalProperties
);

// Electron Diffusion coefficient
dimensionedScalar De
(
    "De",
    dimensionSet(0, 2, -1, 0, 0, 0, 0),
    physicalProperties
);

// Ion Diffusion coefficient
dimensionedScalar Di
(
    "Di",
    dimensionSet(0, 2, -1, 0, 0, 0, 0),
    physicalProperties
);

// Electron mobility
dimensionedScalar mue
(
    "mue",
    dimensionSet(-1, 0, 2, 0, 0, 1, 0),
    physicalProperties
);
volScalarField mueF
(IOobject("mueF",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
    mesh,
    mue * dimensionedScalar("temp",dimless,scalar(1))
);

// Positive Ion mobility
dimensionedScalar mup
(
    "mup",
    dimensionSet(-1, 0, 2, 0, 0, 1, 0),
    physicalProperties
);
volScalarField mupF(IOobject("mupF",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
    mesh,
    mup * dimensionedScalar("temp",dimless,scalar(1))
);

// Negative Ion mobility
dimensionedScalar mun
(
    "mun",
    dimensionSet(-1, 0, 2, 0, 0, 1, 0),
    physicalProperties
);
volScalarField munF(IOobject("munF",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
    mesh,
    mun * dimensionedScalar("temp",dimless,scalar(1))
);

// Recombination coefficients
dimensionedScalar r("r", dimensionSet(-1, 0, 2, 0, 0, 1, 0), physicalProperties);
dimensionedScalar rnp("rnp", dimensionSet(-1, 0, 2, 0, 0, 1, 0), physicalProperties);
dimensionedScalar rep("rep", dimensionSet(-1, 0, 2, 0, 0, 1, 0), physicalProperties);

// Ionization coefficient
dimensionedScalar alpha("alpha", dimensionSet(0, -1, 0, 0, 0, 0, 0), physicalProperties);
Info<< "Calculating field alphaField\n" << endl;
volScalarField alphaF(IOobject("alphaF", runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
    mesh,
    alpha * dimensionedScalar("temp",dimensionSet(0,0,0,0,0,0,0),scalar(1))
);

// Recombinations rate
volScalarField etaF(IOobject("etaF", runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
    mesh,
    alpha * dimensionedScalar("temp",dimensionSet(0,0,0,0,0,0,0),scalar(1))
);

// Dielectric elative permittivity
dimensionedScalar erD
(
    "erD",
    dimensionSet(0, 0, 0, 0, 0, 0, 0),
    physicalProperties
);

// Convergence variables
int convVoltAext=1;
int convVoltDext=1;
int convVoltArho=1;
int convVoltDrho=1;
int counter=0;

// Electron number density
Info<< "Reading field ne\n" << endl;
volScalarField ne(IOobject("ne",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),
    mesh
);

// Ion number density
Info<< "Reading field np\n" << endl;
volScalarField np(IOobject("np",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),
    mesh
);

// Negative ion number density
Info<< "Reading field nn\n" << endl;
volScalarField nn(IOobject("nn",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),
    mesh
);

// Electron number density in the background
volScalarField ne0(IOobject("ne0",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),
    ne
);


// Air external electric potential
Info<< "Reading field voltAext\n" << endl;
volScalarField voltAext(IOobject("voltAext",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),
    mesh
);

// Dielectric external electric potential
Info<< "Reading field voltDext\n" << endl;
volScalarField voltDext(IOobject("voltDext",runTime.timeName(),meshD,IOobject::MUST_READ,IOobject::AUTO_WRITE),
    meshD
);

// Air induced electric potential
Info<< "Reading field voltArho\n" << endl;
volScalarField voltArho(IOobject("voltArho",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),
    mesh
);

// Dielectric electric potential
Info<< "Reading field voltDrho\n" << endl;
volScalarField voltDrho(IOobject("voltDrho",runTime.timeName(),meshD,IOobject::MUST_READ,IOobject::AUTO_WRITE),
    meshD
);

// Air total electric potential
Info<< "Reading field voltA\n" << endl;
volScalarField voltA(IOobject("voltA",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("temp",dimensionSet(1,2,-3,0,0,-1,0),scalar(0))
);

// Dielectric total potential
Info<< "Reading field voltD\n" << endl;
volScalarField voltD(IOobject("voltD",runTime.timeName(),meshD,IOobject::NO_READ,IOobject::AUTO_WRITE),
    meshD,
    dimensionedScalar("temp",dimensionSet(1,2,-3,0,0,-1,0),scalar(0))
);

// Electron flux
Info<< "Calculating field eleFlux\n" << endl;
surfaceScalarField eleFlux(IOobject("eleFlux",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),
    -mue*mesh.magSf()*fvc::snGrad(voltA)
);

//Positive Ion flux
surfaceScalarField posFlux(IOobject("posFlux",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),
    -mup*mesh.magSf()*fvc::snGrad(voltA)
);

//Negative Ion flux
surfaceScalarField negFlux(IOobject("negFlux",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),
    -mun*mesh.magSf()*fvc::snGrad(voltA)
);

// Air external electric field
Info<< "Calculating field EAext\n" << endl;
volVectorField EAext(IOobject("EAext",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
    -fvc::grad(voltAext)
);

// Dielectric external electric field
Info<< "Calculating field EDext\n" << endl;
volVectorField EDext(IOobject("EDext",runTime.timeName(),meshD,IOobject::NO_READ,IOobject::AUTO_WRITE),
    -fvc::grad(voltDext)
);

// Air induced electric field
Info<< "Calculating field EArho\n" << endl;
volVectorField EArho(IOobject("EArho",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
    -fvc::grad(voltArho)
);

// Dielectric induced electric field
Info<< "Calculating field EDrho\n" << endl;
volVectorField EDrho(IOobject("EDrho",runTime.timeName(),meshD,IOobject::NO_READ,IOobject::AUTO_WRITE),
    -fvc::grad(voltDrho)
);

// Air total electric field
Info<< "Calculating field EA\n" << endl;
volVectorField EA(IOobject("EA",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
    -fvc::grad(voltA)
);

// Dielectric total electric field
Info<< "Calculating field ED\n" << endl;
volVectorField ED(IOobject("ED",runTime.timeName(),meshD,IOobject::NO_READ,IOobject::AUTO_WRITE),
    -fvc::grad(voltD)
);

// Air relative permittivity field
volScalarField epsilonA(IOobject("epsilonA",runTime.timeName(),mesh,IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("temp",dimless,scalar(1))
);

// Dielectric relative permittivity field
volScalarField epsilonD(IOobject("epsilonD",runTime.timeName(),meshD,IOobject::NO_READ,IOobject::NO_WRITE),
    meshD,
    erD
);

// Electron diffusion coefficient
volScalarField DeF(IOobject("DeF",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("temp",dimensionSet(0,2,-1,0,0,0,0),scalar(0))
);

// Positive Ion diffusion coefficient
volScalarField DpF(IOobject("DpF",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("temp",dimensionSet(0,2,-1,0,0,0,0),scalar(0))
);

// Negative Ion diffusion coefficient
volScalarField DnF(IOobject("DnF",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("temp",dimensionSet(0,2,-1,0,0,0,0),scalar(0))
);

// Electron flux vector
volVectorField GammaEle(IOobject("GammaEle",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
    -mueF*ne*EA - DeF*fvc::grad(ne)
);

// Positive ion flux vector
volVectorField GammaPos(IOobject("GammaPos",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
    mupF*np*EA - DpF*fvc::grad(np)
);

// Negative ion flux vector
volVectorField GammaNeg(IOobject("GammaNeg",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
    -munF*nn*EA - DnF*fvc::grad(nn)
);

// Electron sourse term
volScalarField SLe(IOobject("SLe",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
    alphaF*mag(GammaEle) - (rep*ne*np + etaF*mag(GammaEle))
);
// Positive ion sourse term
volScalarField SLp(IOobject("SLp",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
    alphaF*mag(GammaEle) - (rep*ne*np + rnp*nn*np)
);
// Negative ion sourse term
volScalarField SLn(IOobject("SLn",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
    etaF*mag(GammaEle) - rnp*nn*np
);

dimensionedScalar dimAlpha("dimAlpha",alpha.dimensions(),scalar(1));
dimensionedScalar dimE("dimE",EA.dimensions(),scalar(1));
dimensionedScalar dimMu("dimMu",mue.dimensions(),scalar(1));

#include "createFvOptions.H"

// Space charge density
Info<< "Reading field rhoq\n" << endl;
volScalarField rhoq(IOobject("rhoq",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
    e*(np-ne-nn)
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
