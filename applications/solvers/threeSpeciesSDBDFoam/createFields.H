// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Create electrical fields
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "createSolidFields.H"

Info<< "Reading voltageProperties\n" << endl;
IOdictionary voltageProperties
(
    IOobject
    (
        "voltageProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
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

// Gas Pressure
const dimensionedScalar press
(
    "press",
    // dimensionSet(1, -1, -2, 0, 0, 0, 0),
    dimless,
    physicalProperties
);

// Gas Temperature
const dimensionedScalar Tgas
(
    "Tgas",
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

// Gas number density
volScalarField N
(
    IOobject
    (
        "N",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    press/(constant::physicoChemical::k * Tgas)
);

// Ambient pressure field
volScalarField P
(
    IOobject
    (
        "P",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    press
);

volScalarField ones
(
    IOobject
    (
        "ones",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("temp",dimensionSet(0,0,0,0,0,0,0),scalar(1))
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
    dimensionedScalar("temp",dimless,scalar(1))
);

// Region convergence variables
int regionLoopCounter=0;
std::vector<int> voltEqnIter = {1};

// Read controlDict for maximum region loops
int maxRegionLoops =
    runTime.controlDict().lookupOrDefault("maxRegionLoops", 10000);

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
    mesh,
    dimensionedScalar("temp", dimensionSet(1,2,-3,0,0,-1,0), scalar(0))
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

Info<< "Calculating field EExtAmp\n" << endl;
volVectorField EExtAmp
(
    IOobject
    (
        "EExtAmp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(voltExtAmp)
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

Info<< "Calculating field E_N\n" << endl;
volScalarField E_N
(
    IOobject
    (
        "E_N",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    1.0e21*mag(E/N)
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Create species fields
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //




// // Define species flux
// forAll(composition.species(), i)
// {
//     YFlux.set
//     (
//         i,
//         new surfaceScalarField
//         (
//             IOobject
//             (
//                 "YFlux",
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE
//             ),
//             linearInterpolate(mobilityCoeffSpecies[i])*mesh.magSf()*fvc::snGrad(volt)
//         )
//     );
// }


Info<< "Calculating field surfaceCharge\n" << endl;
volScalarField surfaceCharge
(
    IOobject
    (
        "surfaceCharge",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("temp", dimensionSet(0,-3,1,0,0,1,0), scalar(1.0e5))
    // linearInterpolate(rhoq)//*mesh.magSf()
);

PtrList<volScalarField> surfaceChargeSolid(solidRegions.size());
forAll(solidRegions, i)
{
    Info<< "    Adding to surfaceChargeSolid\n" << endl;
    surfaceChargeSolid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "surfaceCharge",
                runTime.timeName(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            solidRegions[i],
            dimensionedScalar("temp", dimensionSet(0,-3,1,0,0,1,0), scalar(1.0e5))
        )
    );
}

// Electron Diffusion coefficient
dimensionedScalar De("De", dimensionSet(0, 2, -1, 0, 0, 0, 0), physicalProperties);

// Ion Diffusion coefficient
dimensionedScalar Di("Di", dimensionSet(0, 2, -1, 0, 0, 0, 0), physicalProperties);

// Electron mobility scalar
dimensionedScalar mue("mue", dimensionSet(-1, 0, 2, 0, 0, 1, 0), physicalProperties);
volScalarField mueF
(
    IOobject
    (
        "mueF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    mue * dimensionedScalar("temp",dimless,scalar(1))
);

// Positive Ion mobility
dimensionedScalar mup("mup", dimensionSet(-1, 0, 2, 0, 0, 1, 0), physicalProperties);
volScalarField mupF
(
    IOobject
    (
        "mupF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    mup * dimensionedScalar("temp",dimless,scalar(1))
);

// Negative Ion mobility
dimensionedScalar mun("mun", dimensionSet(-1, 0, 2, 0, 0, 1, 0), physicalProperties);
volScalarField munF
(
    IOobject
    (
        "munF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    mun * dimensionedScalar("temp",dimless,scalar(1))
);

// Recombination coefficients
dimensionedScalar r("r", dimensionSet(-1, 0, 2, 0, 0, 1, 0), physicalProperties);
dimensionedScalar rnp("rnp", dimensionSet(0, 3, -1, 0, 0, 0, 0), physicalProperties);
dimensionedScalar rep("rep", dimensionSet(0, 3, -1, 0, 0, 0, 0), physicalProperties);

// Ionization coefficient
Info<< "Calculating field alpha\n" << endl;
dimensionedScalar alpha("alpha", dimensionSet(0, -1, 0, 0, 0, 0, 0), physicalProperties);
volScalarField alphaF
(
    IOobject
    (
        "alphaF", 
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    alpha * dimensionedScalar("temp",dimensionSet(0,0,0,0,0,0,0),scalar(1))
);

// Recombinations rate
Info<< "Calculating field eta\n" << endl;
volScalarField etaF
(
    IOobject
    (
        "etaF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    alpha * dimensionedScalar("temp",dimensionSet(0,0,0,0,0,0,0),scalar(1))
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

// Ion number density
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
    -mue*mesh.magSf()*fvc::snGrad(volt)
);

//Positive Ion flux
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
    -mup*mesh.magSf()*fvc::snGrad(volt)
);

//Negative Ion flux
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
    -mun*mesh.magSf()*fvc::snGrad(volt)
);

// Electron diffusion coefficient
Info<< "Calculating field DeF" << endl;
volScalarField DeF
(
    IOobject
    (
        "DeF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("temp",dimensionSet(0,2,-1,0,0,0,0),scalar(0))
);

// Positive Ion diffusion coefficient
Info<< "Calculating field DpF" << endl;
volScalarField DpF
(
    IOobject
    (
        "DpF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("temp",dimensionSet(0,2,-1,0,0,0,0),scalar(0))
);

// Negative Ion diffusion coefficient
Info<< "Calculating field DnF\n" << endl;
volScalarField DnF
(
    IOobject
    (
        "DnF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("temp",dimensionSet(0,2,-1,0,0,0,0),scalar(0))
);

// Electron flux vector
Info<< "Calculating field GammaEle" << endl;
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
    -mueF*ne*E - DeF*fvc::grad(ne)
);

// Info<< "TEEEEEEEEEEEEEEEEEEEEEEEEEEEEEST" << endl;

// Positive ion flux vector
Info<< "Calculating field GammaPos" << endl;
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
    mupF*np*E - DpF*fvc::grad(np)
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
    -munF*nn*E - DnF*fvc::grad(nn)
);

// Electron sourse term
Info<< "Calculating field SLe" << endl;
volScalarField SLe
(
    IOobject
    (
        "SLe",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (alphaF-etaF)*mag(GammaEle) - rep*ne*np // *****************************************************************
);

// Positive ion sourse term
Info<< "Calculating field SLp" << endl;
volScalarField SLp
(
    IOobject
    (
        "SLp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alphaF*mag(GammaEle) - (rep*ne*np + rnp*nn*np)
);

// Negative ion sourse term
Info<< "Calculating field SLn\n" << endl;
volScalarField SLn
(
    IOobject
    (
        "SLn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    etaF*mag(GammaEle) - rnp*nn*np
);

// Dimension units
dimensionedScalar dimAlpha("dimAlpha",alpha.dimensions(),scalar(1));
dimensionedScalar dimE("dimE",E.dimensions(),scalar(1));
dimensionedScalar dimMu("dimMu",mue.dimensions(),scalar(1));

// Space charge density
Info<< "Calculating field rhoq\n" << endl;
volScalarField rhoq
(
    IOobject
    (
        "rhoq",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    constant::electromagnetic::e*(np-ne-nn)
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