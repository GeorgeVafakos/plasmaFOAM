// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Create general fields
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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

// Gas pressure
const dimensionedScalar pGas
(
    "pGas",
    dimensionSet(1, -1, -2, 0, 0, 0, 0),
    physicalProperties
);

// Gas temperature
const dimensionedScalar TGas
(
    "TGas",
    dimensionSet(0, 0, 0, 1, 0, 0, 0),
    physicalProperties
);

// Electron temperature
const dimensionedScalar TEle
(
    "TEle",
    dimensionSet(0, 0, 0, 1, 0, 0, 0),
    physicalProperties
);

// Gas Density
const dimensionedScalar rhoGas
(
    "rhoGas",
    dimDensity,
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
    pGas/(constant::physicoChemical::k * TGas)
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
    pGas
);

// Electron temperature
volScalarField Te
(
    IOobject
    (
        "Te",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    TEle
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Create electrical fields
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "createSolidFields.H"

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

Info<< "Calculating field EN\n" << endl;
volScalarField EN
(
    IOobject
    (
        "EN",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    1.0e21*mag(E/N)
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Create species
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Instantiate plasma chemistry model
plasmaChemistryModel composition(mesh);
composition.printStoichiometricMatrices();



// composition.setRateCalculator
// (
    //     autoPtr<reactionRateCoeffsBase>(new Foam::constantRateCoeff())
    // );
    
    
    // composition.setRateCalculator
    // (
        //     autoPtr<reactionRateCoeffsBase>(new Foam::ArrheniusRateCoeff())
        // );
        
        
composition.updateReactionRateCoeffs();
Info<< "----------------------------- Reading field mobilityCoeffSpecies\n" << endl;



int inertIndex = composition.inertIndex();

#include "createTransportCoeffs.H"

Info<< "Reading field rhoq\n" << endl;
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
    mesh,
    dimensionedScalar("temp", dimensionSet(0,-3,1,0,0,1,0), scalar(0))
);

PtrList<surfaceScalarField> YFlux(composition.species().size());
PtrList<volScalarField> n(composition.species().size());

// Create a temporary fallback field
tmp<volScalarField> tnDefault;

// Define name for default field
word nDefaultName("nDefault");

// Create an IOobject to check if "nDefault" exists in the case
IOobject nDefaultIO
(
    nDefaultName,
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::NO_WRITE
);

// Check if nDefault exists
if (nDefaultIO.typeHeaderOk<volScalarField>(true))
{    
    tnDefault = new volScalarField
    (
        nDefaultIO, 
        mesh
    );
}
else
{
    // If not found, create a default empty field
    tnDefault = new volScalarField
    (
        IOobject
        (
            nDefaultName, 
            runTime.timeName(), 
            mesh, 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("temp", dimensionSet(0, -3, 0, 0, 0, 0, 0), scalar(0))
    );
}

forAll(composition.species(), i)
{
    // // Define name for species
    // const word& nName("n"+composition.species()[i]);

    IOobject header
    (
        composition.species()[i],
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    
    if (header.typeHeaderOk<volScalarField>(true))
    {
        n.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    composition.species()[i],
                    runTime.timeName(),
                    mesh, 
                    IOobject::MUST_READ, 
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }
    else
    {
        n.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    composition.species()[i], 
                    runTime.timeName(), 
                    mesh, 
                    IOobject::NO_READ, 
                    IOobject::AUTO_WRITE
                ),
                tnDefault()
            )
        );
    }
}

// Define species flux
forAll(composition.species(), i)
{
    YFlux.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "YFlux",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            linearInterpolate(mobilityCoeffSpecies[i])*mesh.magSf()*fvc::snGrad(volt)
        )
    );
}
