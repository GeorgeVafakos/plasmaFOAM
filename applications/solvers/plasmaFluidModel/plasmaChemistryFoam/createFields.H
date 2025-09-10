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

// Gas pressure value
const dimensionedScalar pGasValue
(
    "pGasValue",
    dimensionSet(1, -1, -2, 0, 0, 0, 0),
    physicalProperties
);

// Gas temperature value
const dimensionedScalar TGasValue
(
    "TGasValue",
    dimensionSet(0, 0, 0, 1, 0, 0, 0),
    physicalProperties
);

// Electron temperature value
const dimensionedScalar TEleValue
(
    "TEleValue",
    dimensionSet(0, 0, 0, 1, 0, 0, 0),
    physicalProperties
);

// Gas density value
const dimensionedScalar rhoGasValue
(
    "rhoGasValue",
    dimDensity,
    physicalProperties
);

// Gas number density field
volScalarField NGas
(
    IOobject
    (
        "NGas",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    pGasValue/(constant::physicoChemical::k * TGasValue)
);

// Ambient pressure field
volScalarField pGas
(
    IOobject
    (
        "pGas",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    pGasValue
);

// Electron temperature field
volScalarField TEle
(
    IOobject
    (
        "TEle",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    TEleValue
);

// Gas temperature field
volScalarField TGas
(
    IOobject
    (
        "TGas",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    TGasValue
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
    1.0e21*mag(E/NGas)
);

dimensionedScalar maxVal("maxVal", EInd.dimensions(), 5.0e6);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Create species
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Instantiate plasma chemistry model
plasmaChemistryModel composition(mesh);
int inertIndex = composition.inertIndex();

PtrList<volScalarField> N(composition.species().size());

// Create a temporary fallback field
tmp<volScalarField> tNdefault;

// Define name for default field
word NdefaultName("Ndefault");

// Create an IOobject to check if "nDefault" exists in the case
IOobject NdefaultIO
(
    NdefaultName,
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::NO_WRITE
);

// Check if nDefault exists
if (NdefaultIO.typeHeaderOk<volScalarField>(true))
{    
    tNdefault = new volScalarField
    (
        NdefaultIO, 
        mesh
    );
}
else
{
    // If not found, create a default empty field
    tNdefault = new volScalarField
    (
        IOobject
        (
            NdefaultName, 
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
    IOobject header
    (
        composition.species()[i],
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    
    if (header.typeHeaderOk<volScalarField>(true))
    {
        N.set
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
        N.set
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
                tNdefault()
            )
        );
    }
}

// Initalise all reaction rate coefficients
composition.defineReactionRateCoeffs();

// Define and initialize species transport coefficients
#include "createTransportCoeffs.H"

PtrList<surfaceScalarField> NFlux(composition.species().size());

// Define species flux
forAll(composition.species(), i)
{
    NFlux.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "NFlux",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            linearInterpolate(mobilityCoeffSpecies[i])*mesh.magSf()*fvc::snGrad(volt)
        )
    );
}

// Read mass fractions
Info<< "Reading field YHe\n" << endl;
volScalarField YHe
(
    IOobject
    (
        "YHe",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
N[composition.species().find("He")] = YHe*NGas;

Info<< "Reading field YN2\n" << endl;
volScalarField YN2
(
    IOobject
    (
        "YN2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
N[composition.species().find("N2")] = YN2*NGas;

Info<< "Reading field YO2\n" << endl;
volScalarField YO2
(
    IOobject
    (
        "YO2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
N[composition.species().find("O2")] = YO2*NGas;

composition.updateReactionRateCoeffs();

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
