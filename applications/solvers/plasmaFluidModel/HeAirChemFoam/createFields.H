#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;
autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
psiReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

const word inertSpecie(thermo.get<word>("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species() << exit(FatalIOError);
}

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField& p = thermo.p();

#include "compressibleCreatePhi.H" // NOTE TO MYSELF: This is important for calculating Ji

pressureControl pressureControl(p, rho, pimple.dict(), false);

mesh.setFluxRequired(p.name());

Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating reaction model\n" << endl;
autoPtr<CombustionModel<psiReactionThermo>> reaction
(
    CombustionModel<psiReactionThermo>::New(thermo, turbulence())
);

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

volScalarField Qdot
(
    IOobject
    (
        "Qdot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
);

#include "createDpdt.H"

#include "createK.H"

#include "createMRF.H"
#include "createFvOptions.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Import mass fractions from the transport model
// const volScalarField& Y = thermo.Y();
// 
// #include "chemistryModel.H"
// autoPtr<chemistryModel> chemistry
// (
//     chemistryModel::New
//     (
//         thermo
//     )
// );
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Create electrical fields
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "createSolidFields.H"

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
    dimensionSet(1, -1, -2, 0, 0, 0, 0),
    physicalProperties
);

// Gas Temperature
const dimensionedScalar Tgas
(
    "Tgas",
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

// Info<< "Calculating field YiFlux\n" << endl;
// surfaceScalarField YiFlux
// (
//     IOobject
//     (
//         "YiFlux",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     linearInterpolate(mobilityCoeffSpecies[1])*mesh.magSf()*fvc::snGrad(volt)
// ); 


dimensionedScalar dimRho("dimAlpha",dimDensity,scalar(1));
dimensionedScalar dimmass("dimnass",dimMass,scalar(1));



#include "createTransportCoeffs.H"

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
    // Define name for species
    const word& nName("n"+composition.species()[i]);

    IOobject header
    (
        nName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );

    // Info<< "Species " << composition.species()[i] << " TRUE OR FALSE ========== "<< header.typeHeaderOk<volScalarField>(true) << endl;

    // This is for NEqn
    // ------------------
    
    if (header.typeHeaderOk<volScalarField>(true))
    {
        n.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    nName,
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
                    nName, 
                    runTime.timeName(), 
                    mesh, 
                    IOobject::NO_READ, 
                    IOobject::AUTO_WRITE
                ),
                tnDefault()
            )
        );
    }




    // This is for YEqn
    // ------------------

    // n.set
    // (
    //     i,
    //     new volScalarField
    //     (
    //         IOobject
    //         (
    //             nName,
    //             runTime.timeName(),
    //             mesh,
    //             IOobject::NO_READ,
    //             IOobject::AUTO_WRITE
    //         ),
    //         mesh,
    //         dimensionedScalar("temp", dimensionSet(0, -3, 0, 0, 0, 0, 0), scalar(0))
    //     )
    // );
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





// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// He rates
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
const wordList reactionNames({"R1","R2","R3","R4","R5","R6","R7","R8","R9"});
PtrList<dimensionedScalar> rate(reactionNames.size());
dimensionedScalar dimRate("dimRate",dimensionSet(0,3,-1,0,0,0,0),scalar(1));
dimensionedScalar dimRate2("dimRate2",dimensionSet(0,6,-1,0,0,0,0),scalar(1));

// Populate pointer lists
forAll(reactionNames, i)
{
    rate.set
    (
        i,
        new dimensionedScalar
        (
            "rate",
            dimless,
            scalar(1)
        )
    );
}


PtrList<volScalarField> Ri_;
Ri_.setSize(n.size());

forAll(n, i)
{
	Ri_.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"R_" + n[i].name(),
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionedScalar("Ri", dimensionSet(0, -3, -1, 0, 0, 0, 0), 0)
			//reaction->R(Y[i]) & Y[i]
		)
	);
	// if I write Ri_[i] = reaction->R(Y[i]); the compiler complains
	// Ri_[i] = reaction->R(Y[i]) & Y[i]; 
	// Ri_[i].write(); 
}


Info<< "Reading field HeF\n" << endl;
volScalarField HeF
(
    IOobject
    (
        "HeF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("temp", dimensionSet(0,0,0,0,0,0,0), scalar(1.0))
);