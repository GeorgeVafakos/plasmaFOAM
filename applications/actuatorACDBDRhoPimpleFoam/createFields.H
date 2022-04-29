#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;

autoPtr<fluidThermo> pThermo
(
    fluidThermo::New(mesh)
);
fluidThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField& p = thermo.p();

// // Set the forceSwitch
// scalar forceSwitch = 0;
// scalar Tp = runTime.controlDict().lookupOrDefault("Period", 1.0e-3);
// scalar startForcePerc = runTime.controlDict().lookupOrDefault("startForcePerc", 0.0);
// scalar stopForcePerc = runTime.controlDict().lookupOrDefault("stopForcePerc", 0.25);
// scalar numPeriods = floor(runTime.value()/Tp);

// Read times where the force acts
// Info<< "Reading forceTimesDict\n" << endl;
// IOdictionary forceTimesDict
// (
//     IOobject
//     (
//         "forceTimesDict",
//         runTime.constant(),
//         mesh,
//         IOobject::MUST_READ_IF_MODIFIED,
//         IOobject::NO_WRITE
//     )
// );

// scalarList forceTimes
// (
//     forceTimesDict.lookup("forceTimes")
// );

// int iterCount = 0;
// int numSteps = forceTimes.size()-1;
// scalarList timeSteps[NumTimes] = {2,3,4,3};

// for (int i = 0; i < NumTimes; ++i)
// {
//     timeSteps[i] = timeSteps[i+1]-timeSteps[i]
// }


// Info<< "Teeeeeeeeeeeeeeeeeeeest = " << timeSteps.size() << nl << endl;


volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
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

#include "compressibleCreatePhi.H"

pressureControl pressureControl(p, rho, pimple.dict(), false);

mesh.setFluxRequired(p.name());

Info<< "Creating turbulence model\n" << endl;
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

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

#include "createMRF.H"
#include "createFvOptions.H"



// Info<< "Calculating field Fc\n" << endl;
// volVectorField Fc
// (
//     IOobject
//     (
//         "Fc",
//         runTime.timeName(),
//         mesh,
//         IOobject::MUST_READ,
//         IOobject::NO_WRITE
//     ),
//     mesh
// );
