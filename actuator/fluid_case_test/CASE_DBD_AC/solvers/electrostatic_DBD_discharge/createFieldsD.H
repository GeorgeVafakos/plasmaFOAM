Info<< "Reading field phiElDiel\n" << endl;
volScalarField phiElDiel
(
    IOobject
    (
        "phiElDiel",
        runTime.timeName(),
        meshD,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshD
);

volScalarField k_Diel
(
    IOobject
    (
        "k_Diel",
        runTime.constant(),
        meshD,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    meshD
);

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
    -fvc::grad(phiElDiel)
); 


// Variables for check of the convergence
int convphiElDiel=1;
