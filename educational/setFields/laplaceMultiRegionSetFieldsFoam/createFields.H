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

dimensionedScalar kAir
(
    physicalProperties.lookup("kAir")
);

Info<< "Reading field Tair\n" << endl;
volScalarField Tair
(
    IOobject
    (
        "Tair",
        runTime.timeName(),
        meshA,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshA
);

volScalarField kappaAir
(
    IOobject
    (
        "kappaAir",
        runTime.timeName(),
        meshA,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    meshA,
    kAir
);

Info<< "Calculating field qAir\n" << endl;
volVectorField qAir
(
    IOobject
    (
        "qAir",
        runTime.timeName(),
        meshA,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -kAir*fvc::grad(Tair)
);

Info<< "Calculating field S\n" << endl;
volScalarField S
(
    IOobject
    (
        "S",
        runTime.timeName(),
        meshA,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshA
);
