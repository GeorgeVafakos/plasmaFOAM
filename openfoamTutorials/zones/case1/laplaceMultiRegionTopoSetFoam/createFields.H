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
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField kappaAir
(
    IOobject
    (
        "kappaAir",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    kAir
);

Info<< "Calculating field qAir\n" << endl;
volVectorField qAir
(
    IOobject
    (
        "qAir",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -kAir*fvc::grad(Tair)
);
