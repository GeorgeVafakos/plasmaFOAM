// ************************************************************************* //
//
// Create fields in the solid regions
//
// ************************************************************************* //


dimensionedScalar kSolid
(
    physicalProperties.lookup("kSolid")
);

Info<< "Reading field Tsolid\n" << endl;
volScalarField Tsolid
(
    IOobject
    (
        "Tsolid",
        runTime.timeName(),
        meshS,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshS
);

volScalarField kappaSolid
(
    IOobject
    (
        "kappaSolid",
        runTime.timeName(),
        meshS,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    meshS,
    kSolid
);

Info<< "Calculating field qSolid\n" << endl;
volVectorField qSolid
(
    IOobject
    (
        "qSolid",
        runTime.timeName(),
        meshS,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -kSolid*fvc::grad(Tsolid)
);

