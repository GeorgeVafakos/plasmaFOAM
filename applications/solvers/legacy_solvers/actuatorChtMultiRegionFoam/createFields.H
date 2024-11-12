#include "createFluidFields.H"
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

dimensionedScalar kappaA
(
    physicalProperties.lookup("kappaA")
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
    kappaA
);
