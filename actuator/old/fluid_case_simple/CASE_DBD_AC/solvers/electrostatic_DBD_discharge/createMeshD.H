Foam::fvMesh meshD
(
    IOobject
    (
        "Dielectric",
        runTime.timeName(),
        runTime,
        IOobject::MUST_READ
    )
);
