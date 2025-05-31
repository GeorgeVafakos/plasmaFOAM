Foam::fvMesh meshD
(
    IOobject
    (
        "dielectric",
        runTime.timeName(),
        runTime,
        IOobject::MUST_READ
    )
);
