Foam::fvMesh meshI
(
    IOobject
    (
        "insulator",
        runTime.timeName(),
        runTime,
        IOobject::MUST_READ
    )
);
