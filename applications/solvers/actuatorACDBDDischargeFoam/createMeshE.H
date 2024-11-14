Foam::fvMesh meshE
(
    IOobject
    (
        "encapsulator",
        runTime.timeName(),
        runTime,
        IOobject::MUST_READ
    )
);
