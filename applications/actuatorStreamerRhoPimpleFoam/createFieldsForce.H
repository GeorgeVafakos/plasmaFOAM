// ************************************************************************* //
// Create Force
// ************************************************************************* //
Info<< "Reading forceDict\n" << endl;
IOdictionary forceDict
(
    IOobject
    (
        "forceDict",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

scalarList forceTimes
(
    forceDict.lookup("forceTimes")
);

int iterCount = 0;
int numSteps = forceTimes.size()-1;

wordList forcesNames
(
    forceDict.lookup("forceNames")
);

PtrList <volVectorField> Fc(forcesNames.size());

forAll(forcesNames, s)
{
    Fc.set
    (
        s,
        new volVectorField
        (
            IOobject
            (
                forcesNames[s],
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
        mesh
        )
    );
}