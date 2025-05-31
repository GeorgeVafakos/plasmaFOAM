PtrList<IOdictionary> physicalPropertiesDiel(solidRegions.size());
PtrList<dimensionedScalar> erDiel(solidRegions.size());
PtrList<volScalarField> epsilonDiel(solidRegions.size());
PtrList<volScalarField> voltIndDiel(solidRegions.size());
PtrList<volVectorField> EIndDiel(solidRegions.size());

// Populate solid field pointer lists
forAll(solidRegions, i)
{
    Info<< "    Adding to physicalPropertiesDiel\n" << endl;
    physicalPropertiesDiel.set
    (
        i,
        new IOdictionary
        (
            IOobject
            (
                "physicalProperties",
                runTime.constant(),
                solidRegions[i],
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        )
    );

    Info<< "    Adding to erDiel\n" << endl;
    erDiel.set
    (
        i,
        new dimensionedScalar
        (
            "er",
            dimless,
            physicalPropertiesDiel[i]
        )
    );

    Info<< "    Adding to epsilonDiel\n" << endl;
    epsilonDiel.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "epsilon",
                runTime.timeName(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            solidRegions[i],
            erDiel[i]
        )
    );

    Info<< "    Adding to voltIndDiel\n" << endl;
    voltIndDiel.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "voltInd",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
    );

    Info<< "    Adding to EIndDiel\n" << endl;
    EIndDiel.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "EInd",
                runTime.timeName(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            -fvc::grad(voltIndDiel[i])
        )
    );
}

