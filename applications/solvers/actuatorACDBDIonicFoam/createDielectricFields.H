PtrList<IOdictionary> physicalPropertiesDiel(solidRegions.size());
PtrList<dimensionedScalar> erDiel(solidRegions.size());
PtrList<volScalarField> voltExtAmpDiel(solidRegions.size());
PtrList<volScalarField> voltIndDiel(solidRegions.size());
PtrList<volScalarField> voltIndDischargeDiel(solidRegions.size());
PtrList<volScalarField> voltDiel(solidRegions.size());
PtrList<volVectorField> EExtAmpDiel(solidRegions.size());
PtrList<volVectorField> EIndDiel(solidRegions.size());
PtrList<volVectorField> EDiel(solidRegions.size());
PtrList<volScalarField> epsilonDiel(solidRegions.size());

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

    Info<< "    Adding to voltExtAmpDiel\n" << endl;
    voltExtAmpDiel.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "voltExtAmp",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
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

    Info<< "    Adding to voltIndDischargeDiel\n" << endl;
    voltIndDischargeDiel.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "voltIndDischargeDiel",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
    );

    Info<< "    Adding to voltDiel\n" << endl;
    voltDiel.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "volt",
                runTime.timeName(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i],
            dimensionedScalar("temp", dimensionSet(1,2,-3,0,0,-1,0), scalar(0))
        )
    );

    Info<< "    Adding to EExtAmpDiel\n" << endl;
    EExtAmpDiel.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "EExtAmp",
                runTime.timeName(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            -fvc::grad(voltExtAmpDiel[i])
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

    Info<< "    Adding to EDiel\n" << endl;
    EDiel.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "E",
                runTime.timeName(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            -fvc::grad(voltDiel[i])
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
}

