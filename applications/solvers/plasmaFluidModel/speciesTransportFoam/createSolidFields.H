PtrList<IOdictionary> dielectricPropertiesSolid(solidRegions.size());
PtrList<dimensionedScalar> erSolid(solidRegions.size());
PtrList<volScalarField> epsilonSolid(solidRegions.size());
PtrList<volScalarField> voltExtAmpSolid(solidRegions.size());
PtrList<volScalarField> voltExtSolid(solidRegions.size());
PtrList<volScalarField> voltIndSolid(solidRegions.size());
PtrList<volScalarField> voltSolid(solidRegions.size());
PtrList<volVectorField> EExtAmpSolid(solidRegions.size());
PtrList<volVectorField> EExtSolid(solidRegions.size());
PtrList<volVectorField> EIndSolid(solidRegions.size());
PtrList<volVectorField> ESolid(solidRegions.size());

// Populate solid field pointer lists
forAll(solidRegions, i)
{
    Info<< "    Adding to dielectricPropertiesSolid\n" << endl;
    dielectricPropertiesSolid.set
    (
        i,
        new IOdictionary
        (
            IOobject
            (
                "dielectricProperties",
                runTime.constant(),
                solidRegions[i],
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        )
    );

    Info<< "    Adding to erSolid\n" << endl;
    erSolid.set
    (
        i,
        new dimensionedScalar
        (
            "er",
            dimless,
            dielectricPropertiesSolid[i]
        )
    );

    Info<< "    Adding to epsilonSolid\n" << endl;
    epsilonSolid.set
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
            erSolid[i]
        )
    );
    
    Info<< "    Adding to voltExtAmpSolid\n" << endl;
    voltExtAmpSolid.set
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

    Info<< "    Adding to voltExtSolid\n" << endl;
    voltExtSolid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "voltExt",
                runTime.timeName(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i],
            dimensionedScalar("temp", dimensionSet(1,2,-3,0,0,-1,0), scalar(0))
        )
    );

    Info<< "    Adding to voltIndSolid\n" << endl;
    voltIndSolid.set
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

    Info<< "    Adding to voltSolid\n" << endl;
    voltSolid.set
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

    Info<< "    Adding to EExtAmpSolid\n" << endl;
    EExtAmpSolid.set
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
            -fvc::grad(voltExtAmpSolid[i])
        )
    );

    Info<< "    Adding to EExtSolid\n" << endl;
    EExtSolid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "EExt",
                runTime.timeName(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            -fvc::grad(voltExtSolid[i])
        )
    );

    Info<< "    Adding to EIndSolid\n" << endl;
    EIndSolid.set
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
            -fvc::grad(voltIndSolid[i])
        )
    );

    Info<< "    Adding to ESolid\n" << endl;
    ESolid.set
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
            -fvc::grad(voltSolid[i])
        )
    );

}
