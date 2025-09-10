// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Define transport coefficient related objects
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
Info<< "Reading species properties\n" << endl;
IOdictionary speciesProperties
(
    IOobject
    (
        "speciesProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Define pointer lists
scalarList chargeNumber(composition.species().size());
PtrList<dimensionedScalar> molarMass(composition.species().size());
PtrList<volScalarField> mobilityCoeffSpecies(composition.species().size());
PtrList<volScalarField> diffusionCoeffSpecies(composition.species().size());
PtrList<transportCoeffsBase> mobilityCoeffHandler(composition.species().size());
PtrList<transportCoeffsBase> diffusionCoeffHandler(composition.species().size());

// Loop over species
forAll(composition.species(), i)
{
    const word& specie = composition.species()[i];
    const dictionary& specieDict(speciesProperties.subDict(specie));

    // Read properties
    chargeNumber[i] = specieDict.get<int>("chargeNumber");

    molarMass.set
    (
        i,
        new dimensionedScalar
        (
            "molarMass",
            dimensionSet(1,0,0,0,-1,0,0),
            1.0e-3 * specieDict.get<scalar>("molarMass")
        )
    );

    // Info<< "Species: " << specie;

    // --- Mobility
    const dictionary& mobDict = specieDict.subDict("mobility");

    mobilityCoeffSpecies.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "mobilityCoeffSpecies_" + specie,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimensionSet(-1, 0, 2, 0, 0, 1, 0), scalar(0))
        )
    );

    autoPtr<transportCoeffsBase> model = transportCoeffsBase::NewMobility
    (
        mesh, 
        mobDict, 
        mobilityCoeffSpecies[i]
    );

    mobilityCoeffHandler.set
    (
        i, 
        model.release()
    );  

    // mobilityCoeffHandler[i].calcMobilityCoeffs(); 

    // Create diffusion coefficient
    const dictionary& diffDict = specieDict.subDict("diffusion");

    diffusionCoeffSpecies.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "diffusionCoeffSpecies_" + specie,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimViscosity, scalar(0))
        )
    );

    autoPtr<transportCoeffsBase> diffModel = transportCoeffsBase::NewDiffusion
    (
        mesh, 
        diffDict, 
        diffusionCoeffSpecies[i]
    );

    diffusionCoeffHandler.set
    (
        i, 
        diffModel.release()
    );

    // Update diffusion coefficient
    // diffusionCoeffHandler[i].calcDiffusionCoeffs();
}

// // Output info
// forAll(composition.species(), i)
// {
//     Info<< "Species " << composition.species()[i]
//         << ": Mobility max = " << gMax(mobilityCoeffSpecies[i])
//         << ", Diffusion max = " << gMax(diffusionCoeffSpecies[i]) << endl;
// }














/*
Info<< "Reading species properties\n" << endl;
IOdictionary speciesProperties
(
    IOobject
    (
        "speciesProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Define pointer lists
PtrList<volScalarField> mobilityCoeffSpecies(composition.species().size());
PtrList<volScalarField> diffusionCoeffSpecies(composition.species().size());
scalarList chargeNumber(composition.species().size());
PtrList<dimensionedScalar> molarMass(composition.species().size());
PtrList<transportCoeffsHandler> mobilityHandler(composition.species().size());
PtrList<transportCoeffsHandler> diffusionHandler(composition.species().size());

// Populate pointer lists
forAll(composition.species(), i)
{
    mobilityCoeffSpecies.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "mobilityCoeffSpecies",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("temp", dimensionSet(-1, 0, 2, 0, 0, 1, 0), scalar(1))
        )
    );

    diffusionCoeffSpecies.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "diffusionCoeffSpecies",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("temp", dimViscosity, scalar(1))
        )
    );
}

// Read speciesProperties dict
forAll(composition.species(), i)
{
    // Define species name and dict
    const word& specie = composition.species()[i];
    const dictionary& specieDict(speciesProperties.subDict(specie));

    // Get species charge number
    int chargeNumberValue(specieDict.get<int>("chargeNumber"));
    chargeNumber[i] = chargeNumberValue;

    // Get species molar mass (in SI)
    scalar molarMassValue(specieDict.get<scalar>("molarMass"));
    molarMass.set
    (
        i,
        new dimensionedScalar
        (
            "molarMass",
            dimensionSet(1,0,0,0,-1,0,0),
            1.0e-3*molarMassValue
        )
    );

    Info<< "Species: " << specie << endl;
    for (const word& transportCoeff : std::initializer_list<word>{"mobility", "diffusion"})
    {
        const dictionary& coeffDict(specieDict.subDict(transportCoeff));
        word coeffType(coeffDict.get<word>("type"));
        if (transportCoeff == "mobility")
        {
            // Populate mobility coefficients pointer lists
            if (coeffType == "transportCoeffsConstValue")
            {
                mobilityHandler.set
                (
                    i,
                    new transpCoeffsConstValue
                    (
                        mesh,
                        coeffDict,
                        mobilityCoeffSpecies[i]
                    )
                );
            }
            else if (coeffType == "transportCoeffsNone")
            {
                mobilityHandler.set
                (
                    i,
                    new transpCoeffsNone
                    (
                        mesh,
                        coeffDict,
                        mobilityCoeffSpecies[i]
                    )
                );
            }
            else if (coeffType =
            else
            {
                FatalErrorInFunction
                    << "mobility coefficient type: " << coeffType
                    << " is not a valid choice. "
                    << abort(FatalError);
            }
        }
        else if (transportCoeff == "diffusion")
        {
            // Populate diffusion coefficients pointer lists
            if (coeffType == "transportCoeffsConstValue")
            {
                diffusionHandler.set
                (
                    i,
                    new transpCoeffsConstValue
                    (
                        mesh,
                        coeffDict,
                        diffusionCoeffSpecies[i]
                    )
                );
            }
            else if (coeffType == "transportCoeffsNone")
            {
                diffusionHandler.set
                (
                    i,
                    new transpCoeffsNone
                    (
                        mesh,
                        coeffDict,
                        diffusionCoeffSpecies[i]
                    )
                );
            }
            else
            {
                FatalErrorInFunction
                    << "diffusion coefficient type: " << coeffType
                    << " is not a valid choice. "
                    << abort(FatalError);
            }
        }
    }
}

forAll(composition.species(), i)
{
    // mobilityCoeffSpecies[i] = transpCoeffsHandler[i].calcTranspCoeffs();
    Info<< "Species " << composition.species()[i] << ": Mobility = " << gMax(mobilityCoeffSpecies[i]) << " Diffusion = " << gMax(diffusionCoeffSpecies[i]) << endl;
}
*/
