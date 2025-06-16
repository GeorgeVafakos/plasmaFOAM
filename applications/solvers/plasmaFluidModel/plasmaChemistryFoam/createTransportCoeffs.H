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
            if (coeffType == "uniformValue")
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
            else if (coeffType == "mobilityHeIonLXCat")
            {
                mobilityHandler.set
                (
                    i,
                    new transpCoeffsHeIonLXCat
                    (
                        mesh,
                        coeffDict,
                        mobilityCoeffSpecies[i]
                    )
                );
            }
            else if (coeffType == "mobilityElectronsRelation")
            {
                mobilityHandler.set
                (
                    i,
                    new transpCoeffsElectronsRelation
                    (
                        mesh,
                        coeffDict,
                        mobilityCoeffSpecies[i]
                    )
                );
            }
            else if (coeffType == "electronMobilityBolsig")
            {
                mobilityHandler.set
                (
                    i,
                    new transportCoeffsBolsig
                    (
                        mesh,
                        coeffDict,
                        mobilityCoeffSpecies[i]
                    )
                );
            }
            else
            {
                FatalErrorInFunction
                    << "mobility coefficient type: " << coeffType
                    << " is not a valid choice. "
                    << "Options are: uniformValue, mobilityHeIonLXCat, mobilityElectronsRelation"
                    << abort(FatalError);
            }
        }
        else if (transportCoeff == "diffusion")
        {
            // Populate diffusion coefficients pointer lists
            if (coeffType == "uniformValue")
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
            else if (coeffType == "diffusionHeIonLXCat")
            {
                diffusionHandler.set
                (
                    i,
                    new transpCoeffsHeIonLXCat
                    (
                        mesh,
                        coeffDict,
                        diffusionCoeffSpecies[i]
                    )
                );
            }
            else if (coeffType == "EinsteinRelation")
            {
                diffusionHandler.set
                (
                    i,
                    new diffusionCoeffsEinsteinRelation
                    (
                        mesh,
                        coeffDict,
                        diffusionCoeffSpecies[i],
                        mobilityCoeffSpecies[i],
                        chargeNumber[i]
                    )
                );
            }
            else
            {
                FatalErrorInFunction
                    << "diffusion coefficient type: " << coeffType
                    << " is not a valid choice. "
                    << "Options are: uniformValue, diffusionHeIonLXCat, EinsteinRelation"
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

