forAll(composition.species(), i)
{
    mobilityHandler[i].calcMobilityCoeffs();
    diffusionHandler[i].calcDiffusionCoeffs();

    Info<< "Species " << i << " " << composition.species()[i] << ": Mobility = " << gMax(mobilityCoeffSpecies[i]) << " Diffusion = " << gMax(diffusionCoeffSpecies[i]) << endl;
    
}

