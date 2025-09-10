// solverPerformance::debug = 0;
{
    volScalarField Nt(0.0*N[0]);
    
    // Update only the non constant reaction rate coefficients at selected time steps
    if (customUpdatePlasmaCoeffs && runTime.timeIndex() % plasmaCoeffsUpdateInterval == 0)
    {
        composition.updateReactionRateCoeffs();
    }

    forAll(N, i)
    {
        
        if (i != inertIndex)
        {
            if (runTime.timeIndex() % plasmaCoeffsUpdateInterval == 0)
            {
                // Update mobility and diffusion coefficients
                if (!mobilityCoeffHandler[i].isConstant())
                {
                    mobilityCoeffHandler[i].calcMobilityCoeffs();
                }

                if (!diffusionCoeffHandler[i].isConstant())
                {
                    diffusionCoeffHandler[i].calcDiffusionCoeffs();
                }
            }

            volScalarField& Ni = N[i];
            surfaceScalarField& NiFlux = NFlux[i];
            volScalarField& Di = diffusionCoeffSpecies[i];

            // Info<< "Species " << i << " " << composition.species()[i] << ": Mobility = " << gMax(mobilityCoeffSpecies[i]) << " Diffusion = " << gMax(diffusionCoeffSpecies[i]) << endl;
            NiFlux = linearInterpolate(chargeNumber[i]*mobilityCoeffSpecies[i])*mesh.magSf()*fvc::snGrad(volt);


            fvScalarMatrix NiEqn
            (
                fvm::ddt(Ni)
                + fvm::div(NiFlux, Ni)
                - fvm::laplacian(Di, Ni)
                ==
                composition.R(i,N)
            );   
            
            NiEqn.relax();
            
            NiEqn.solve("Ni");
            
            Ni.clamp_min(0);
            Nt += Ni;
        }
    }

    N[inertIndex] = NGas - Nt;
    N[inertIndex].clamp_min(0);

}
