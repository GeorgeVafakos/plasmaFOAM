
{
    volScalarField nt(0.0*n[0]);

    forAll(n, i)
    {
        
        if (i != inertIndex)
        {
            mobilityHandler[i].calcMobilityCoeffs();
            diffusionHandler[i].calcDiffusionCoeffs();

            volScalarField& ni = n[i];
            surfaceScalarField& YiFlux = YFlux[i];
            volScalarField& Di = diffusionCoeffSpecies[i];

            Info<< "Species " << i << " " << composition.species()[i] << ": Mobility = " << gMax(mobilityCoeffSpecies[i]) << " Diffusion = " << gMax(diffusionCoeffSpecies[i]) << endl;
            YiFlux = linearInterpolate(chargeNumber[i]*mobilityCoeffSpecies[i])*mesh.magSf()*fvc::snGrad(volt);


            fvScalarMatrix niEqn
            (
                fvm::ddt(ni)
                + fvm::div(YiFlux, ni)
                - fvm::laplacian(Di, ni)
                ==
                composition.R(i,n)
            );   
            
            niEqn.relax();
            
            niEqn.solve("ni");
            
            ni.clamp_min(0);
            nt += ni;
        }
    }

    n[inertIndex] = N - nt;
    n[inertIndex].clamp_min(0);

}
