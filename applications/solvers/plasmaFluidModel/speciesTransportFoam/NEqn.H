
{
    reaction->correct();
    Qdot = reaction->Qdot();
    volScalarField nt(0.0*n[0]);

    forAll(n, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            // Info << "Species: " << i << " Inert: " << inertSpecie << endl;
            // Y[i] =  n[i]*molarMass[i]/(rhoGas*constant::physicoChemical::NA);
            // volScalarField& Yi = Y[i];
            volScalarField& ni = n[i];
            surfaceScalarField& YiFlux = YFlux[i];
            volScalarField& Di = diffusionCoeffSpecies[i];

            YiFlux = linearInterpolate(chargeNumber[i]*mobilityCoeffSpecies[i])*mesh.magSf()*fvc::snGrad(volt);

            
            fvScalarMatrix niEqn
            (
                // fvm::ddt(rho, Yi)
                fvm::ddt(ni)
                //   + mvConvection->fvmDiv(phi, Yi)
                + fvm::div(YiFlux, ni)
                //   - fvm::laplacian(turbulence->muEff(), Yi)
                - fvm::laplacian(Di, ni)
                //  ==
                //     reaction->R(Yi)
                //   + fvOptions(rho, Yi)
            );
            
            
            
            niEqn.relax();
            
            fvOptions.constrain(niEqn);
            
            niEqn.solve("ni");
            
            fvOptions.correct(ni);
            
            ni.clamp_min(0);
            // Info<< "TEEEEEEEEEEEEEEEST" << endl;
            nt += ni;
            
            // volScalarField& ni = n[i];
            // n[i] = rhoGas*constant::physicoChemical::NA*Yi/(molarMass[i]);
        }
    }

    n[inertIndex] = N - nt;
    n[inertIndex].clamp_min(0);

}
