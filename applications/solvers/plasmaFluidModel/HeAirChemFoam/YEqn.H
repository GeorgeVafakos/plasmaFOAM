tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

{
    reaction->correct();
    Qdot = reaction->Qdot();
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            // Info << "Species: " << i << " Inert: " << inertSpecie << endl;
            // Y[i] =  n[i]*molarMass[i]/(rhoGas*constant::physicoChemical::NA);
            volScalarField& Yi = Y[i];
            // volScalarField& Yi = n[i];
            surfaceScalarField& YiFlux = YFlux[i];
            volScalarField& Di = diffusionCoeffSpecies[i];

            YiFlux = linearInterpolate(chargeNumber[i]*mobilityCoeffSpecies[i])*mesh.magSf()*fvc::snGrad(volt);

            // bool flag;

            #include "reactionTerm.H"

            Info<< "Species: "<< i << "if He then: "<<  (i==composition.species().find("e")) << endl;

            Info<< "Reaction species count: " << reaction->size() << endl;
            fvScalarMatrix YiEqn
            (
                // fvm::ddt(rho, Yi)
                fvm::ddt(Yi)
                //   + mvConvection->fvmDiv(phi, Yi)
                + fvm::div(YiFlux, Yi)
                //   - fvm::laplacian(turbulence->muEff(), Yi)
                - fvm::laplacian(Di, Yi)
                //  ==
                //     (1.0/dimRho)*reaction->R(Yi)
                    // reaction->R(Yi)
                //   + fvOptions(rho, Yi)
            );

            // Ri_[i] = reaction->R(Y[i]) & Y[i]; 
            
            YiEqn.relax();
            
            fvOptions.constrain(YiEqn);
            
            YiEqn.solve("Yi");
            
            fvOptions.correct(Yi);
            
            Yi.clamp_min(0);
            // Info<< "TEEEEEEEEEEEEEEEST" << endl;
            Yt += Yi;
            
            // volScalarField& ni = n[i];
            // n[i] = rhoGas*constant::physicoChemical::NA*Yi/(molarMass[i]);
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].clamp_min(0);

}

// tmp<fv::convectionScheme<scalar>> mvConvection
// (
//     fv::convectionScheme<scalar>::New
//     (
//         mesh,
//         fields,
//         phi,
//         mesh.divScheme("div(phi,Yi_h)")
//     )
// );


// {
//     reaction->correct();
//     Qdot = reaction->Qdot();
//     volScalarField Yt(0.0*Y[0]);

//     forAll(Y, i)
//     {
//         if (i != inertIndex && composition.active(i))
//         {
//             volScalarField& Yi = Y[i];

//             fvScalarMatrix YiEqn
//             (
//                 fvm::ddt(rho, Yi)
//               + mvConvection->fvmDiv(phi, Yi)
//               - fvm::laplacian(turbulence->muEff(), Yi)
//              ==
//                 reaction->R(Yi)
//               + fvOptions(rho, Yi)
//             );
//             // Info<< "TEEEEEEEEEEEEEEEEEEEEEEEEST" << endl;

//             YiEqn.relax();

//             fvOptions.constrain(YiEqn);

//             YiEqn.solve("Yi");

//             fvOptions.correct(Yi);

//             Yi.clamp_min(0);
//             Yt += Yi;
//         }
//     }

//     Y[inertIndex] = scalar(1) - Yt;
//     Y[inertIndex].clamp_min(0);
// }