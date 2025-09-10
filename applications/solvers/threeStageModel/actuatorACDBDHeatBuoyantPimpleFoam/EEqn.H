{
    volScalarField& he = thermo.he();
    //he += heStr;
    
    /*
    //if ((runTime.timeIndex())%numSteps==1)
    if (runTime.time().value() >= (int(runTime.value()/forceTimes[forceTimes.size()-1])+0)*Tp + (startTime+timeBetweenDischarges*dischargeCounter) && dischargeCounter<numDischarges)
    {
        he += heStr;
        dischargeCounter++;
        Info<< "TRUUUUUEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << endl;
    }
    */
    
    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvm::div(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian(turbulence->alphaEff(), he)
     ==
        rho*(U&g)
      + radiation->Sh(thermo, he)
      + heStr*fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
    radiation->correct();
}
