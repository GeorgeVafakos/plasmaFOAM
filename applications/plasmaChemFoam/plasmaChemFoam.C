/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    plasmaChemFoam

Description
    Solver streamer plasma chemistry.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createMeshD.H"
    #include "createMeshI.H"
    #include "createFields.H"
    #include "createTimeControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting iteration loop\n" << endl;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // Equations for the External Electric Field
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Reset counters
    counter = 0;
    solverPerformance::debug = 1;

    while ((convVoltAext>0 || convVoltDext>0 || convVoltIext>0 ) && counter<10000)
    {
        // Air
        Foam::solverPerformance solvPerfVoltAext = solve 
        (
            fvm::laplacian(voltAext)
        );
        convVoltAext = solvPerfVoltAext.nIterations();

        // Dielectric
        Foam::solverPerformance solvPerfVoltDext = solve
        (
            fvm::laplacian(voltDext)
        );
        convVoltDext = solvPerfVoltDext.nIterations();

        // Insulator
        Foam::solverPerformance solvPerfVoltIext = solve
        (
            fvm::laplacian(voltIext)
        );
        convVoltIext = solvPerfVoltIext.nIterations();

        counter++;
        solverPerformance::debug = 0;
        if (counter % 50 == 0)
        {
            Info<< "Current Loop = " << counter << endl;
            solverPerformance::debug = 1;
        }
    }

    Info<< "External Field: Region Inner Loops = " << counter << endl;

    // Calculate the external electric field
    EAext = -fvc::grad(voltAext);
    EDext = -fvc::grad(voltDext);
    EIext = -fvc::grad(voltIext);


    while (runTime.loop())
    {
        // alpha = dimAlpha*( (voltAext/dimVolt)*pos(voltAext/dimVolt-2.0e3) +  (voltAext/dimVolt)*pos(10.0e3 - voltAext/dimVolt));
        // alpha = dimAlpha*(  (voltAext/dimVolt)*pos(10.0e3-voltAext/dimVolt)*pos((voltAext/dimVolt)*pos(10.0e3-voltAext/dimVolt) - 5.0e3)    );
        // alpha = press*dimAlpha*sqrt(mag(mag(EAext)/(press*dimE)));
        // alpha = voltAext*pos(voltAext-5.0e3*cond);

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Control time step according to Co num
        #include "CourantNo.H"
        #include "setDeltaT.H" 

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Equations for the Induced Electric Field
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // Reset counters
        counter = 0;
        convVoltAind=1;
        convVoltDind=1;
        convVoltIind=1;
        solverPerformance::debug = 1;

        while ((convVoltAind>0 || convVoltDind>0 || convVoltIind>0 ) && counter<50)
        {
            // Air
            Foam::solverPerformance solvPerfVoltAind = solve 
            (
                fvm::laplacian(voltAind) + (e/epsilon0)*(np-ne-nn)
            );
            convVoltAind = solvPerfVoltAind.nIterations();

            // Dielectric
            Foam::solverPerformance solvPerfVoltDind = solve
            (
                fvm::laplacian(voltDind)
            );
            convVoltDind = solvPerfVoltDind.nIterations();

            // Insulator
            Foam::solverPerformance solvPerfVoltIind = solve
            (
                fvm::laplacian(voltIind)
            );
            convVoltIind = solvPerfVoltIind.nIterations();

            counter++;
            solverPerformance::debug = 0;
            if (counter % 50 == 0)
            {
                Info<< "Current Loop = " << counter << endl;
                solverPerformance::debug = 1;
            }
        }

        Info<< "Induced Field: Region Inner Loops = " << counter << endl;
        solverPerformance::debug = 1;

        // Calculate the induced electric field
        EAind = -fvc::grad(voltAind);
        EDind = -fvc::grad(voltDind);
        EIind = -fvc::grad(voltIind);

        // Calculate total electric potential
        voltA = voltAext + voltAind;
        voltD = voltDext + voltDind;
        voltI = voltIext + voltIind;

        // Total electric field
        EA = EAext + EAind;
        ED = EDext + EDind;
        EI = EIext + EIind;

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Calculation of coefficients
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        alpha = press*dimAlpha*( 0.21*(4.71e-11*pow(mag(EA)/(press*dimE),3.0)*pos(1.4e4 - (mag(EA)/(press*dimE)))
                               + 3.32*sqrt(mag(mag(EA)/(press*dimE) - 12500.0))*pos((mag(EA)/(press*dimE))-1.4e4))

                               + 0.79*(1.17e-10*pow(mag(EA)/(press*dimE),3.0)*pos(1.1e4 - (mag(EA)/(press*dimE)))
                               + (mag(0.0319*mag(EA)/(press*dimE)-211.0)*pos(2.1e4-mag(0.0319*mag(EA)/(press*dimE)-211.0)))*pos(mag(0.0319*mag(EA)/(press*dimE)-211.0)*pos(2.1e4-mag(0.0319*mag(EA)/(press*dimE)-211.0)) - 1.1e4)
                               + 6.32*sqrt(mag(mag(EA)/(press*dimE)-16300.0))*pos((mag(EA)/(press*dimE))-2.1e4))
                               );
                
        mue = (1.0/press)*dimMu*(0.21*(24.32*exp(-mag(EA)/(1057*press*dimE)) + 19.38*exp(-mag(EA)/(23430*press*dimE)) + 14.45)
                                +0.79*(173.1*exp(-mag(EA)/(195.1*press*dimE)) + 36.19*exp(-mag(EA)/(12763*press*dimE)) + 31.73)
                                );

        mup = (1.0/press)*dimMu*(0.21*(0.05492*exp(-mag(EA)/(6858*press*dimE)) + 0.07509*exp(-mag(EA)/(38175*press*dimE)) + 0.0308)
                                +0.79*(0.06841*exp(-mag(EA)/(59678*press*dimE)) + 0.09194*exp(-mag(EA)/(12763*press*dimE)) + 0.0320)
                                );

        mun = (1.0/press)*dimMu*(0.181225);

        heta = press*dimAlpha*( 0.21*(1.307 + (33200/(mag(EA)/(press*dimE)))*exp(pow(log(mag(EA)/(press*dimE))-9.04,2)/2.53)) );
        
        De = kB*Te*mue/e;

        Dp = kB*Tamb*mup/e;

        Dn = kB*Tamb*mun/e;

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Species continuity equation
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // Update face fluxes
        eleFlux =  linearInterpolate(mue)*mesh.magSf()*fvc::snGrad(voltA);
        posFlux = -linearInterpolate(mup)*mesh.magSf()*fvc::snGrad(voltA);
        negFlux =  linearInterpolate(mun)*mesh.magSf()*fvc::snGrad(voltA);

        GammaEle = -mue*ne*EA - De*fvc::grad(ne);
        GammaPos = mup*np*EA;// - Dp*fvc::grad(np);
        GammaNeg = -mun*nn*EA;// - Dp*fvc::grad(nn);

        // Solve the electron continuity
        fvScalarMatrix eleEqn
        (
            fvm::ddt(ne) + fvm::div(eleFlux, ne) - fvm::laplacian(De, ne) 
            ==
            alpha*mag(GammaEle) - rep*ne*np - heta*mag(GammaEle)
        );
        eleEqn.solve();

        // Solve the positive ion continuity
        fvScalarMatrix posEqn
        (
            fvm::ddt(np) + fvm::div(posFlux, np)// - fvm::laplacian(Dp, np) 
            ==
            alpha*mag(GammaEle) - rep*ne*np - rnp*nn*np
        );
        posEqn.solve();

        // Solve the negative ion continuity
        fvScalarMatrix negEqn
        (
            fvm::ddt(nn) + fvm::div(negFlux, nn)// - fvm::laplacian(Dn, nn) 
            ==
            heta*mag(GammaEle) - rnp*nn*np
        );
        negEqn.solve();


        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    runTime.writeAndEnd();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
