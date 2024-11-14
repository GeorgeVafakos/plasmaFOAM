/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    actuatorACDBDIonicFoam

Description
    Solver for the ionic motion for a single period of an AC-DBD actuator.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createMeshD.H"
    #include "createMeshE.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createFieldsSolid.H"
    

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting iteration loop\n" << endl;

    while (runTime.loop())
    {
        if (runTime.timeIndex() % printScreenResults == 0 || runTime.timeIndex() == 1)
        {
            Info<< "Time = " << runTime.timeName() << "  Time step = " << runTime.timeIndex() << nl << endl;
            solverPerformance::debug = 1;
        }

        // Control time step according to Co num
        #include "CourantNo.H"
        #include "setDeltaT.H" 

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Equations for the Induced Electric Field
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // Reset counters
        counter = 0;
        // solverPerformance::debug = 1;

        while ((convVoltArho>0 || convVoltDrho>0 || convVoltErho>0 || convRhoq>0) && counter<1)
        {
            // Air
            Foam::solverPerformance solvPerfVoltArho = solve 
            (
                fvm::laplacian(voltArho) + rhoq/e0
            );
            convVoltArho = solvPerfVoltArho.nIterations();

            // Dielectric
            Foam::solverPerformance solvPerfVoltDrho = solve
            (
                fvm::laplacian(voltDrho)
            );
            convVoltDrho = solvPerfVoltDrho.nIterations();

            // Encapsulator
            Foam::solverPerformance solvPerfVoltErho = solve
            (
                fvm::laplacian(voltErho)
            );
            convVoltErho = solvPerfVoltErho.nIterations();

            counter++;
            solverPerformance::debug = 0;
            if (counter % 50000 == 0 && (runTime.timeIndex() % printScreenResults == 0 || runTime.timeIndex() == 1))
            {
                Info<< "Current Loop = " << counter << endl;
                solverPerformance::debug = 1;
            }
        }

        if (runTime.timeIndex() % printScreenResults == 0 || runTime.timeIndex() == 1)
        {
            Info<< "Region Inner Loops = " << counter << endl;
            solverPerformance::debug = 1;
        }
        
        // Calculate total electric potential in all regions
        voltAext = voltAextMag*Foam::sin(2*M_PI*(1.0/endTime)*runTime.value());
        voltDext = voltDextMag*Foam::sin(2*M_PI*(1.0/endTime)*runTime.value());
        voltEext = voltEextMag*Foam::sin(2*M_PI*(1.0/endTime)*runTime.value());
        
        voltA = voltAext + voltArho;
        voltD = voltDext + voltDrho;
        voltE = voltEext + voltErho;


        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Charge Transport
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // Update rhoFlux
        rhoqFlux = -sign(linearInterpolate(rhoq))*muc*mesh.magSf()*fvc::snGrad(voltA);

        // Solve the charge transport equation
        Foam::solverPerformance solvPerfRhoq = solve
        (
            fvm::ddt(rhoq) + fvm::div(rhoqFlux, rhoq) - fvm::laplacian(Dc, rhoq)
        );
        convRhoq = solvPerfRhoq.nIterations();


        // Calculate the electric field
        EAext = -fvc::grad(voltAext);
        EArho = -fvc::grad(voltArho);
        EA    = -fvc::grad(voltA);
        EDext = -fvc::grad(voltDext);
        EDrho = -fvc::grad(voltDrho);
        ED    = -fvc::grad(voltD);
        EEext = -fvc::grad(voltEext);
        EErho = -fvc::grad(voltErho);
        EE    = -fvc::grad(voltE);
        Fc = rhoq*EA;


        #include "resetStreamer.H"
        #include "writeCustomTime.H"

        //runTime.write();

        solverPerformance::debug = 0;
        if (runTime.timeIndex() % printScreenResults == 0 || runTime.timeIndex() == 1)
        {
            // Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            //     << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            //     << nl << nl << endl;
            runTime.printExecutionTime(Info);
        }
    }

    runTime.writeAndEnd();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
