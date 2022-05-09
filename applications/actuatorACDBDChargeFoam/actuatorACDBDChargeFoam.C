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
    actuatorACDBDChargeFoam

Description
    Solver for charge propagation for a single period in an AC-DBD actuator.

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
    #include "createFieldsSolid.H"
    #include "createTimeControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting iteration loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Control time step according to Co num
        #include "CourantNo.H"
        #include "setDeltaT.H" 

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Equations for the Induced Electric Field
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // Reset counters
        counter = 0;
        solverPerformance::debug = 1;

        while ((convVoltArho>0 || convVoltDrho>0 || convVoltIrho>0 || convRhoq>0) && counter<1)
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

            // Insulator
            Foam::solverPerformance solvPerfVoltIrho = solve
            (
                fvm::laplacian(voltIrho)
            );
            convVoltIrho = solvPerfVoltIrho.nIterations();

            counter++;
            solverPerformance::debug = 0;
            if (counter % 5000 == 0)
            {
                Info<< "Current Loop = " << counter << endl;
                solverPerformance::debug = 1;
            }
        }

        Info<< "Region Inner Loops = " << counter << endl;
        solverPerformance::debug = 1;

        // Calculate total electric potential in all regions
        voltA = voltAext*Foam::sin(2*M_PI*(1.0/endTime)*runTime.value()) + voltArho;
        voltD = voltDext*Foam::sin(2*M_PI*(1.0/endTime)*runTime.value()) + voltDrho;
        voltI = voltIext*Foam::sin(2*M_PI*(1.0/endTime)*runTime.value()) + voltIrho;


        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Charge Transport
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // Update rhoFlux
        rhoqFlux = -k*mesh.magSf()*fvc::snGrad(voltA);

        // Solve the charge transport equation
        Foam::solverPerformance solvPerfRhoq = solve
        (
            fvm::ddt(rhoq) + fvm::div(rhoqFlux, rhoq)
        );
        convRhoq = solvPerfRhoq.nIterations();


        // Calculate the electric field
        EA = -fvc::grad(voltA);
        ED = -fvc::grad(voltD);
        EI = -fvc::grad(voltI);
        Fc = rhoq*EA;


        #include "resetStreamer.H"
        #include "writeCustomTime.H"


        // runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << nl << endl;
    }

    runTime.writeAndEnd();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
