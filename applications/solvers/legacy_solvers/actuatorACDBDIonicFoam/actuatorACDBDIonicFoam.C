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
#include "regionProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "addCheckCaseOptions.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    

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
        regionLoopCounter = 0;
        // solverPerformance::debug = 1;

        while (*std::max_element(voltIndIter.begin(), voltIndIter.end()) && regionLoopCounter<1)
        {
            voltIndIter.clear();

            // Solve gas regions
            Foam::solverPerformance solvPerfVolt = solve
            (
                fvm::laplacian(voltInd) + rhoq/e0
            );
            voltIndIter.push_back(solvPerfVolt.nIterations());

            // Solve dielectric regions
            forAll(solidRegions, i)
            {
                #include "setRegionDielectricFields.H"
                Foam::solverPerformance solvPerfVolt = solve 
                (
                    fvm::laplacian(voltInd)
                );
                voltIndIter.push_back(solvPerfVolt.nIterations());
            }

            // Print performance at custom iteration intervals
            regionLoopCounter++;
            solverPerformance::debug = 0;
            int printPerformance = 500;
            if (regionLoopCounter % printPerformance == 0 && (runTime.timeIndex() % printScreenResults == 0 || runTime.timeIndex() == 1))
            {
                Info<< "Current Loop = " << regionLoopCounter << endl;
                solverPerformance::debug = 1;
            }
        }

        if (runTime.timeIndex() % printScreenResults == 0 || runTime.timeIndex() == 1)
        {
            Info<< "Region Inner Loops = " << regionLoopCounter << endl;
            solverPerformance::debug = 1;
        }
        

        // Calculate total electric potential and field in all regions
        voltExt = voltExtAmp*Foam::sin(2*M_PI*(1.0/endTime)*runTime.value());
        volt = voltExt + voltInd;
        EExt = -fvc::grad(voltExt);
        EInd = -fvc::grad(voltInd);
        E = -fvc::grad(volt);
        forAll(solidRegions, i)
        {
            #include "setRegionDielectricFields.H"
            voltExt = voltExtAmp*Foam::sin(2*M_PI*(1.0/endTime)*runTime.value());
            volt = voltExt + voltInd;
            EExt = -fvc::grad(voltExt);
            EInd = -fvc::grad(voltInd);
            E = -fvc::grad(volt);
        }


        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Charge Transport
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Update rhoFlux
        rhoqFlux = -sign(linearInterpolate(rhoq))*muc*mesh.magSf()*fvc::snGrad(volt);

        // Solve the charge transport equation
        solve
        (
            fvm::ddt(rhoq) + fvm::div(rhoqFlux, rhoq) - fvm::laplacian(Dc, rhoq)
        );

        // Calculate the EHD force
        Fc = rhoq*E;

        #include "addNewDischarge.H"
        #include "writeCustomTime.H"

        solverPerformance::debug = 0;
        if (runTime.timeIndex() % printScreenResults == 0 || runTime.timeIndex() == 1)
        {
            runTime.printExecutionTime(Info);
        }
    }

    runTime.writeAndEnd();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
