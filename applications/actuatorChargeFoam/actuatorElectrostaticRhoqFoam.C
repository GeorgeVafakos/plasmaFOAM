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
    actuatorFoamElStatic

Description
    Solver for electrostatics in a DBD actuator.

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
    // #include "createTimeControls.H"
    #include "createFields.H"
    #include "createFieldsSolid.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting iteration loop\n" << endl;

    while (runTime.loop() && (convVoltA>0 || convVoltD>0 || convVoltI>0 || convVoltR>0 ))
    {
        Info<< "Iteration = " << runTime.timeName() << nl << endl;

        // Iteration counter
        iterCount++;

        // Control time step according to Co num
        // #include "CourantNo.H"
        // #include "setDeltaT.H" 


        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Poisson Equations
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // Air
        Foam::solverPerformance solvPerfVoltA = solve 
        (
            fvm::laplacian(voltA) 
        );
        convVoltA = solvPerfVoltA.nIterations();

        // Dielectric
        Foam::solverPerformance solvPerfVoltD = solve
        (
            fvm::laplacian(voltD)
        );
        convVoltD = solvPerfVoltD.nIterations();

        // Insulator
        Foam::solverPerformance solvPerfVoltI = solve
        (
            fvm::laplacian(voltI)
        );
        convVoltI = solvPerfVoltI.nIterations();

        // Field from rhoq
        Foam::solverPerformance solvPerfVoltR = solve 
        (
            fvm::laplacian(voltR) + rhoq/e0
        );
        convVoltR = solvPerfVoltR.nIterations();


        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Charge Transport
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // // Update rhoFlux
        // rhoqFlux = -k*mesh.magSf()*fvc::snGrad(voltA);

        // // Solve the charge transport equation
        // Foam::solverPerformance solvPerfRhoq = solve
        // (
        //     fvm::ddt(rhoq) + fvm::div(rhoqFlux, rhoq)
        // );
        // convRhoq = solvPerfRhoq.nIterations();


        // Calculate the electric field
        EA = -fvc::grad(voltA);
        ED = -fvc::grad(voltD);
        EI = -fvc::grad(voltI);
        ER = -fvc::grad(voltR);
        // Fc = rhoq*EA;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << nl << endl;
    }

    runTime.writeAndEnd();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
