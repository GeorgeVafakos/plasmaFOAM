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
    actuatorElectrostatic2RegionsFoam

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
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createFieldsSolid.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting iteration loop\n" << endl;

    while (runTime.loop() && (convVoltA>0 || convVoltD>0 || convRhoq>0))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Iteration counter
        iterCount++;

        // Automatic calculation of DeltaT
        #include "CourantNo.H"
        #include "setDeltaT.H" 

        // Calculate the surface charge on the dielectric surface
        gradVoltA=fvc::snGrad(voltA);
        label patchID = mesh.boundaryMesh().findPatchID("region0_to_dielectric");
        if (iterCount>2)
            epsilonA.boundaryFieldRef()[patchID] == 1./erD.value()*(unit.boundaryField()[patchID]-1/e0.value()/mesh.surfaceInterpolation::deltaCoeffs()[patchID]*rhoq.boundaryField()[patchID]/gradVoltA.boundaryField()[patchID]);

        // Solve the Poisson-Boltzmann equation in air
        Foam::solverPerformance solvPerfVoltA = solve 
        (
            fvm::laplacian(voltA) + rhoq/e0
        );
        convVoltA = solvPerfVoltA.nIterations();

        // Solve the Laplace equation in the solid
        Foam::solverPerformance solvPerfVoltD = solve
        (
            fvm::laplacian(voltD)
        );
        convVoltD = solvPerfVoltD.nIterations();
        
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
