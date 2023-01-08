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
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting iteration loop\n" << endl;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // Equations for the External Electric Field
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Reset counters
    counter = 0;
    solverPerformance::debug = 1;

    while ((convVoltAext>0 || convVoltDext>0 ) && counter<100000)
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

    while (runTime.loop())
    {
        Info<< "Iteration = " << runTime.timeName() << nl << endl;


        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Equations for the Induced Electric Field
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // Reset counters
        counter = 0;
        solverPerformance::debug = 1;

        while ((convVoltArho>0 || convVoltDrho>0 ) && counter<100000)
        {
            // Air
            Foam::solverPerformance solvPerfVoltArho = solve 
            (
                fvm::laplacian(voltArho) + (e/epsilon0)*(ni-ne)
            );
            convVoltArho = solvPerfVoltArho.nIterations();

            // Dielectric
            Foam::solverPerformance solvPerfVoltDrho = solve
            (
                fvm::laplacian(voltDrho)
            );
            convVoltDrho = solvPerfVoltDrho.nIterations();

            counter++;
            solverPerformance::debug = 0;
            if (counter % 50 == 0)
            {
                Info<< "Current Loop = " << counter << endl;
                solverPerformance::debug = 1;
            }
        }

        Info<< "Induced Field: Region Inner Loops = " << counter << endl;

        // Calculate total electric potential
        voltA = voltAext + voltArho;
        voltD = voltDext + voltDrho;

        // Calculate the induced electric field
        EArho = -fvc::grad(voltArho);
        EDrho = -fvc::grad(voltDrho);

        // Total electric field
        EA = EAext + EArho;
        ED = EDext + EDrho;

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Specied continuity equation
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // Update face fluxes
        eleFlux =  mue*mesh.magSf()*fvc::snGrad(voltA);
        ionFlux = -mui*mesh.magSf()*fvc::snGrad(voltA);

        // Solve the electron continuity
        fvScalarMatrix eleEqn
        (
            fvm::ddt(ne) + fvm::div(eleFlux, ne) - fvm::laplacian(De, ne) 
            ==
            fvm::Sp(alpha*mag(-mue*EA), ne) - fvm::Sp(r*ni, ne)
        );
        eleEqn.solve();

        // Solve the ion continuity
        fvScalarMatrix ionEqn
        (
            fvm::ddt(ni) + fvm::div(ionFlux, ni) - fvm::laplacian(Di, ni) 
            ==
            alpha*mag(-mue*EA)*ne - fvm::Sp(r*ne, ni)
        );
        ionEqn.solve();


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
