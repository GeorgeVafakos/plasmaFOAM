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
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createMeshD.H"
    #include "createControl.H"
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

    while ((convVoltAext>0 || convVoltDext>0 ) && counter<100)
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
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Control time step according to Co num
        #include "CourantNo.H"
        #include "setDeltaT.H" 


        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Equations for the Induced Electric Field
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // Reset counters
        counter = 0;
        convVoltArho=1;
        convVoltDrho=1;
        solverPerformance::debug = 1;
        

        while ((convVoltArho>0 || convVoltDrho>0 ) && counter<100)
        {
            // Air
            Foam::solverPerformance solvPerfVoltArho = solve 
            (
                fvm::laplacian(voltArho) + (e/epsilon0)*(np-nn-ne)
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
	
	//scalar tt=runTime.value();
	scalar ff= Foam::sin(2.0*14.0e3*3.14159*runTime.value());
        
        // Calculate the induced electric field
        EArho = -fvc::grad(voltArho);
        EDrho = -fvc::grad(voltDrho);

        // Calculate total electric potential
        voltA = voltAext*ff + voltArho;
        voltD = voltDext*ff + voltDrho;

        // Total electric field
        EA = EAext*ff + EArho;
        ED = EDext*ff + EDrho;

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Specied continuity equation
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // Update face fluxes
        #include "coeffs.H"
        eleFlux =  linearInterpolate(mueF)*mesh.magSf()*fvc::snGrad(voltA);
        posFlux = -linearInterpolate(mupF)*mesh.magSf()*fvc::snGrad(voltA);
        negFlux =  linearInterpolate(munF)*mesh.magSf()*fvc::snGrad(voltA);
        
        GammaEle = -mueF*ne*EA - DeF*fvc::grad(ne);
        GammaPos =  mupF*np*EA - DpF*fvc::grad(np);
        GammaNeg = -munF*nn*EA - DnF*fvc::grad(nn);
        SLe = alphaF*mag(GammaEle) - (rep*ne*np + etaF*mag(GammaEle));
        SLp = alphaF*mag(GammaEle) - (rep*ne*np + rnp*nn*np);
        SLn = etaF*mag(GammaEle) - rnp*nn*np;
        

        // Solve the electron continuity
        fvScalarMatrix eleEqn
        (
            fvm::ddt(ne) + fvm::div(eleFlux, ne) - fvm::laplacian(DeF, ne) 
            ==
            SLe
        );
        eleEqn.solve();
        
        // Solve the positive ion continuity
        fvScalarMatrix posEqn
        (
            fvm::ddt(np) + fvm::div(posFlux, np) - fvm::laplacian(DpF, np) 
            ==
            SLp
        );       
        posEqn.solve();

        // Solve the negative ion continuity
        fvScalarMatrix negEqn
        (
            fvm::ddt(nn) + fvm::div(negFlux, nn) - fvm::laplacian(DnF, nn)
            ==
            SLn 
        );
        negEqn.solve();

        
        ne = max(ne,ne0);
        np = max(np,0*ne0);
        nn = max(nn,0*ne0);
        // T.clamp_min(dimensionedScalar("Tmin", dimTemperature, 300.0));
        
        rhoq = e*(np-ne-nn);
        Fc = rhoq*EA;

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
