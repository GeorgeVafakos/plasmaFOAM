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
    plasmaChemFoam

Description
    Solver for the local electric field.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

#include "regionProperties.H"
#include "voltageHandler.H"
#include "transportCoeffsHandler.H"
#include "transportCoeffsBolsig.H"
#include "interpolation2DTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for the local electric field."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createVoltExtObjs.H"
    // #include "createTransportCoeffs.H"
    // #include "readBolsigData.H"

    // Calculate amplitude of potenetial of external field
    #include "extFieldEqn.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;
    while( runTime.loop() )
    {

        
        // Calculate local electric field
        #include "localFieldEqn.H"

        // Control time step according to Co num
        #include "CourantNo.H"
        #include "setDeltaT.H" 

        ++runTime;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "calcCoeffs.H"

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Equations for species continuity
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        /*

        #include "reactionTerms.H"
        
        // #include "YEqn.H"
        #include "NEqn.H"

        
        // // Calculate rhoq
        // rhoq = 0.0*rhoq;
        // forAll(composition.species(), i)
        // {
        //     // rhoq += chargeNumber[i]*constant::electromagnetic::e*rhoGas*constant::physicoChemical::NA*Y[i]/(molarMass[i]);
        //     rhoq += chargeNumber[i]*constant::electromagnetic::e*n[i];
        // }

        

        // surfaceCharge = linearInterpolate(rhoq);// & mesh.Sf();

        // volt.write();

        */
        // break;

        // label eIndex = composition.species().find("e");
        // volScalarField& mob_e = mobilityCoeffSpecies[eIndex];

        // forAll(mob_e, celli)
        // {
        //     scalar E_value  = E_N[celli];
        //     scalar He_value = HeF[celli];

        //     // // Optional: Clamp inputs to interpolation table limits
        //     // E_value  = max(min(E_value,  mobilityTable.xMax()), mobilityTable.xMin());
        //     // He_value = max(min(He_value, mobilityTable.yMax()), mobilityTable.yMin());

        //     mob_e[celli] = mobilityTable(He_value, E_value);
        // }

        mobilityCoeffSpecies[composition.species().find("e")].write();
        

        runTime.printExecutionTime(Info);
        runTime.write();
        // mobilityCoeffSpecies[composition.species().find("e")].write();
        // diffusionCoeffSpecies[composition.species().find("e")].write();
    }

    runTime.writeAndEnd();

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //