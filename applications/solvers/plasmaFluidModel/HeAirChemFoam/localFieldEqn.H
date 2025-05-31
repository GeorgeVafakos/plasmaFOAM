// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Equations for the time varying Induced Electric Field 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Reset counters
regionLoopCounter = 0;
solverPerformance::debug = 1;
Info<< "Current Iteration = " << 1 << endl;

// Loop through different regions to solve the Poisson equations for voltInd
while (*std::max_element(voltEqnIter.begin(), voltEqnIter.end()) && regionLoopCounter<maxRegionLoops)
{
    voltEqnIter.clear();

    // Solve gas region
    Foam::solverPerformance solvPerfVolt = solve
    (
        fvm::laplacian(voltInd) + rhoq/constant::electromagnetic::epsilon0
    );
    voltEqnIter.push_back(solvPerfVolt.nIterations());

    // Solve solid regions
    forAll(solidRegions, i)
    {
        #include "setRegionSolidFields.H"
        Foam::solverPerformance solvPerfVoltEqn = solve 
        (
            fvm::laplacian(voltInd)
        );
        voltEqnIter.push_back(solvPerfVoltEqn.nIterations());
    }

    // Print performance at custom iteration intervals
    regionLoopCounter++;
    solverPerformance::debug = 0;
    int printPerformance = 500;
    if (regionLoopCounter % printPerformance == 0)
    {
        Info<< "Current Iteration = " << regionLoopCounter << endl;
        solverPerformance::debug = 1;
    }
}
Info<< "Induced Field: Region Inner Loops = " << regionLoopCounter << endl;

// Calculate the induced electric field
EInd = -fvc::grad(voltInd);
forAll(solidRegions, i)
{
    #include "setRegionSolidFields.H"
    EInd = -fvc::grad(voltInd);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Calculate for the time varying External Electric Field
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Calculate the potential of the external electric field
voltExt = voltHandler->calcVoltage();
forAll(solidRegions, i)
{
    #include "setRegionSolidFields.H"
    voltExtSolid[i] = voltHandlerSolid[i].calcVoltage();
}

// Calculate the external electric field
EExt = -fvc::grad(voltExt);
forAll(solidRegions, i)
{
    #include "setRegionSolidFields.H"
    EExt = -fvc::grad(voltExt);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Calculate local electric field
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
volt = voltExt + voltInd;
E = EExt + EInd;
forAll(solidRegions, i)
{
    #include "setRegionSolidFields.H"
    volt = voltExt + voltInd;
    E = EExt + EInd;
}

// Calculate reduced electrin field in Td
E_N = 1.0e21*mag(E/N);