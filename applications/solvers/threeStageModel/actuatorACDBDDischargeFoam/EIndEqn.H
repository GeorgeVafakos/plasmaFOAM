// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Equations for the Induced Electric Field 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Reset counters
regionLoopCounter = 0;
solverPerformance::debug = 1;
Info<< "Current Iteration = " << 1 << endl;

// Loop through different regions to solve equations for voltInd
while (*std::max_element(voltEqnIter.begin(), voltEqnIter.end()) && regionLoopCounter<maxRegionLoops)
{
    voltEqnIter.clear();

    // Solve gas region
    Foam::solverPerformance solvPerfVolt = solve
    (
        fvm::laplacian(voltInd) + rhoq/e0
    );
    voltEqnIter.push_back(solvPerfVolt.nIterations());

    // Calculate the induced electric field
    EInd = -fvc::grad(voltInd);

    // Solve solid regions
    forAll(solidRegions, i)
    {
        #include "setRegionDielectricFields.H"
        Foam::solverPerformance solvPerfVoltEqn = solve 
        (
            fvm::laplacian(voltInd)
        );
        voltEqnIter.push_back(solvPerfVoltEqn.nIterations());

        // Calculate the induced electric field
        EInd = -fvc::grad(voltInd);
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

    // Write current solution
    if (regionLoopCounter % writeRegionLoop == 0)
        runTime.writeNow();
}
Info<< "Induced Field: Region Inner Loops = " << regionLoopCounter << endl;

