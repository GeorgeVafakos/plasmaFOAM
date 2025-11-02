// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Equations for the Induced Electric Field 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Reset counters
regionLoopCounter = 0;
if (runTime.timeIndex() % printScreenResults == 0 || runTime.timeIndex() == 1)
{
    Info<< "Current Iteration = " << 1 << endl;
    solverPerformance::debug = 1;
}

voltEqnIter.clear();
voltEqnIter.push_back(1);

// Loop through different regions to solve equations for voltInd
while (*std::max_element(voltEqnIter.begin(), voltEqnIter.end()) && regionLoopCounter<maxRegionLoopsIndEqn)
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
        #include "setRegionDielectricFields.H"
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
    if (regionLoopCounter % printPerformance == 0 && (runTime.timeIndex() % printScreenResults == 0 || runTime.timeIndex() == 1))
    {
        Info<< "Current Iteration = " << regionLoopCounter << endl;
        solverPerformance::debug = 1;
    }    
}

if (runTime.timeIndex() % printScreenResults == 0 || runTime.timeIndex() == 1)
{
    Info<< "Region Inner Loops = " << regionLoopCounter << endl;
    solverPerformance::debug = 1;
}

