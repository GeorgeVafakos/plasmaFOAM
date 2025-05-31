// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Equations for the External Electric Field
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//- Calculate the amplitude field value across all regions.
//- volFields fields ending in <Amp> denote the amplitude.

Info<< "\nStarting external electric field loop\n" << endl;

// Reset counters
regionLoopCounter = 0;
solverPerformance::debug = 1;
Info<< "Current Iteration = " << 1 << endl;

// Loop through different regions to solve the Laplace equations for voltExtAmp
while (*std::max_element(voltEqnIter.begin(), voltEqnIter.end()) && regionLoopCounter<maxRegionLoops)
{
    voltEqnIter.clear();

    // Solve gas region
    Foam::solverPerformance solvPerfVoltEqn = solve
    (
        fvm::laplacian(voltExtAmp)
    );
    voltEqnIter.push_back(solvPerfVoltEqn.nIterations());

    // Solve solid regions
    forAll(solidRegions, i)
    {
        #include "setRegionSolidFields.H"
        Foam::solverPerformance solvPerfVoltEqn = solve 
        (
            fvm::laplacian(voltExtAmp)
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
Info<< "External Field Amplitude: Region Inner Loops = " << regionLoopCounter << endl;

// Calculate the external electric field amplitude
EExtAmp = -fvc::grad(voltExtAmp);
forAll(solidRegions, i)
{
    #include "setRegionSolidFields.H"
    EExtAmp = -fvc::grad(voltExtAmp);
}

