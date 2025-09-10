// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Equations for the time varying Induced Electric Field 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Reset counters
regionLoopCounter = 0;
if (runTime.timeIndex() % printScreenResults == 0 || runTime.timeIndex() == 1)
{
    solverPerformance::debug = 1;
    Info<< "Current Iteration = " << 1 << endl;
}

// If rhoq is zero, skip the Poisson equation
if (max(mag(rhoq)).value() > SMALL)
{
    // Ensure that the while loop will run at least once
    voltEqnIter.clear();
    voltEqnIter.push_back(1);

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

    // Calculate the induced electric field
    EInd = -fvc::grad(voltInd);
    forAll(solidRegions, i)
    {
        #include "setRegionSolidFields.H"
        EInd = -fvc::grad(voltInd);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Clamp Induced Electric Field 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
EInd = min(scalar(1.0), maxVal/(mag(EInd)+dimensionedScalar("smallVal", EInd.dimensions(), SMALL))) * EInd;
forAll(solidRegions, i)
{
    #include "setRegionSolidFields.H"
    EInd = min(scalar(1.0), maxVal/(mag(EInd)+dimensionedScalar("smallVal", EInd.dimensions(), SMALL))) * EInd;
}














// EInd = min(mag(U), 50) * (U / (mag(U) + SMALL));

// volScalarField magE("magE", mag(EInd));

// dimensionedScalar E_min("E_min", magE.dimensions(), 1e4);
// dimensionedScalar E_max("E_max", magE.dimensions(), 5e6);

// volScalarField magE_limited("magE_limited", min(max(magE, E_min), E_max));

// Calculate the induced electric field
dimensionedScalar maxVal("maxVal", EInd.dimensions(), 50.0);
EInd = min(mag(EInd), maxVal) * (EInd / (mag(EInd) + SMALL));
forAll(solidRegions, i)
{
    #include "setRegionSolidFields.H"
    EInd = min(mag(EInd), maxVal) * (EInd / (mag(EInd) + SMALL));
}


// volVectorField EInd_limited
// (
//     IOobject
//     (
//         "E_limited", 
//         runTime.timeName(), 
//         mesh
//     ),
//     E * (magE_limited / (magE ))
// );

// dimensionedScalar E_dimless_min("E_dimless_min", dimless, 1e-2);
// dimensionedScalar E_dimless_max("E_dimless_max", dimless, 5e6);