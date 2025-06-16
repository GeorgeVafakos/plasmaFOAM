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

// Calculate reduced electric field in Td
EN = 1.0e21*mag(E/N);