// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Define voltage related objects
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Read voltageType from dictionary
word voltageType = word(physicalProperties.lookup("voltageType")); 

// Pointers for the voltageHandler class objects
autoPtr<voltageHandler> voltHandler;
PtrList<voltageHandler> voltHandlerSolid(solidRegions.size());

// Populate region0 pointer and solid field pointer lists
if (voltageType == "ac" || voltageType == "AC") 
{
    voltHandler.reset
    (
        new voltageAC
        (
            mesh, 
            physicalProperties, 
            voltExt, 
            voltExtAmp
        ) 
    ); 
    
    forAll(solidRegions, i)
    {
        voltHandlerSolid.set
        (
            i,
            new voltageAC
            (
                solidRegions[i],
                physicalProperties,
                voltExtSolid[i],
                voltExtAmpSolid[i]
            )
        );
    }
}
else if (voltageType == "pulsed")
{
    voltHandler.reset
    (
        new voltagePulsed
        (
            mesh, 
            physicalProperties, 
            voltExt, 
            voltExtAmp
        ) 
    ); 

    forAll(solidRegions, i)
    {
        voltHandlerSolid.set
        (
            i,
            new voltagePulsed
            (
                solidRegions[i],
                physicalProperties,
                voltExtSolid[i],
                voltExtAmpSolid[i]
            )
        );
    }
}
else if (voltageType == "dc" || voltageType == "DC")
{
    voltHandler.reset
    (
        new voltageDC
        (
            mesh, 
            physicalProperties, 
            voltExt, 
            voltExtAmp
        ) 
    ); 

    forAll(solidRegions, i)
    {
        voltHandlerSolid.set
        (
            i,
            new voltageDC
            (
                solidRegions[i],
                physicalProperties,
                voltExtSolid[i],
                voltExtAmpSolid[i]
            )
        );
    }
}
else
{
    FatalErrorInFunction
        << "voltageType: " << voltageType
        << " is not a valid choice. "
        << "Options are: pulsed, ac, dc"
        << abort(FatalError);
}
