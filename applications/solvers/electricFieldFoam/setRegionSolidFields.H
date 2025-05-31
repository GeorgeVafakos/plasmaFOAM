fvMesh& mesh = solidRegions[i];
IOdictionary& physicalProperties = physicalPropertiesSolid[i];
volScalarField& epsilon = epsilonSolid[i];
volScalarField& voltExtAmp = voltExtAmpSolid[i];
volScalarField& voltExt = voltExtSolid[i];
volScalarField& voltInd = voltIndSolid[i];
volScalarField& volt = voltSolid[i];
volVectorField& EExtAmp = EExtAmpSolid[i];
volVectorField& EExt = EExtSolid[i];
volVectorField& EInd = EIndSolid[i];
volVectorField& E = ESolid[i];

// Suppress unused variable warnings
(void)mesh;
(void)physicalProperties;
(void)epsilon;
(void)voltExtAmp;
(void)voltExt;
(void)voltInd;
(void)volt;
(void)EExtAmp;
(void)EExt;
(void)EInd;
(void)E;
