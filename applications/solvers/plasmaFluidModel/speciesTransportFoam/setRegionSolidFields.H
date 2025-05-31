fvMesh& mesh = solidRegions[i];
IOdictionary& dielectricProperties = dielectricPropertiesSolid[i];
volScalarField& epsilon = epsilonSolid[i];
volScalarField& voltExtAmp = voltExtAmpSolid[i];
volScalarField& voltExt = voltExtSolid[i];
volScalarField& voltInd = voltIndSolid[i];
volScalarField& volt = voltSolid[i];
volVectorField& EExtAmp = EExtAmpSolid[i];
volVectorField& EExt = EExtSolid[i];
volVectorField& EInd = EIndSolid[i];
volVectorField& E = ESolid[i];

// Suppress unused variable warnings for the variables that are not used
// using C-like void type casting. The variables that are used are not affected.
(void)mesh;
(void)dielectricProperties;
(void)epsilon;
(void)voltExtAmp;
(void)voltExt;
(void)voltInd;
(void)volt;
(void)EExtAmp;
(void)EExt;
(void)EInd;
(void)E;
