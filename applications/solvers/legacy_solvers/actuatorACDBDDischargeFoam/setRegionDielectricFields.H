fvMesh& mesh = solidRegions[i];
IOdictionary& physicalProperties = physicalPropertiesDiel[i];
volScalarField& epsilon = epsilonDiel[i];
volScalarField& voltExtAmp = voltExtAmpDiel[i];
volScalarField& voltInd = voltIndDiel[i];
volScalarField& volt = voltDiel[i];
volVectorField& EExtAmp = EExtAmpDiel[i];
volVectorField& EInd = EIndDiel[i];
volVectorField& E = EDiel[i];

// Suppress unused variable warnings
(void)mesh;
(void)physicalProperties;
(void)epsilon;
(void)voltExtAmp;
(void)voltInd;
(void)volt;
(void)EExtAmp;
(void)EInd;
(void)E;
