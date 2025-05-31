fvMesh& mesh = solidRegions[i];
IOdictionary& physicalProperties = physicalPropertiesDiel[i];
volScalarField& epsilon = epsilonDiel[i];
volScalarField& voltExtAmp = voltExtAmpDiel[i];
volScalarField& voltExt = voltExtDiel[i];
volScalarField& voltInd = voltIndDiel[i];
volScalarField& voltIndDischarge = voltIndDischargeDiel[i];
volScalarField& volt = voltDiel[i];
volVectorField& EExtAmp = EExtAmpDiel[i];
volVectorField& EExt = EExtDiel[i];
volVectorField& EInd = EIndDiel[i];
volVectorField& E = EDiel[i];

// Suppress unused variable warnings
(void)mesh;
(void)physicalProperties;
(void)epsilon;
(void)voltExtAmp;
(void)voltExt;
(void)voltInd;
(void)voltIndDischarge;
(void)volt;
(void)EExtAmp;
(void)EExt;
(void)EInd;
(void)E;
