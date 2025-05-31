fvMesh& mesh = solidRegions[i];
IOdictionary& physicalProperties = physicalPropertiesDiel[i];
volScalarField& epsilon = epsilonDiel[i];
volScalarField& voltInd = voltIndDiel[i];
volVectorField& EInd = EIndDiel[i];

// Suppress unused variable warnings
(void)mesh;
(void)physicalProperties;
(void)epsilon;
(void)voltInd;
(void)EInd;