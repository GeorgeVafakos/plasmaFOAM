bool LTS = fv::localEulerDdt::enabled(meshAir);

tmp<volScalarField> trDeltaT;

if (LTS)
{
    Info<< "Using LTS" << endl;

    trDeltaT = tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                fv::localEulerDdt::rDeltaTName,
                runTime.timeName(),
                meshAir,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            meshAir,
            dimensionedScalar(dimless/dimTime, 1),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );
}
