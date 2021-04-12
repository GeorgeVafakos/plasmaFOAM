// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

Info<< "force name = " << forcesNames[metritis]  << nl << endl; 

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho, U) + fvm::div(phi, U)
  + MRF.DDt(rho, U)
  + turbulence->divDevRhoReff(U)
  + fCoulomb[metritis]

 ==
    fvOptions(rho, U)
);


fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
    K = 0.5*magSqr(U);
}
