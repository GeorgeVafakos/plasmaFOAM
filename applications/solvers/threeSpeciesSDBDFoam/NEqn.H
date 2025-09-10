// Update coefficients
#include "coeffs.H"

// Calculate species fluxes
eleFlux =  linearInterpolate(mueF)*mesh.magSf()*fvc::snGrad(volt);
posFlux = -linearInterpolate(mupF)*mesh.magSf()*fvc::snGrad(volt);
negFlux =  linearInterpolate(munF)*mesh.magSf()*fvc::snGrad(volt);

// Calculate reaction source terms
GammaEle = -mueF*ne*E - DeF*fvc::grad(ne);
GammaPos =  mupF*np*E - DpF*fvc::grad(np);
GammaNeg = -munF*nn*E - DnF*fvc::grad(nn);
SLe = (alphaF-etaF)*mag(GammaEle) - rep*ne*np;
SLp = alphaF*mag(GammaEle) - (rep*ne*np + rnp*nn*np);
SLn = etaF*mag(GammaEle) - rnp*nn*np;

// Solve the electron continuity
fvScalarMatrix eleEqn
(
    fvm::ddt(ne) + fvm::div(eleFlux, ne) - fvm::laplacian(DeF, ne) 
    ==
    SLe
);
eleEqn.solve();

// Solve for the positive ion continuity
fvScalarMatrix posEqn
(
    fvm::ddt(np) + fvm::div(posFlux, np) - fvm::laplacian(DpF, np) 
    ==
    SLp
);       
posEqn.solve();

// Solve for the negative ion continuity
fvScalarMatrix negEqn
(
    fvm::ddt(nn) + fvm::div(negFlux, nn) - fvm::laplacian(DnF, nn)
    ==
    SLn 
);
negEqn.solve();

// Ensure that species don't fall under a certain threshold
ne.clamp_min(dimensionedScalar("neMin", ne.dimensions(), 1.0e+13));

// Calculate the space charge density
rhoq = constant::electromagnetic::e*(np-ne-nn);

// Calculate the force
Fc = rhoq*E;
