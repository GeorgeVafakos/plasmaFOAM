/*
volScalarField magE = mag(EA);

forAll(EA, celli)
{
scalar	EtoPs = magE[celli]/760.0;
alphaFF[celli] = 760.0*( 0.21 * ((4.71e-11*Foam::pow(EtoPs,3.0)*(EtoPs<1.4e4))+
				3.32*Foam::sqrt(mag(EtoPs-1.25e4))*(EtoPs>1.4e4)) +
			   0.79 * ((1.17e-10*Foam::pow(EtoPs,3.0))*(EtoPs<1.1e4)+
				(0.0319*EtoPs - 211.0)*(EtoPs>1.1e4)*(EtoPs<2.1e4)+
				6.32*Foam::sqrt(mag(EtoPs-1.63e4))*(EtoPs>2.1e4)) );
				
mueFF[celli] = (1.0/760.0)*(0.21*(24.32*Foam::exp(-EtoPs/1.057e3) + 19.38*Foam::exp(-EtoPs/2.3430e4) + 14.45)+
                          0.79*(173.1*Foam::exp(-EtoPs/195.1)   + 36.19*Foam::exp(-EtoPs/1.2763e4) + 31.73)
                                );
mupFF[celli] = (1.0/760.0)*(0.79*(0.05492*Foam::exp(-EtoPs/6.858e3)  + 0.07509*Foam::exp(-EtoPs/3.8175e4) + 0.0308)+
                          0.21*(0.06841*Foam::exp(-EtoPs/5.9678e4) + 0.09194*Foam::exp(-EtoPs/1.2763e4) + 0.0320)
                                );                                
etaFF[celli] = 760.0*0.21*(1.307 + 3.32e4/EtoPs*Foam::exp(Foam::pow(Foam::log(EtoPs)-9.04,2)/2.53) );                              
        
}*/

volScalarField EtoP = mag(E)/(press*dimE);
alphaF = press*dimAlpha*( 0.21 * ((4.71e-11*Foam::pow(EtoP,3.0)*(EtoP<1.4e4*ones))+
				3.32*Foam::sqrt(mag(EtoP-1.25e4*ones))*(EtoP>1.4e4*ones)) +
			   0.79 * ((1.17e-10*Foam::pow(EtoP,3.0))*(EtoP<1.1e4*ones)+
				(0.0319*EtoP - 211.0*ones)*(EtoP>1.1e4*ones)*(EtoP<2.1e4*ones)+
				6.32*Foam::sqrt(mag(EtoP-1.63e4*ones))*(EtoP>2.1e4*ones)) );
mueF = (1.0/press)*dimMu*(0.21*(24.32*Foam::exp(-EtoP/1.057e3) + 19.38*Foam::exp(-EtoP/2.3430e4) + 14.45*ones)+
                          0.79*(173.1*Foam::exp(-EtoP/195.1)   + 36.19*Foam::exp(-EtoP/1.2763e4) + 31.73*ones)
                                );
mupF = (1.0/press)*dimMu*(0.79*(0.05492*Foam::exp(-EtoP/6.858e3)  + 0.07509*Foam::exp(-EtoP/3.8175e4) + 0.0308*ones)+
                          0.21*(0.06841*Foam::exp(-EtoP/5.9678e4) + 0.09194*Foam::exp(-EtoP/1.2763e4) + 0.0320*ones)
                                );                               
etaF = press*dimAlpha*0.21*(1.307*ones + 3.32e4/EtoP*Foam::exp(-Foam::pow(Foam::log(EtoP)-9.04*ones,2)/2.53) );                            
// etaF = 0.1*alphaF;
       
DeF = constant::physicoChemical::k*Te*mueF/constant::electromagnetic::e;

DpF = constant::physicoChemical::k*Tgas*mupF/constant::electromagnetic::e;

DnF = constant::physicoChemical::k*Tgas*munF/constant::electromagnetic::e;


