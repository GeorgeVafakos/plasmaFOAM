// Here are the rate coefficients

scalar Te = 2*11606;

rate[0] = 2.584e-18*Foam::pow(Te,0.68)*Foam::exp(-285409/Te);
rate[1] = 2.308e-16*Foam::pow(Te,0.31)*Foam::exp(-229700/Te);
rate[2] = 4.661e-18*Foam::pow(Te,0.60)*Foam::exp(-55460/Te);
rate[3] = 2.7e-16;
rate[4] = 1.099e-17*Foam::pow(Te,0.31);
rate[5] = 1.268e-18*Foam::pow(Te,0.71)*Foam::exp(-39450/Te);
rate[6] = 5.386e-13*Foam::pow(Te,-0.5);
rate[7] = 1.3e-45;
rate[8] = 1.0e-43;

Ri_[composition.species().find("e")] = 
        + 1.0*dimRate*rate[1]*n[composition.species().find("e")] 
        + 1.0*dimRate*rate[18]*n[composition.species().find("Hem")]*n[composition.species().find("Hem")] 
        + 1.0*dimRate*rate[29]*n[composition.species().find("He2m")]*n[composition.species().find("N2")] 
        + 1.0*dimRate*rate[27]*n[composition.species().find("Hem")]*n[composition.species().find("N2")] 
        + 1.0*dimRate*rate[30]*n[composition.species().find("He2m")]*n[composition.species().find("O2")] 
        + 1.0*dimRate*rate[28]*n[composition.species().find("Hem")]*n[composition.species().find("O2")] 
        - 1.0*dimRate*rate[12]*n[composition.species().find("He2plus")]*n[composition.species().find("e")] 
        - 1.0*dimRate*rate[14]*n[composition.species().find("He2plus")]*n[composition.species().find("e")] 
        + 1.0*dimRate*rate[10]*n[composition.species().find("Hem")]*n[composition.species().find("e")] 
        - 1.0*dimRate*rate[11]*n[composition.species().find("Heplus")]*n[composition.species().find("e")] 
        - 1.0*dimRate*rate[13]*n[composition.species().find("Heplus")]*n[composition.species().find("e")] 
        + 1.0*dimRate*rate[4]*n[composition.species().find("N2")]*n[composition.species().find("e")] 
        - 1.0*dimRate*rate[5]*n[composition.species().find("O2")]*n[composition.species().find("e")] 
        - 1.0*dimRate*rate[8]*n[composition.species().find("O2")]*n[composition.species().find("e")] 
        + 1.0*dimRate*rate[9]*n[composition.species().find("O2")]*n[composition.species().find("e")] 
        - 1.0*dimRate*rate[17]*n[composition.species().find("O2plus")]*n[composition.species().find("e")] 
        - 1.0*dimRate2*rate[16]*n[composition.species().find("N2plus")]*n[composition.species().find("e")]*n[composition.species().find("e")] 
        - 1.0*dimRate2*rate[15]*n[composition.species().find("N2")]*n[composition.species().find("N2plus")]*n[composition.species().find("e")];

Ri_[composition.species().find("Heplus")] = 
        + 1.0*dimRate*rate[1]*n[composition.species().find("e")] 
        - 1.0*dimRate*rate[20]*n[composition.species().find("Heplus")] 
        + 1.0*dimRate*rate[18]*n[composition.species().find("Hem")]*n[composition.species().find("Hem")] 
        - 1.0*dimRate*rate[25]*n[composition.species().find("Heplus")]*n[composition.species().find("N2")] 
        + 1.0*dimRate*rate[10]*n[composition.species().find("Hem")]*n[composition.species().find("e")] 
        - 1.0*dimRate*rate[11]*n[composition.species().find("Heplus")]*n[composition.species().find("e")] 
        - 1.0*dimRate*rate[13]*n[composition.species().find("Heplus")]*n[composition.species().find("e")];

Ri_[composition.species().find("Hem")] = 
        + 1.0*dimRate*rate[0]*n[composition.species().find("e")] 
        - 1.0*dimRate*rate[19]*n[composition.species().find("Hem")] 
        - 2*dimRate*rate[18]*n[composition.species().find("Hem")]*n[composition.species().find("Hem")] 
        - 1.0*dimRate*rate[27]*n[composition.species().find("Hem")]*n[composition.species().find("N2")] 
        - 1.0*dimRate*rate[28]*n[composition.species().find("Hem")]*n[composition.species().find("O2")] 
        + 1.0*dimRate*rate[12]*n[composition.species().find("He2plus")]*n[composition.species().find("e")] 
        + 1.0*dimRate*rate[14]*n[composition.species().find("He2plus")]*n[composition.species().find("e")] 
        - 1.0*dimRate*rate[10]*n[composition.species().find("Hem")]*n[composition.species().find("e")] 
        + 1.0*dimRate*rate[11]*n[composition.species().find("Heplus")]*n[composition.species().find("e")] 
        + 1.0*dimRate*rate[13]*n[composition.species().find("Heplus")]*n[composition.species().find("e")];

Ri_[composition.species().find("He2plus")] = 
        + 1.0*dimRate*rate[20]*n[composition.species().find("Heplus")] 
        - 1.0*dimRate*rate[26]*n[composition.species().find("He2plus")]*n[composition.species().find("N2")] 
        - 1.0*dimRate*rate[12]*n[composition.species().find("He2plus")]*n[composition.species().find("e")] 
        - 1.0*dimRate*rate[14]*n[composition.species().find("He2plus")]*n[composition.species().find("e")];

Ri_[composition.species().find("He2m")] = 
        + 1.0*dimRate*rate[19]*n[composition.species().find("Hem")] 
        - 1.0*dimRate*rate[29]*n[composition.species().find("He2m")]*n[composition.species().find("N2")] 
        - 1.0*dimRate*rate[30]*n[composition.species().find("He2m")]*n[composition.species().find("O2")];

Ri_[composition.species().find("N2")] = 
        + 1.0*dimRate*rate[31]*n[composition.species().find("N2C3p")] 
        - 1.0*dimRate*rate[29]*n[composition.species().find("He2m")]*n[composition.species().find("N2")] 
        - 1.0*dimRate*rate[26]*n[composition.species().find("He2plus")]*n[composition.species().find("N2")] 
        - 1.0*dimRate*rate[27]*n[composition.species().find("Hem")]*n[composition.species().find("N2")] 
        - 1.0*dimRate*rate[25]*n[composition.species().find("Heplus")]*n[composition.species().find("N2")] 
        - 1.0*dimRate*rate[2]*n[composition.species().find("N2")]*n[composition.species().find("e")] 
        - 1.0*dimRate*rate[3]*n[composition.species().find("N2")]*n[composition.species().find("e")] 
        - 1.0*dimRate*rate[4]*n[composition.species().find("N2")]*n[composition.species().find("e")] 
        + 1.0*dimRate2*rate[16]*n[composition.species().find("N2plus")]*n[composition.species().find("e")]*n[composition.species().find("e")] 
        + 1.0*dimRate2*rate[15]*n[composition.species().find("N2")]*n[composition.species().find("N2plus")]*n[composition.species().find("e")];

Ri_[composition.species().find("N2plus")] = 
        + 1.0*dimRate*rate[29]*n[composition.species().find("He2m")]*n[composition.species().find("N2")] 
        + 1.0*dimRate*rate[26]*n[composition.species().find("He2plus")]*n[composition.species().find("N2")] 
        + 1.0*dimRate*rate[27]*n[composition.species().find("Hem")]*n[composition.species().find("N2")] 
        + 1.0*dimRate*rate[25]*n[composition.species().find("Heplus")]*n[composition.species().find("N2")] 
        + 1.0*dimRate*rate[4]*n[composition.species().find("N2")]*n[composition.species().find("e")] 
        - 1.0*dimRate2*rate[16]*n[composition.species().find("N2plus")]*n[composition.species().find("e")]*n[composition.species().find("e")] 
        - 1.0*dimRate2*rate[15]*n[composition.species().find("N2")]*n[composition.species().find("N2plus")]*n[composition.species().find("e")];

Ri_[composition.species().find("N2C3p")] = 
        + 1.0*dimRate*rate[2]*n[composition.species().find("N2")]*n[composition.species().find("e")] 
        - 1.0*dimRate*rate[31]*n[composition.species().find("N2C3p")];

Ri_[composition.species().find("NA3Su")] = 
        + 1.0*dimRate*rate[3]*n[composition.species().find("N2")]*n[composition.species().find("e")];

Ri_[composition.species().find("O2")] = 
        + 2*dimRate*rate[21]*n[composition.species().find("O2minus")]*n[composition.species().find("O2plus")] 
        - 1.0*dimRate*rate[28]*n[composition.species().find("Hem")]*n[composition.species().find("O2")] 
        - 1.0*dimRate*rate[30]*n[composition.species().find("He2m")]*n[composition.species().find("O2")] 
        + 2*dimRate*rate[22]*n[composition.species().find("O2minus")]*n[composition.species().find("O2plus")] 
        - 1.0*dimRate*rate[5]*n[composition.species().find("O2")]*n[composition.species().find("e")] 
        - 1.0*dimRate*rate[6]*n[composition.species().find("O2")]*n[composition.species().find("e")] 
        - 1.0*dimRate*rate[7]*n[composition.species().find("O2")]*n[composition.species().find("e")] 
        - 1.0*dimRate*rate[8]*n[composition.species().find("O2")]*n[composition.species().find("e")] 
        - 1.0*dimRate*rate[9]*n[composition.species().find("O2")]*n[composition.species().find("e")] 
        + 2*dimRate2*rate[23]*n[composition.species().find("N2")]*n[composition.species().find("O2minus")]*n[composition.species().find("O2plus")] 
        + 2*dimRate2*rate[24]*n[composition.species().find("O2")]*n[composition.species().find("O2minus")]*n[composition.species().find("O2plus")];

Ri_[composition.species().find("O2minus")] = 
        + 1.0*dimRate*rate[8]*n[composition.species().find("O2")]*n[composition.species().find("e")] 
        - 1.0*dimRate*rate[22]*n[composition.species().find("O2minus")]*n[composition.species().find("O2plus")] 
        - 1.0*dimRate*rate[21]*n[composition.species().find("O2minus")]*n[composition.species().find("O2plus")] 
        - 1.0*dimRate2*rate[23]*n[composition.species().find("N2")]*n[composition.species().find("O2minus")]*n[composition.species().find("O2plus")] 
        - 1.0*dimRate2*rate[24]*n[composition.species().find("O2")]*n[composition.species().find("O2minus")]*n[composition.species().find("O2plus")];

Ri_[composition.species().find("O2plus")] = 
        + 1.0*dimRate*rate[30]*n[composition.species().find("He2m")]*n[composition.species().find("O2")] 
        + 1.0*dimRate*rate[28]*n[composition.species().find("Hem")]*n[composition.species().find("O2")] 
        - 1.0*dimRate*rate[21]*n[composition.species().find("O2minus")]*n[composition.species().find("O2plus")] 
        - 1.0*dimRate*rate[22]*n[composition.species().find("O2minus")]*n[composition.species().find("O2plus")] 
        + 1.0*dimRate*rate[9]*n[composition.species().find("O2")]*n[composition.species().find("e")] 
        - 1.0*dimRate*rate[17]*n[composition.species().find("O2plus")]*n[composition.species().find("e")] 
        - 1.0*dimRate2*rate[23]*n[composition.species().find("N2")]*n[composition.species().find("O2minus")]*n[composition.species().find("O2plus")] 
        - 1.0*dimRate2*rate[24]*n[composition.species().find("O2")]*n[composition.species().find("O2minus")]*n[composition.species().find("O2plus")];

Ri_[composition.species().find("O2a1D")] = 
        + 1.0*dimRate*rate[5]*n[composition.species().find("O2")]*n[composition.species().find("e")];

Ri_[composition.species().find("O")] = 
        + 2*dimRate*rate[6]*n[composition.species().find("O2")]*n[composition.species().find("e")] 
        + 1.0*dimRate*rate[7]*n[composition.species().find("O2")]*n[composition.species().find("e")] 
        + 2*dimRate*rate[17]*n[composition.species().find("O2plus")]*n[composition.species().find("e")];

Ri_[composition.species().find("O1D")] = 
        + 1.0*dimRate*rate[7]*n[composition.species().find("O2")]*n[composition.species().find("e")];