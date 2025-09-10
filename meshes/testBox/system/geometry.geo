// Gmsh geometry for the plasma jet reactor
R_farfield = 5e-3;
R_dielOut = 3e-3;
R_dielIn = 2e-3;
R_ele = 0.5e-3;
L_diel = 10e-3;
L_ele = 3e-3;
L = 20e-3;
h=0.001;

Point(1) = {0, 0, 0, h};
Point(2) = {0, 0, 1, h};
Point(3) = {0, 1, 1, h};
Point(4) = {0, 1, 0, 0.001*h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Line {1} = 10 Using Progression 1;
Transfinite Line {2} = 10 Using Progression 1;
Transfinite Line {3} = 10 Using Progression 1;
// Transfinite Line {4} = 20 Using Progression 1;
Transfinite Surface {1};
Recombine Surface {1};

