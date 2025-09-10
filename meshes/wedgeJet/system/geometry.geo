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
Point(2) = {0, 0, L_ele, h};
Point(3) = {0, 0, L, h};
Point(4) = {0, R_ele, 0, h};
Point(5) = {0, R_ele, L_ele, h};
Point(6) = {0, R_dielIn, 0, h};
Point(7) = {0, R_dielIn, L_diel, h};
Point(8) = {0, R_dielOut, 0, h};
Point(9) = {0, R_dielOut, L_diel, h};
Point(10) = {0, R_farfield, 0, h};
Point(11) = {0, R_farfield, L, h};
Point(12) = {0, R_farfield, L_ele, h};

Line(1) = {4, 5};
Line(2) = {5, 2};
Line(3) = {2, 3};
Line(4) = {3, 11};
Line(5) = {12, 10};
Line(6) = {10, 8};
Line(7) = {8, 6};
Line(8) = {6, 4};
Line(9) = {5, 12};
Line(10) = {12, 11};

Curve Loop(1) = {1, 8, 7, 6, 5, 9};
Plane Surface(1) = {1};

Curve Loop(2) = {2, 9, 10, 4, 3};
Plane Surface(2) = {2};

Transfinite Line {1, 2, 9, 10} = 20 Using Progression 1;
Transfinite Surface {1};
Recombine Surface {1};

// For the outer layer, either:
Transfinite Line {5, 6, 7, 8} = 20 Using Progression 1;
Transfinite Surface {2};
Recombine Surface {2};

Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;

// Physical Volume("region0") = {1:4};
