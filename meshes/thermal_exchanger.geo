// Heat Exchanger mesh

SetFactory("OpenCASCADE");

// ---------------- PARAMETERS ----------------
n = 10;

length = 10;
h_insul = 0.5;
h_fluid = 1;
h_membrane = 0.1;
// number of divisions
nx = Ceil(n * length);  // along pipe
n_insul = Ceil(n * h_insul);
n_fluid = Ceil(n * h_fluid);
n_membrane = Ceil(n * h_membrane);

// ---------------- POINTS ----------------
Point(1) = {0, 0, 0};
Point(2) = {length, 0, 0};
Point(3) = {length, h_insul, 0};
Point(4) = {length, h_insul + h_fluid, 0};
Point(5) = {length, h_insul + h_fluid + h_membrane, 0};
Point(6) = {length, h_insul + 2 * h_fluid + h_membrane, 0};
Point(7) = {length, 2* h_insul + 2 * h_fluid + h_membrane, 0};
Point(8) = {0, 2* h_insul + 2 * h_fluid + h_membrane, 0};
Point(9) = {0, h_insul + 2 * h_fluid + h_membrane, 0};
Point(10) = {0, h_insul + h_fluid + h_membrane, 0};
Point(11) = {0, h_insul + h_fluid, 0};
Point(12) = {0, h_insul, 0};

// ---------------- LINES ----------------
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 1};

Line(13) = {12, 3};
Line(14) = {11, 4};
Line(15) = {10, 5};
Line(16) = {9, 6};

// ---------------- CURVES and SURFACES ----------------
Curve Loop(100) = {1, 2, -13, 12}; // bottom wall
Plane Surface(100) = {100};
Curve Loop(200) = {13, 3, -14, 11}; // bottom fluid
Plane Surface(200) = {200};
Curve Loop(300) = {14, 4, -15, 10}; // membrane
Plane Surface(300) = {300};
Curve Loop(400) = {15, 5, -16, 9}; // top fluid
Plane Surface(400) = {400};
Curve Loop(500) = {16, 6, 7, 8}; // top wall
Plane Surface(500) = {500};


// --- transfinite curves ---
Transfinite Curve {1, 13, 14, 15, 16, 7} = nx Using Progression 1;
Transfinite Curve {2, 12, 6, 8} = n_insul Using Progression 1;
Transfinite Curve {3, 11, 5, 9} = n_fluid Using Progression 1;
Transfinite Curve {4, 10} = n_membrane Using Progression 1;

// --- transfinite surfaces ---
Transfinite Surface {100};
Transfinite Surface {200};
Transfinite Surface {300};
Transfinite Surface {400};
Transfinite Surface {500};


// keep quads
Recombine Surface {100, 200, 300, 400, 500};

// ---------------- PHYSICAL GROUPS ----------------

// boundaries
Physical Curve(10) = {1}; // bottom
Physical Curve(20) = {7}; // top
Physical Curve(30) = {11}; // lower inlet
Physical Curve(40) = {5}; // upper inlet
// Physical Curve(70) = {9, 10}; // fluid - walls interfaces


// Materials
Physical Surface(0) = {200, 400}; // fluid
Physical Surface(1) = {100, 500, 300}; // walls
// Physical Surface(2) = {300}; // membrane

// ---------------- MESH ----------------
Mesh.RecombineAll = 1;
Mesh.RecombinationAlgorithm = 2;
Mesh.SaveAll = 0;
