SetFactory("OpenCASCADE");

// ---------------- PARAMETERS ----------------
pipe_diameter = 1;
pipe_length = 3;
walls_thickness = 1;

// number of divisions
nx = 15;  // along pipe
ny = 5;  // across thickness

// ---------------- POINTS ----------------
Point(1) = {0, 0, 0};
Point(2) = {pipe_length, 0, 0};
Point(3) = {pipe_length, walls_thickness, 0};
Point(4) = {pipe_length, walls_thickness + pipe_diameter, 0};
Point(5) = {pipe_length, 2 * walls_thickness + pipe_diameter, 0};
Point(6) = {0, 2 * walls_thickness + pipe_diameter, 0};
Point(7) = {0, walls_thickness + pipe_diameter, 0};
Point(8) = {0, walls_thickness, 0};

// ---------------- LINES ----------------
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};

Line(9) = {3, 8};
Line(10) = {4, 7};

// ---------------- CURVES and SURFACES ----------------
Curve Loop(100) = {1, 2, 9, 8}; // top wall
Plane Surface(100) = {100};
Curve Loop(200) = {-9, 3, 10, 7}; // fluid
Plane Surface(200) = {200};
Curve Loop(300) = {-10, 4, 5, 6}; // bottom wall
Plane Surface(300) = {300};

// --- transfinite curves ---
Transfinite Curve {1, 9, 10, 5} = nx Using Progression 1;
Transfinite Curve {2, 3, 4, 6, 7, 8} = ny Using Progression 1;

// --- transfinite surfaces ---
Transfinite Surface {100};
Transfinite Surface {200};
Transfinite Surface {300};

// keep quads
Recombine Surface {100, 200, 300};

// ---------------- PHYSICAL GROUPS ----------------

// boundaries
Physical Curve(10) = {1}; // top
Physical Curve(20) = {5}; // bottom
Physical Curve(30) = {2, 4}; // back of walls
Physical Curve(40) = {8, 6}; // front of walls
Physical Curve(50) = {3}; // outlet
Physical Curve(60) = {7}; // inlet
// Physical Curve(70) = {9, 10}; // fluid - walls interfaces


// Materials
Physical Surface(1) = {100, 300}; // walls
Physical Surface(2) = {200}; // fluid

// ---------------- MESH ----------------
Mesh.RecombineAll = 1;
Mesh.RecombinationAlgorithm = 2;
Mesh.SaveAll = 0;
