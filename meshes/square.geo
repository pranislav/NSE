SetFactory("OpenCASCADE");

// ---------------- PARAMETERS ----------------
a = 1;

// number of divisions
n = 5;

// ---------------- POINTS ----------------
Point(1) = {0, 0, 0};
Point(2) = {a, 0, 0};
Point(3) = {a, a, 0};
Point(4) = {0, a, 0};

// ---------------- LINES ----------------
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// ---------------- SURFACE ----------------
Curve Loop(10) = {1, 2, 3, 4};
Plane Surface(100) = {10};

// ---------------- MESH ----------------
Transfinite Curve {1, 2, 3, 4} = n + 1 Using Progression 1;
Transfinite Surface {100};

// ---------------- Boundaries ----------------
Physical Curve(10) = {1};
Physical Curve(20) = {2};
Physical Curve(30) = {3};
Physical Curve(40) = {4};

Physical Surface(100) = {100};

// ---------------- MESH ----------------
Mesh.RecombineAll = 1;
Mesh.RecombinationAlgorithm = 2;
Mesh.SaveAll = 0;