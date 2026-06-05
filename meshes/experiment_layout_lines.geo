SetFactory("OpenCASCADE");

// ---------------- PARAMETERS ----------------
system_width  = 15;
system_height = 10;

insulator_top = 2;
brass_top     = 4;

gap = 1;
blade_top = system_height - gap;

blade_tip_width = 0.2;  // small nonzero for top
blade_width     = 2;    // width at brass_top

W = system_width/2;
xb = blade_width/2;
xt = blade_tip_width/2;

// ---------------- POINTS ----------------
// Outer rectangle corners
Point(1) = {-W, 0, 0};
Point(2) = { W, 0, 0};
Point(3) = { W, system_height, 0};
Point(4) = {-W, system_height, 0};

// Insulator top (horizontal line)
Point(5)  = {-W, insulator_top, 0};
Point(6)  = { W, insulator_top, 0};

// Brass top (horizontal line)
Point(7)  = {-W, brass_top, 0};
Point(8)  = { W, brass_top, 0};

// Blade corners
Point(9)  = {-xb, brass_top, 0};   // left at brass_top
Point(10) = { xb, brass_top, 0};   // right at brass_top
Point(11) = {0, blade_top, 0};   // tip
Point(13) = {-xb, insulator_top, 0}; // left bottom of blade
Point(14) = { xb, insulator_top, 0}; // right bottom of blade

// ---------------- CURVES ----------------

// Outer rectangle
Line(1) = {1,2}; // bottom
Line(2) = {2,3}; // right
Line(3) = {3,4}; // up
Line(4) = {4,1}; // left
// Curve Loop(400) = {1,2,3,4};
// Plane Surface(400) = {400};

// Blade + insulator surface (one connected area)
Line(6) = {11,10};     // right blade side (sloped)
Line(7) = {10,14};    // right blade side (vertical)
Line(8) = {14,6};      // horizontal right to system border
Line(9) = {6,2};       // right side of insulator
Line(10) = {2,1};      // bottom
Line(11) = {1,5};    // left side of insulator
Line(12) = {5,13};   // system border to horizontal left
Line(13) = {13,9};   // vertical left
Line(14) = {9,11};    // left blade side (sloped)

// Curve loop and surface for insulator+blade
Curve Loop(100) = {6,7,8,9,10, 11, 12, 13, 14};
Plane Surface(100) = {100};

// Brass left
Line(15) = {7,9};      // top
Line(18) = {5,7};      // left side
Curve Loop(200) = {15, -13, -12, 18};
Plane Surface(200) = {200};

// Brass right
Line(19) = {10,8};     // top
Line(20) = {8,6};     // right side
Curve Loop(300) = {19, 20, -8, -7};
Plane Surface(300) = {300};

// Water
Line(23) = {4,3};
Line(24) = {3,8};
Line(30) = {7,4};
Curve Loop(500) = {23, 24, -19, -6, -14, -15, 30};
Plane Surface(500) = {500};


// ---------------- PHYSICAL GROUPS ----------------

// Boundaries
Physical Curve(10) = {30};              // left inlet, fluid segment only
Physical Curve(20) = {10, 11, 18, 20, 9}; // exterior solid boundary
Physical Curve(30) = {23};              // top
Physical Curve(40) = {24};              // right outlet, fluid segment only

// Materials
Physical Surface(1) = {100};         // insulator + blade
Physical Surface(2) = {200,300};     // brass left + right
//BooleanDifference(3) = { Surface{400}; Delete; }{ Surface{100,200,300};};
//Physical Surface(3) = {3};
Physical Surface(3) = {500};

