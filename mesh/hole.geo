SetFactory("OpenCASCADE");

L = 5.;
H = 5.;
R = 1.;

M_PI = 3.14159265;

Rectangle(1) = {0., 0., 0., L/2., H/2., 0.};
Disk(2) = {0., 0., 0., R, R};
BooleanDifference{ Surface{1}; Delete; }{Surface{2}; Delete; }


Physical Curve("force_x", 1) = {4};
Physical Curve("fix_x", 2) = {2};
Physical Curve("fix_y", 3) = {5};
Physical Surface("hole", 1) = {1};

// Mesh settings
Mesh.MeshSizeFactor = 0.15;
Mesh.Algorithm = 6;
Mesh 2;
