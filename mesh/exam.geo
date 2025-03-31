SetFactory("OpenCASCADE");

cl = 1.;
L = 2.;
H = 1.;

Rectangle(1) = {0., 0., 0., L, H, 0.};

Physical Curve("fix_x", 1) = {1, 3, 4};
Physical Curve("fix_y", 2) = {1, 3};

// Mesh settings
eps = 1e-5;
p() = Point In BoundingBox{L-eps, H-eps, -eps, L+eps, H+eps, eps};
MeshSize { PointsOf{ Surface{1}; }} = 1.0;
MeshSize {p()} = 0.1;

Mesh.MeshSizeFactor = 0.10;
Mesh.Algorithm = 11;
Mesh 2;
