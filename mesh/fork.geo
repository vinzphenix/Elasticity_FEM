SetFactory("OpenCASCADE");

// Parameters
w = 0.03;
hl = 0.0956;
hr = 0.0956;
hb = 0.040;
dt = 0.008;
db = 0.007;
rr = 0.0015;

r = w / 2;
R = r + dt;
hmin = Min(hl, hr);

delta = Sqrt((R + rr)^2 - (db / 2. + rr)^2) - R;

// Bottom cut tools
For i In {0:1}
    sign = (1 - 2*i);
    tx = sign * db / 2. - R * i;
    ty = -R - hb;
    obj = 2+i;
    Rectangle(obj) = {tx, ty, 0, R, hb - delta};
    If (0. < rr)
        tool = 4;
        tx = sign * (db / 2. + rr);
        Disk(tool) = {tx, -R - delta, 0, rr, rr};
        BooleanUnion{Surface{obj}; Delete; }{Surface{tool}; Delete; }
    EndIf
    tx = sign * (db / 2. + rr) * R / (R + rr) - R * i;
    tool = 5;
    Rectangle(tool) = {tx, ty, 0., R, hb + R};
    BooleanUnion{Surface{obj}; Delete; } {Surface{tool}; Delete ;}
    Disk(tool) = {0., 0., 0., R, R};
    BooleanDifference{Surface{obj}; Delete; } {Surface{tool}; Delete; }
    tools~{i+1} = obj;
EndFor

obj = 1;
Rectangle(obj) = {
    -dt - w/2., -r - db - hb/2., 0., 
    2*dt+w, hb / 2. + db + r + hmin / 2.
};
BooleanDifference
    {Surface{obj}; Delete; }
    {Surface{tools~{1}}; Surface{tools~{2}}; Delete; }

// Full bottom branch
tool = 2;
Rectangle(tool) = {-db/2., -R-hb, 0., db, hb, rr};
BooleanUnion{Surface{obj}; Delete; } {Surface{tool}; Delete ;}

// Cut top branches
tool = 2;
Rectangle(tool) = {-w/2., 0, 0, w, hmin};
BooleanDifference {Surface{obj}; Delete; } {Surface{tool}; Delete; }
Disk(tool) = {0, 0, 0, r, r};
BooleanDifference {Surface{obj}; Delete; } {Surface{tool}; Delete; }

// Add top branches
tool = 2;
Rectangle(tool) = {-w/2. - dt, 0, 0, dt, hl, rr};
BooleanUnion{ Surface{obj}; Delete; } {Surface{tool}; Delete;}
Rectangle(tool) = {w/2, 0, 0, dt, hr, rr};
BooleanUnion{ Surface{obj}; Delete; } {Surface{tool}; Delete;}


// Size field
Field[1] = Distance;
If (0. < rr)
    Field[1].CurvesList = {60, 51, 54};
    Field[1].NodesList = {53, 54, 55, 60};
Else
    Field[1].CurvesList = {34, 38, 30};
    Field[1].NodesList = {31, 32, 35, 40};
EndIf
Field[1].Sampling = 50;
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = dt/8.;
Field[2].LcMax = dt/3.;
Field[2].DistMin = dt/8.;
Field[2].DistMax = dt/4.;

Field[3] = Distance;
If (0. < rr)
    Field[3].CurvesList = {44, 50};
Else
    Field[3].NodesList = {36, 39};
EndIf
Field[3].Sampling = 50;
Field[4] = Threshold;
Field[4].IField = 3;
Field[4].LcMin = dt/8.;
Field[4].LcMax = dt/3.;
Field[4].DistMin = dt/4.;
Field[4].DistMax = dt/2.;

Field[5] = Min;
Field[5].FieldsList = {2, 4};

Background Field = 5;

// Mesh settings
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFactor = 0.50;
Mesh.Algorithm = 5;
//Mesh.RecombinationAlgorithm = 3;
//Mesh.RecombineAll = 1;

Mesh 2;

// Physical groups
If (0. < rr)
    Physical Curve ("fix_n", 1) = {45, 46, 47, 48, 49};
    Physical Curve ("fix_t", 2) = {45, 46, 47, 48, 49};
    Physical Curve ("force_n", 3) = {53, 55};
Else
    Physical Curve ("fix_n", 1) = {35, 36, 37};
    Physical Curve ("fix_t", 2) = {35, 36, 37};
    Physical Curve ("force_n", 3) = {31, 40};
EndIf

Physical Surface("fork") = {1};
