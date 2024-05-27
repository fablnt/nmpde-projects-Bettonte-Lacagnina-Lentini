// Gmsh project created on Mon May 13 20:23:08 2024
SetFactory("OpenCASCADE");
//+
Ellipse(1) = {0, 0, 0, 70, 50, 0, 2*Pi};
//+
Mesh.Algorithm = 6;
//+
Mesh.ElementOrder = 1;
//+
Curve Loop(1) = {1};
//+
Surface(1) = {1};
//+
Physical Surface(16) = {1};
//+
Physical Curve(17) = {1};

