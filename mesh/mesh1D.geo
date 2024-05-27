SetFactory("OpenCASCADE");
Point(1) = {-1, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Line(1) = {1, 2};
Transfinite Line {1} = 200 Using Progression 1;
Mesh.MshFileVersion = 2.2;

