// Gmsh project created on Thu Oct 28 02:25:34 2021
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 1, 0, 1.0};
//+
Point(2) = {0, 0.5, 0, 1.0};
//+
Point(3) = {0.75, 0.5, 0, 1.0};
//+
Point(4) = {0.75, 0, 0, 1.0};
//+
Point(5) = {3, 0, 0, 1.0};
//+
Point(6) = {3, 1, 0, 1.0};
//+
Line(1) = {2, 1};
//+
Line(2) = {1, 6};
//+
Line(3) = {6, 5};
//+
Line(4) = {5, 4};
//+
Line(5) = {4, 3};
//+
Line(6) = {3, 2};
//+
Line Loop(1) = {2, 3, 4, 5, 6, 1};
//+
Plane Surface(1) = {1};
