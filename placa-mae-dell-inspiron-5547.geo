// Gmsh project created on Thu Oct 21 13:46:51 2021
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {27, 0, 0, 1.0};
//+
Point(3) = {27, 26.3, 0, 1.0};
//+
Point(4) = {0, 26.3, 0, 1.0};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Point(5) = {20, 20.5, 0, 1.0};
//+
Point(6) = {23.5, 20.5, 0, 1.0};
//+
Point(7) = {23.5, 24, 0, 1.0};
//+
Point(8) = {20, 24, 0, 1.0};
//+
Point(9) = {1.5, 21.5, 0, 1.0};
//+
Point(10) = {18.5, 21.5, 0, 1.0};
//+
Point(11) = {18.5, 22.5, 0, 1.0};
//+
Point(12) = {1.5, 22.5, 0, 1.0};
//+
Point(13) = {1.5, 23, 0, 1.0};
//+
Point(14) = {18.5, 23, 0, 1.0};
//+
Point(15) = {18.5, 24, 0, 1.0};
//+
Point(16) = {1.5, 24, 0, 1.0};
//+
Line(5) = {9, 10};
//+
Line(6) = {10, 11};
//+
Line(7) = {11, 12};
//+
Line(8) = {12, 9};
//+
Line(9) = {13, 14};
//+
Line(10) = {14, 15};
//+
Line(11) = {15, 16};
//+
Line(12) = {16, 13};
//+
Line(13) = {5, 6};
//+
Line(14) = {6, 7};
//+
Line(15) = {7, 8};
//+
Line(16) = {8, 5};
//+
Point(17) = {8.25, 12, 0, 1.0};
//+
Point(18) = {13.75, 12, 0, 1.0};
//+
Point(19) = {13.75, 17.5, 0, 1.0};
//+
Point(20) = {8.25, 17.5, 0, 1.0};
//+
Point(21) = {8.25, 18, 0, 1.0};
//+
Point(22) = {8.25, 19, 0, 1.0};
//+
Point(23) = {13.75, 19, 0, 1.0};
//+
Point(24) = {13.75, 18, 0, 1.0};
//+
Line(17) = {17, 20};
//+
Line(18) = {20, 19};
//+
Line(19) = {19, 18};
//+
Line(20) = {18, 17};
//+
Line(21) = {21, 22};
//+
Line(22) = {22, 23};
//+
Line(23) = {23, 24};
//+
Line(24) = {24, 21};
//+
Point(25) = {19, 5, 0, 1.0};
//+
Point(26) = {20, 5, 0, 1.0};
//+
Point(27) = {20, 17.75, 0, 1.0};
//+
Point(28) = {19, 17.75, 0, 1.0};
//+
Line(25) = {25, 28};
//+
Line(26) = {28, 27};
//+
Line(27) = {27, 26};
//+
Line(28) = {26, 25};
//+
Point(29) = {2.5, 0, 0, 1.0};
//+
Point(30) = {2.5, 2.5, 0, 1.0};
//+
Point(31) = {0, 2.5, 0, 1.0};
//+
Line(29) = {1, 31};
//+
Line(30) = {31, 30};
//+
Line(31) = {30, 29};
//+
Line(32) = {29, 1};
//+
Line(33) = {31, 4};
//+
Line(34) = {29, 2};
//+ Base-Placa
Line Loop(1) = {33, 30,31, 34, 2, 3};
Line Loop(2) = {11, 12, 9, 10};
Line Loop(3) = {7, 8, 5, 6};
Line Loop(4) = {22, 23, 24, 21};
Line Loop(5) = {18, 19, 20, 17};
Line Loop(6) = {25, 26, 27, 28};
Line Loop(7) = {16, 13, 14, 15};
Plane Surface(1) = {1, 2, 3, 4, 5, 6, 7};
//+ Alimentacao
Line Loop(9) = {30, 31, 32, 29};
Plane Surface(2) = {9};
//+ Processador AMD4
Line Loop(10) = {17, 18, 19, 20};
Plane Surface(3) = {10};
Line Loop(11) = {22, 23, 24, 21};
Plane Surface(4) = {11};
//+ SSD M2
Line Loop(12) = {25, 26, 27, 28};
Plane Surface(5) = {12};
//+ RAM DDR4 DIMM
Line Loop(13) = {7, 8, 5, 6};
Plane Surface(6) = {13};
Line Loop(14) = {11, 12, 9, 10};
Plane Surface(7) = {14};
//+ Placa de Video SATA 6Gbs
Line Loop(15) = {16, 13, 14, 15};
Plane Surface(8) = {15};
