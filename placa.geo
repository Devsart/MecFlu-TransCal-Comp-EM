lc = 0.1;
lt = 0.01;
ls = 0.05;
//+
Point(1) = {-0.5, 0, 0, lc};
Point(2) = {0.5, 0, 0, lc};
Point(3) = {0.5, 1, 0, lc};
Point(4) = {-0.5, 1, 0, lc};
//+
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
//+ barramento
Point(5) = {-0.4, 0.8, 0, lt};
Point(6) = {-0.2, 0.8, 0, lt};
Point(7) = {-0.2, 0.2, 0, lt};
Point(8) = {-0.4, 0.2, 0, lt};
Point(9) = {-0.1, 0.8, 0, lt};
Point(10) = {0.1, 0.8, 0, lt};
Point(11) = {0.1, 0.2, 0, lt};
Point(12) = {-0.1, 0.2, 0, lt};
//+ cpus
Point(13) = {0.3, 0.8, 0, ls};
Point(14) = {0.3, 0.7, 0, ls};
Point(15) = {0.4, 0.7, 0, ls};
Point(16) = {0.4, 0.8, 0, ls};
Point(17) = {0.3, 0.5, 0, ls};
Point(18) = {0.4, 0.5, 0, ls};
Point(19) = {0.4, 0.4, 0, ls};
Point(20) = {0.3, 0.4, 0, ls};
Point(21) = {0.3, 0.2, 0, ls};
Point(22) = {0.4, 0.2, 0, ls};
Point(23) = {0.4, 0.1, 0, ls};
Point(24) = {0.3, 0.1, 0, ls};
//+
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 9};
//+
Line(13) = {13, 16};
Line(14) = {16, 15};
Line(15) = {15, 14};
Line(16) = {14, 13};
Line(17) = {17, 18};
Line(18) = {18, 19};
Line(19) = {19, 20};
Line(20) = {20, 17};
Line(21) = {21, 22};
Line(22) = {22, 23};
Line(23) = {23, 24};
Line(24) = {24, 21};
//+
Curve Loop(1) = {3, 4, 1, 2};
Curve Loop(2) = {5, 6, 7, 8};
Curve Loop(3) = {9, 10, 11, 12};
Curve Loop(4) = {13, 14, 15, 16};
Curve Loop(5) = {17, 18, 19, 20};
Curve Loop(6) = {21, 22, 23, 24};
//+
Plane Surface(1) = {1, 2, 3, 4, 5, 6};
//+
Plane Surface(2) = {2};
Plane Surface(3) = {3};
//+
Plane Surface(4) = {4};
//+
Plane Surface(5) = {5};
//+
Plane Surface(6) = {6};
