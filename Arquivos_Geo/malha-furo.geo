//+ 1o. quadrado
lc1 = 0.05;
Point(1) = {-1, -0, 0, lc1};
Point(2) = {0 , -0, 0, lc1};
Point(3) = {-0, 1 , 0, lc1};
Point(4) = {-1, 1 , 0, lc1};
//+
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
//+ 2o. quadrado
lc2 = 0.1;
Point(5) = {-0.3, 0.7, 0, lc2};
Point(6) = {-0.7, 0.7, 0, lc2};
Point(7) = {-0.7, 0.3, 0, lc2};
Point(8) = {-0.3, 0.3, 0, lc2};
//+
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
//+
Curve Loop(1) = {3, 4, 1, 2};
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(1) = {1, 2};
//+
Physical Curve("externo") = {3, 2, 1, 4};
Physical Curve("interno") = {5, 8, 7, 6};
Physical Surface(3) = {1};
