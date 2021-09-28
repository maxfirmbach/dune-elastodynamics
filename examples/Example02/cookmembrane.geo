
Point(1)  = {0.0, 0.0, 0.0};
Point(2)  = {48.0, 44.0, 0.0};
Point(3)  = {48.0, 60.0, 0.0};
Point(4)  = {0.0, 44.0, 0.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Transfinite Line{1, 3} = 33;
Transfinite Line{2, 4} = 33;

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Transfinite Surface{1};

Physical Surface(1) = {1};
Physical Line(0) = {1, 3};
Physical Line(1) = {4};
Physical Line(2) = {2};
