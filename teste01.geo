
    lc = 0.0005;
    Point(1) = {0, 0, 0, lc};
    Point(2) = {0.01, 0, 0, lc};
    Point(3) = {0, 0.01, 0, lc};
    Point(4) = {-0.01, 0, 0, lc};
    Point(5) = {0, -0.01, 0, lc};
    Circle(1) = {2, 1, 3};
    Circle(2) = {3, 1, 4};
    Circle(3) = {4, 1, 5};
    Circle(4) = {5, 1, 2};
    Line Loop(5) = {1, 2, 3, 4};
    Plane Surface(6) = {5};
    Physical Curve(100) = {1,2,3,4};
    Physical Surface(200) = {6};
    Mesh 2;
    Save "teste01.msh";
    