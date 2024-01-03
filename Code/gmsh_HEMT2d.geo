device_width =     9.0e-4;
electrode_width=   5.0e-5;
gate_width=        2.0e-4;

air_thickness =        1e-7;
device_thickness =     1.5e-4;
AlGaN_thickness=       2.5e-6;
GaN_thickness=         1.475e-4;

x_nitride_spacing=        2.5e-5;
x_gate_spacing=           1.5e-5;
y_channel_maxspacing=     1.0e-7;
y_channel_minspacing=     1.0e-8;
y_GaN_upspacing=          5.0e-6;
y_GaN_midspacing=         1.0e-5;
y_GaN_bottomspacing=     1.0e-4;


x_bulk_left=     0.0;
x_bulk_right=    x_bulk_left+device_width;
x_center=        0.5 * (x_bulk_left + x_bulk_right);
x_drain_right=   x_bulk_left+electrode_width;
x_source_left=   x_bulk_right-electrode_width;
x_gate_left =    x_center - 0.5 * (gate_width);
x_gate_right =   x_center + 0.5 * (gate_width);
x_device_left =  x_bulk_left - air_thickness;
x_device_right = x_bulk_right + air_thickness;



y_AlGaN_top=       0.0;
y_device_top=      y_AlGaN_top - air_thickness;
y_AlGaN_bottom=    y_AlGaN_top + AlGaN_thickness;
y_GaN_bottom=      y_AlGaN_bottom+GaN_thickness;
y_device_bottom=   y_GaN_bottom+air_thickness;


////AlGaN
Point(1) = {x_bulk_left, y_AlGaN_top, 0, 
x_nitride_spacing};
Point(2) = {x_drain_right, y_AlGaN_top, 0, 
x_nitride_spacing};
Point(3) = {x_gate_left, y_AlGaN_top, 0, 
x_gate_spacing};
Point(4) = {x_gate_right, y_AlGaN_top, 0, 
x_gate_spacing};
Point(5) = {x_source_left, y_AlGaN_top, 0, 
x_nitride_spacing};
Point(6) = {x_bulk_right, y_AlGaN_top, 0, 
x_nitride_spacing};
Point(7) = {x_bulk_right, y_AlGaN_bottom, 0, 
x_gate_spacing};
Point(8) = {x_bulk_left, y_AlGaN_bottom, 0, 
x_gate_spacing};
////GaN
Point(9) = {x_bulk_left, y_GaN_bottom, 0, 
x_gate_spacing};
Point(10) = {x_bulk_right, y_GaN_bottom, 0, 
x_gate_spacing};

Line(1) = {8, 7};
Line(2) = {7, 6};
Line(3) = {6, 5};
Line(4) = {5, 4};
Line(5) = {4, 3};
Line(6) = {3, 2};
Line(7) = {2, 1};
Line(8) = {1, 8};
Line Loop(9) = {1, 2, 3, 4, 5, 6,7,8};
Plane Surface(10) = {9};
Line(11) = {9, 10};
Line(12) = {10, 7};
Line(13) = {8, 9};
Line Loop(14) = {11,12,-1,13};
Plane Surface(15) = {14};

Field[1] = Attractor;
Field[1].NNodesByEdge = 500;
Field[1].EdgesList = {1};

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = y_channel_maxspacing;  // 在吸引子附近网格尺寸最小
Field[2].LcMax = y_GaN_bottomspacing;   // 远离吸引子的网格尺寸最大
Field[2].DistMin = 1e-7;     // 网格开始细化的距离
Field[2].DistMax = GaN_thickness*5;    // 网格恢复粗糙的距离


Field[8] = Min;
Field[8].FieldsList = {2};
Background Field = 8;


Physical Line("gate_contact") = {5};
Physical Line("AlGaN_GaN_interface") = {1};
Physical Line("source_contact") = {3};
Physical Line("drain_contact") = {7};
Physical Line("body_contact") = {11};
Physical Surface("AlGaN") = {10};
Physical Surface("GaN") = {15};
