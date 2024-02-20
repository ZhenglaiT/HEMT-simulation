device_width =     2.5e-4;
electrode_width=   2.5e-5;
gate_width=        1.0e-5;

air_thickness =        1e-7;
device_thickness =     1.02e-4;
electrode_thickness=   5e-6;
AlGaN_thickness=       2e-6;
GaN_thickness=         1e-4;

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



y_electrode_top=       0.0;
y_device_top=    y_electrode_top -electrode_thickness-air_thickness;
y_AlGaN_top=    y_electrode_top -electrode_thickness;
y_AlGaN_bottom=   y_AlGaN_top + AlGaN_thickness;
y_GaN_bottom=      y_electrode_top+GaN_thickness-electrode_thickness;
y_device_bottom=   y_GaN_bottom+air_thickness;


////AlGaN
Point(1) = {x_bulk_left, y_electrode_top, 0, 
x_nitride_spacing};
Point(2) = {x_drain_right, y_electrode_top, 0, 
x_nitride_spacing};
Point(3) = {x_drain_right, y_AlGaN_bottom, 0, 
x_gate_spacing};
Point(4) = {x_drain_right, y_AlGaN_top, 0, 
x_gate_spacing};
Point(5) = {x_gate_left, y_AlGaN_top, 0, 
x_nitride_spacing};
Point(6) = {x_gate_right, y_AlGaN_top, 0, 
x_nitride_spacing};
Point(7) = {x_source_left, y_AlGaN_top, 0, 
x_gate_spacing};
Point(8) = {x_source_left, y_AlGaN_bottom, 0, 
x_gate_spacing};
////GaN
Point(9) = {x_source_left,  y_electrode_top, 0, 
x_gate_spacing};
Point(10) = {x_bulk_right, y_electrode_top, 0, 
x_gate_spacing};
Point(11) = {x_bulk_right, y_GaN_bottom, 0, 
x_gate_spacing};
Point(12) = {x_bulk_left, y_GaN_bottom, 0, 
x_gate_spacing};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 8};
Line(4) = {8, 9};
Line(5) = {9, 10};
Line(6) = {10, 11};
Line(7) = {11, 12};
Line(8) = {12, 1};
Line Loop(9) = {1, 2, 3, 4, 5, 6,7,8};
Plane Surface(10) = {9};
Line(11) = {3, 4};
Line(12) = {4, 5};
Line(13) = {5, 6};
Line(14) = {6, 7};
Line(15) = {7, 8};
Line Loop(16) = {11,12,13,14,15,-3};
Plane Surface(17) = {16};

Field[1] = Attractor;
Field[1].NNodesByEdge = 1000;
Field[1].EdgesList = {3};

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = y_channel_maxspacing/2;  // 在吸引子附近网格尺寸最小
Field[2].LcMax = y_GaN_bottomspacing/2;   // 远离吸引子的网格尺寸最大
Field[2].DistMin = 1e-7;     // 网格开始细化的距离
Field[2].DistMax = GaN_thickness*5;    // 网格恢复粗糙的距离


Field[8] = Min;
Field[8].FieldsList = {2};
Background Field = 8;


Physical Line("gate_contact") = {13};
//Physical Line("AlGaN_GaN_interface") = {3};
//Physical Line("source_contact_AlGaN") = {15};
Physical Line("source_contact") = {4,5};
//Physical Line("drain_contact_AlGaN") = {11};
Physical Line("drain_contact") = {1,2};
Physical Line("body_contact") = {7};
Physical Surface("AlGaN") = {17};
Physical Surface("GaN") = {10};
