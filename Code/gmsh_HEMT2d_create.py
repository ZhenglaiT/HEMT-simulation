# Copyright 2013 Devsim LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from devsim import *
device ="HEMT"

device_width    =9.0e-4
electrode_width =5.0e-5
gate_width      =2.0e-4


air_thickness    =1e-7
device_thickness =1.5e-4
AlGaN_thickness  =2.5e-6
GaN_thickness    =1.475e-4
#这个是参考Silvaco的相关网格参数，而由于Gmsh的不同，所以实际的离散和这个很不同
x_nitride_spacing     =2.5e-5
x_gate_spacing        =1.5e-5
y_channel_maxspacing  =1.0e-7
y_channel_minspacing  =1.0e-8
y_GaN_upspacing       =5.0e-6
y_GaN_midspacing      =1.0e-5
y_GaN_bottomspacing   =1.0e-4


# mesh generation
# gmsh -2 -format msh2 gmsh_HEMT2d.geo -o gmsh_HEMT2d.msh

create_gmsh_mesh(file="gmsh_HEMT2d.msh", mesh="HEMT")
add_gmsh_region    (mesh="HEMT", gmsh_name="AlGaN",    region="AlGaN", material="AlGaN")
add_gmsh_region    (mesh="HEMT", gmsh_name="GaN",    region="GaN", material="GaN")


add_gmsh_contact   (mesh="HEMT", gmsh_name="drain_contact",  region="AlGaN", name="drain", material="metal")
add_gmsh_contact   (mesh="HEMT", gmsh_name="source_contact", region="AlGaN", name="source", material="metal")
add_gmsh_contact   (mesh="HEMT", gmsh_name="body_contact",   region="GaN", name="body", material="metal")
add_gmsh_contact   (mesh="HEMT", gmsh_name="gate_contact",   region="AlGaN", name="gate", material="metal")
# add_gmsh_contact   (mesh="SOI", gmsh_name="Thermalcontact_right",   region="bulk", name="Thermalright", material="metal")
# add_gmsh_contact   (mesh="SOI", gmsh_name="Thermalcontact_left",   region="bulk", name="Thermalleft", material="metal")



add_gmsh_interface (mesh="HEMT", gmsh_name="AlGaN_GaN_interface", region0="AlGaN", region1="GaN", name="AlGaN_GaN")

finalize_mesh(mesh="HEMT")
create_device(mesh="HEMT", device="HEMT")

node_model(name="Donors",    device=device, region="AlGaN", equation="1.001")
node_model(name="Acceptors", device=device, region="AlGaN", equation="1")
node_model(name="NetDoping",    device=device, region="AlGaN", equation="Donors - Acceptors")

node_model(name="Donors",    device=device, region="GaN", equation="1.001")
node_model(name="Acceptors", device=device, region="GaN", equation="1")
node_model(name="NetDoping",    device=device, region="GaN", equation="Donors - Acceptors")
#the little doping is to avoid singularity in the calculation of mobility