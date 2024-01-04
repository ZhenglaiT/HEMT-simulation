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
from simple_physics import *
from devsim.python_packages.ramp import *

import sys
#sys.path.append("C:\programdata\anaconda3\lib\site-packages")
import numpy
import gmsh_HEMT2d_create

if False:
    set_parameter(name = "extended_solver", value=True)
    set_parameter(name = "extended_model", value=True)
    set_parameter(name = "extended_equation", value=True)
set_parameter(name="threads_available", value=4)
set_parameter(name="V_t", value=0.0259)


device = "HEMT"
AlGaN_regions = ("AlGaN",)
GaN_regions = ("GaN",)
regions = ("AlGaN", "GaN")
interfaces = ("AlGaN_GaN",)
SchottkyContacts = ("gate","drain","source")
OhmicContacts =("body",)


set_parameter(device=device, name="initialTemperature", value=300.0)

for i in AlGaN_regions:
    SetAlGaNParameters(device, i, xcomp=0.3)

for i in GaN_regions:
    SetGaNParameters(device, i)

for i in regions:
    CreateSolution(device, i, "Potential")
    CreateSolution(device, i, "Temperature")
    CreateSolution(device, i, "Temperature:Electrons")
    CreateSolution(device, i, "Temperature:Holes")
    CreateTinitial(device,i)
    set_node_values(device=device, region=i, name="Temperature", init_from="Init_temperature")
    CreateDensityofStates(device,i)
    edge_from_node_model(device=device, region=i, node_model="Temperature")
    edge_average_model(device=device, region=i,edge_model="edgeTemp", node_model="Temperature")
    edge_average_model(device=device, region=i, edge_model="edgeTemp", node_model="Temperature", derivative="Temperature")
    # initial the temperature field for equilibrium solution

for i in regions:
    V_t_T = "boltzmannConstant*Temperature/(ElectronCharge)"
    node_model(device=device, region=i, name="V_t_T", equation=V_t_T)
    node_model(device=device, region=i, name="V_t_T:Temperature", equation="boltzmannConstant/ElectronCharge")
    edge_average_model(device=device, region=i,edge_model="edgeV_t_T", node_model="V_t_T")
    edge_average_model(device=device, region=i, edge_model="edgeV_t_T",node_model="V_t_T", derivative="Temperature")
# create model to calculate mobility and thermal voltage

for i in AlGaN_regions:
    CreatAlGaNPotentialOnly(device, i)


for i in GaN_regions:
    CreatGaNPotentialOnly(device, i)

# Set up contacts
for i in ("gate",):
    tmp = get_region_list(device=device, contact=i)
    r = tmp[0]
    print("%s %s" % (r, i))
    CreateSchottkyContacts(device, r, i, workfun=5.0)
    set_parameter(device=device, name="gate_bias", value=0.0)

for i in ("drain", "source"):
    tmp = get_region_list(device=device, contact=i)
    r = tmp[0]
    print("%s %s" % (r, i))
    CreateSchottkyContacts(device, r, i, workfun=3.93)
    set_parameter(device=device, name="drain_bias", value=0.0)
    set_parameter(device=device, name="source_bias", value=0.0)

for i in ("body",):
    tmp = get_region_list(device=device, contact=i)
    r = tmp[0]
    print("%s %s" % (r, i))
    CreateSiliconPotentialOnlyContact(device, r, i)
    set_parameter(device=device, name="body_bias", value=0.0)

for i in interfaces:
    CreateHeterojunctionInterface(device, i)

solve(type="dc", absolute_error=1.0e-13,relative_error=1e-12, maximum_iterations=30)
solve(type="dc", absolute_error=1.0e-13,relative_error=1e-12, maximum_iterations=30)

write_devices(file="gmsh_HEMT2d_potentialonly", type="tecplot")

# T field

for i in regions:
    CreateSolution(device, i, "Electrons")
    CreateSolution(device, i, "Holes")
    set_node_values(device=device, region=i, name="Electrons",init_from="IntrinsicElectrons")
    set_node_values(device=device, region=i, name="Holes",init_from="IntrinsicHoles")
    # TODO

    CreateAlbrechtMobilityModels(device, i)
   
  
CreateAlGaNDriftDiffusion(device, "AlGaN",mu_n="mu_n_Albrecht", mu_p="mu_p_Albrecht")
CreateGaNDriftDiffusion(device, "GaN",mu_n="mu_n_Albrecht", mu_p="mu_p_Albrecht")


for i in interfaces:
    CreateThermionicEmission(device, i)


# for i in regions:
#     # TODO: create derivative at the same time as the temperature
#     # TODO: maybe edge_average_model doesn't need intermediate derivatives?
#     if i in silicon_regions:
#         CreateJouleHeatSilicon(device, i)
#         CreateSiThermalConductivity(device,i)
#     elif i in oxide_regions:
#         CreateJouleHeatOxide(device, i)
#         CreateOxideThermalConductivity(device, i)           
#     else:
#         raise RuntimeError("unknown region type for region " + i)

# CreateTemperaturefield(device,i)

for c in ("gate",):
    tmp = get_region_list(device=device, contact=c)
    r = tmp[0]
    CreateDriftDiffusionAtSchottkyContact(device, r, c, workfun=5.0)

for c in ("drain", "source"):
    tmp = get_region_list(device=device, contact=c)
    r = tmp[0]
    CreateDriftDiffusionAtSchottkyContact(device, r, c, workfun=3.93)

for c in OhmicContacts:
    tmp = get_region_list(device=device, contact=c)
    r = tmp[0]
    CreateSiliconDriftDiffusionAtContact(device, r, c)


solve(type="dc", absolute_error=1.0e30,relative_error=1e-5, maximum_iterations=1000)


for r in regions:
    node_model(device=device, region=r, name="logElectrons",equation="log(Electrons)/log(10)")



for r in regions:
    element_from_edge_model(edge_model="ElectricField",   device=device, region=r)
    element_from_edge_model(edge_model="ElectronCurrent", device=device, region=r)
    element_from_edge_model(edge_model="HoleCurrent",     device=device, region=r)



element_model(device=device, region="AlGaN", name="Emag", equation="(ElectricField_x^2 + ElectricField_y^2)^(0.5)")
element_model(device=device, region="AlGaN", name="mumag", equation="(mu_vsat_e_x^2 + mu_vsat_e_y^2)^(0.5)")
element_model(device=device, region="GaN", name="Emag", equation="(ElectricField_x^2 + ElectricField_y^2)^(0.5)")
element_model(device=device, region="GaN", name="mumag", equation="(mu_vsat_e_x^2 + mu_vsat_e_y^2)^(0.5)")

#
printAllCurrents(device)
rampbias(device, "gate", 0.0, 0.02, 0.001, 100, 1e-10, 1e30, printAllCurrents)
rampbias(device, "drain", 0.4, 0.01, 0.001, 100, 1e-10, 1e30, printAllCurrents)    ##current unit is A/cm


write_devices(file="gmsh_HEMT2d_dd.dat", type="tecplot")

