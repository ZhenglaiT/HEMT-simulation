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
from devsim.python_packages.model_create import *

#Create the drift-diffusion with Possion dependent band structure 
#Calculate the equivalent potential to use same discrete form
def CreateEquivalentPotential(device,region):
    Eff_potential_n="Potential+Affinity+V_t_T*log(nstate_density)"
    Eff_potential_p="Potential+Affinity+Eg-V_t_T*log(pstate_density)"

    CreateNodeModel(device, region, "E_p_n", Eff_potential_n)
    CreateNodeModel(device, region, "E_p_p", Eff_potential_p)
    CreateNodeModelDerivative(device, region, "E_p_n", Eff_potential_n, "Potential")
    CreateNodeModelDerivative(device, region, "E_p_n", Eff_potential_n, "Temperature")
    CreateNodeModelDerivative(device, region, "E_p_p", Eff_potential_p, "Potential")
    CreateNodeModelDerivative(device, region, "E_p_p", Eff_potential_p, "Temperature")

    edge_from_node_model(device=device, region=region, node_model="E_p_n")
    edge_from_node_model(device=device, region=region, node_model="E_p_p")

def CreateBernoulli_Ep (device, region):
    EnsureEdgeFromNodeModelExists(device, region, "Potential")
    EnsureEdgeFromNodeModelExists(device, region, "Temperature")
    EnsureEdgeFromNodeModelExists(device, region, "Potential")
    EnsureEdgeFromNodeModelExists(device, region, "E_p_n")
    EnsureEdgeFromNodeModelExists(device, region, "E_p_p")

    vdiffstr_n="(E_p_n@n0 - E_p_n@n1)/edgeV_t_T"
    CreateEdgeModel(device, region, "vdiffn", vdiffstr_n)
    CreateEdgeModel(device, region, "vdiffn:Potential@n0",  "edgeV_t_T^(-1)")
    CreateEdgeModel(device, region, "vdiffn:Potential@n1",  "-vdiffn:Potential@n0")
    CreateEdgeModel(device, region, "Bern01_n",              "B(vdiffn)")
    CreateEdgeModel(device, region, "Bern01_n:Potential@n0", "dBdx(vdiffn) * vdiffn:Potential@n0")
    CreateEdgeModel(device, region, "Bern01_n:Potential@n1", "-Bern01_n:Potential@n0")

    CreateEdgeModelDerivatives(device, region,"vdiffn",vdiffstr_n,"Temperature")
    CreateEdgeModelDerivatives(device, region,"Bern01_n","B(vdiffn)","Temperature")

    vdiffstr_p="(E_p_p@n0 - E_p_p@n1)/edgeV_t_T"
    CreateEdgeModel(device, region, "vdiffp", vdiffstr_p)
    CreateEdgeModel(device, region, "vdiffp:Potential@n0",  "edgeV_t_T^(-1)")
    CreateEdgeModel(device, region, "vdiffp:Potential@n1",  "-vdiffp:Potential@n0")
    CreateEdgeModel(device, region, "Bern01_p",              "B(vdiffp)")
    CreateEdgeModel(device, region, "Bern01_p:Potential@n0", "dBdx(vdiffp) * vdiffp:Potential@n0")
    CreateEdgeModel(device, region, "Bern01_p:Potential@n1", "-Bern01_p:Potential@n0")

    CreateEdgeModelDerivatives(device, region,"vdiffp",vdiffstr_p,"Temperature")
    CreateEdgeModelDerivatives(device, region,"Bern01_p","B(vdiffp)","Temperature")

def CreateElectronCurrent_BS(device, region,mu_n):
    '''
    Electron current
    '''
    EnsureEdgeFromNodeModelExists(device, region, "Potential")
    EnsureEdgeFromNodeModelExists(device, region, "Electrons")
    EnsureEdgeFromNodeModelExists(device, region, "Holes")
    # Make sure the bernoulli functions exist
    if not InEdgeModelList(device, region, "Bern01_n"):
        CreateBernoulli_Ep(device, region)
    #### test for requisite models here

    Jn = "ElectronCharge*{0}*EdgeInverseLength*edgeV_t_T*kahan3(Electrons@n1*Bern01_n,  Electrons@n1*vdiffn,  -Electrons@n0*Bern01_n)".format(mu_n)

    CreateEdgeModel(device, region, "ElectronCurrent", Jn)
    for i in ("Electrons", "Potential", "Holes", "Temperature"):
        CreateEdgeModelDerivatives(device, region, "ElectronCurrent", Jn, i)



def CreateHoleCurrent_BS(device, region,mu_p):
    '''
    Hole current
    '''
    EnsureEdgeFromNodeModelExists(device, region, "Potential")
    EnsureEdgeFromNodeModelExists(device, region, "Holes")
    # Make sure the bernoulli functions exist

    if not InEdgeModelList(device, region, "Bern01_p"):
        CreateBernoulli_Ep(device, region)
    ##### test for requisite models here


    Jp ="-ElectronCharge*{0}*EdgeInverseLength*edgeV_t_T*kahan3(Holes@n1*Bern01_p, -Holes@n0*Bern01_p, -Holes@n0*vdiffp)".format(mu_p)

    CreateEdgeModel(device, region, "HoleCurrent", Jp)
    for i in ("Holes", "Potential", "Electrons", "Temperature"):
        CreateEdgeModelDerivatives(device, region, "HoleCurrent", Jp, i)

