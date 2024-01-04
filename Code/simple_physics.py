#Copyright 2013 Devsim LLC
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

from simple_dd import *
from devsim import *
#TODO: make this a class so that paramters can be changed
contactcharge_node="contactcharge_node"
contactcharge_edge="contactcharge_edge"
ece_name="ElectronContinuityEquation"
hce_name="HoleContinuityEquation"
celec_model = "(1e-10 + 0.5*abs(NetDoping+(NetDoping^2 + 4 * n_i^2)^(0.5)))"
chole_model = "(1e-10 + 0.5*abs(-NetDoping+(NetDoping^2 + 4 * n_i^2)^(0.5)))"
# celec_model = "(0.5*abs(NetDoping+(NetDoping^2 + 4 * n_i^2)^(0.5)))"
# chole_model = "(0.5*abs(-NetDoping+(NetDoping^2 + 4 * n_i^2)^(0.5)))"

q      = 1.6e-19 # coul
k      = 1.3806503e-23 # J/K
eps_0  = 8.85418e-14 # F/cm^2
#T      = 300 # K
eps_si = 11.3
eps_ox = 3.9
# TODO: make this temperature dependent
n_i    = 1.0e10 # #/cm^3
# constant in our approximation
mu_n   = 400
mu_p   = 200
Sithermalconductivity = 1.48  ##W/K*cm
Oxidethermalconductivity = 0.014 
GaNthermalconductivity = 1.3
AlNthermalconductivity = 2.85
Pi = 3.1415926
h  = 6.62607015e-34 #J·s
m0 =9.10938356e-31  #kg
# heatsource = -1e8    ##w/cm^2
n_i_GaN = 1.0618e-10
n_i_AlGaN = 3.126e-15        #为了简便，对于这两种材料的本征载流子浓度都是基于T=300K计算的
#Boltzmann statistics :ni=sqrt(Nc*Nv)*exp(-Eg/2kT)

ns_AlGaN=2.7946e13        #sheet charge, 正号代表总电荷为正，负号代表总电荷为负
ns_GaN=-2.4690e13        #正号代表总电荷为正，负号代表总电荷为负

# ns_AlGaN=0
# ns_GaN=0


def GetContactBiasName(contact):
    return "{0}_bias".format(contact)

def GetContactTemperatureName(contact):
    return "{0}_biasT".format(contact)

def GetContactNodeModelName(contact):
    return "{0}nodemodel".format(contact)

def GetContactTemperatureNodeModelName(contact):
    return "{0}Tnodemodel".format(contact)

def PrintCurrents(device, contact):
    '''
       print out contact currents
    '''
    # TODO add charge
    contact_bias_name = GetContactBiasName(contact)
    electron_current  = get_contact_current(device=device, contact=contact, equation=ece_name)
    hole_current      = get_contact_current(device=device, contact=contact, equation=hce_name)
    total_current     = electron_current + hole_current                                        
    voltage           = get_parameter(device=device, name=GetContactBiasName(contact))
    print("{0}\t{1}\t{2}\t{3}\t{4}".format(contact, voltage, electron_current, hole_current, total_current))

#in the future, worry about workfunction
def CreateOxideContact(device, region, contact):
    conteq="Permittivity*ElectricField"
    contact_bias_name  = GetContactBiasName(contact)
    contact_model_name = GetContactNodeModelName(contact)
    eq = "Potential - {0}".format(contact_bias_name)
    CreateContactNodeModel(device, contact, contact_model_name, eq)
    CreateContactNodeModelDerivative(device, contact, contact_model_name, eq, "Potential")

    #TODO: make everyone use dfield
    if not InEdgeModelList(device, region, contactcharge_edge):
        CreateEdgeModel(device, region, contactcharge_edge, "Permittivity*ElectricField")
        CreateEdgeModelDerivatives(device, region, contactcharge_edge, "Permittivity*ElectricField", "Potential")

    contact_equation(device=device, contact=contact, name="PotentialEquation",
                     node_model=contact_model_name, edge_charge_model= contactcharge_edge)

def CreateDensityofStates(device,region):
    #记得在计算热电功率时把里面的态密度给取出来
    # Nc = "2*pow((2*pi*DOS_emasses*m0*boltzmannConstant*Temperature/planckconstant^2),1.5)"    内部还需要乘一个电子质量
    # Nv = "2*pow((2*pi*DOS_hmasses*m0*boltzmannConstant*Temperature/planckconstant^2),1.5)"
    #下面化简后的式子来自于“On the hole effective mass and the free hole statistics in wurtzite GaN”
    Nc  ="4.83e15*pow((DOS_emasses*Temperature),1.5)"     #cm-3
    Nv  ="4.83e15*pow((DOS_hmasses*Temperature),1.5)"

    CreateNodeModel           (device,  region, "nstate_density", Nc)
    CreateNodeModel           (device,  region, "pstate_density", Nv)
    CreateNodeModelDerivative (device,  region, "nstate_density", Nc, "Temperature")
    CreateNodeModelDerivative (device,  region, "pstate_density", Nv, "Temperature")
    CreateSolution(device, region, "nstate_density:Electrons")
    CreateSolution(device, region, "pstate_density:Holes")
   


def SetAlGaNParameters(device,region,xcomp):
    set_parameter(device=device, region=region, name="xcomp",   value=xcomp)
    set_parameter(device=device, region=region, name="Permittivity",   value=8.5*xcomp+8.9*(1-xcomp))
    set_parameter(device=device, region=region, name="Eg",   value=3.9727)
    set_parameter(device=device, region=region, name="Affinity",   value=3.7234)
    set_parameter(device=device, region=region, name="DOS_emasses",   value=0.314*xcomp+0.2*(1-xcomp))
    set_parameter(device=device, region=region, name="DOS_hmasses",   value=0.417*xcomp+1.0*(1-xcomp))
    set_parameter(device=device, region=region, name="ElectronCharge", value=q)
    set_parameter(device=device, region=region, name="boltzmannConstant",              value=k)  
    set_parameter(device=device, region=region, name="pi",              value=Pi)    
    set_parameter(device=device, region=region, name="planckconstant",              value=h)    
    set_parameter(device=device, region=region, name="n_i",            value=n_i_AlGaN)
    set_parameter(device=device, region=region, name="GaNthermalconductivity",     value=GaNthermalconductivity)
    set_parameter(device=device, region=region, name="AlNthermalconductivity",     value=AlNthermalconductivity)
    set_parameter(device=device, region=region, name="phix",     value=3.75*xcomp**0.4*(1-xcomp)**0.27/(1-0.31*xcomp)**0.5)
    set_parameter(device=device, region=region, name="holes_mobility",     value=81) #For Albrecht
    set_parameter(device=device, region=region, name="n1", value=n_i_AlGaN)
    set_parameter(device=device, region=region, name="p1", value=n_i_AlGaN)
    set_parameter(device=device, region=region, name="taun", value=1.2e-8)  #According to the Silvaco, tn=tp=t_GaN
    set_parameter(device=device, region=region, name="taup", value=1.2e-8) 
    # set_parameter(device=device, region=region, name="surface_charge", value=ns) 
    # set_parameter(device=device, region=region, name="An", value=24.0351) #Richardson constants
    # set_parameter(device=device, region=region, name="Ap", value=120.1753) #Richardson constants
    set_parameter(device=device, region=region, name="An", value=28.1451) #Richardson constants
    set_parameter(device=device, region=region, name="Ap", value=99.1566) #Richardson constants
    # The choice of An in article or guide books is contradicted. In silvaco, use the An of both sides.
    # In Charon, use the same parameter defined by user for both sides. 
    # Some article use the An of the material on the side with the smallest equivalent mass (Solid-State Electronics, 1993, 36(3): 321-330.)

def SetGaNParameters(device,region):
    set_parameter(device=device, region=region, name="Permittivity",   value=8.9)
    set_parameter(device=device, region=region, name="Eg",   value=3.4346)
    set_parameter(device=device, region=region, name="Affinity",   value=4.1) #eV
    set_parameter(device=device, region=region, name="DOS_emasses",   value=0.2)
    set_parameter(device=device, region=region, name="DOS_hmasses",   value=1.0)
    set_parameter(device=device, region=region, name="ElectronCharge", value=q)
    set_parameter(device=device, region=region, name="boltzmannConstant",              value=k)  
    set_parameter(device=device, region=region, name="pi",              value=Pi)    
    set_parameter(device=device, region=region, name="planckconstant",              value=h)    
    set_parameter(device=device, region=region, name="n_i",            value=n_i_GaN)
    set_parameter(device=device, region=region, name="GaNthermalconductivity",     value=GaNthermalconductivity)
    set_parameter(device=device, region=region, name="holes_mobility",     value=6.44) #For Albrecht
    set_parameter(device=device, region=region, name="n1", value=n_i_GaN)
    set_parameter(device=device, region=region, name="p1", value=n_i_GaN)
    set_parameter(device=device, region=region, name="taun", value=1.2e-8) 
    set_parameter(device=device, region=region, name="taup", value=1.2e-8) 
    # set_parameter(device=device, region=region, name="surface_charge", value=ns) 
    set_parameter(device=device, region=region, name="An", value=24.0351) #Richardson constants
    set_parameter(device=device, region=region, name="Ap", value=120.1753) #Richardson constants
#####
##### Constants
#####
def SetOxideParameters(device, region):
    '''
      Sets physical parameters
    '''
    set_parameter(device=device, region=region, name="Permittivity",   value=eps_ox * eps_0)
    set_parameter(device=device, region=region, name="ElectronCharge", value=q)
    set_parameter(device=device, region=region, name="Oxidethermalconductivity",     value=Oxidethermalconductivity)

def SetSiliconParameters(device, region):
    '''
      Sets physical parameters assuming constants
    '''
    #### TODO: make T a free parameter and T dependent parameters as models
    set_parameter(device=device, region=region, name="Permittivity",   value=eps_si * eps_0)
    set_parameter(device=device, region=region, name="ElectronCharge", value=q)
    set_parameter(device=device, region=region, name="n_i",            value=n_i)
    #set_parameter(device=device, region=region, name="T",              value=T)
    set_parameter(device=device, region=region, name="boltzmannConstant",              value=k)  
    set_parameter(device=device, region=region, name="mu_n",           value=mu_n)
    set_parameter(device=device, region=region, name="mu_p",           value=mu_p)
    set_parameter(device=device, region=region, name="Sithermalconductivity",     value=Sithermalconductivity)
    #set_parameter(device=device, region=region, name="source",     value=heatsource)

    #default SRH parameters
    set_parameter(device=device, region=region, name="n1", value=n_i)
    set_parameter(device=device, region=region, name="p1", value=n_i)
    set_parameter(device=device, region=region, name="taun", value=1e-7)  #The parameter in devsim is 1e-5
    set_parameter(device=device, region=region, name="taup", value=1e-7)  #But the parameter in silvaco is 1e-7


def CreateSiliconPotentialOnly(device, region):
    '''
      Creates the physical models for a Silicon region
    '''
    if not InNodeModelList(device, region, "Potential"):
        print("Creating Node Solution Potential")
        CreateSolution(device, region, "Potential")
    elec_i = "n_i*exp(Potential/V_t_T)"
    hole_i = "n_i^2/IntrinsicElectrons"
    charge_i = "kahan3(IntrinsicHoles, -IntrinsicElectrons, NetDoping)"
    pcharge_i = "-ElectronCharge * IntrinsicCharge"

    # require NetDoping
    for i in (
        ("IntrinsicElectrons", elec_i),
        ("IntrinsicHoles", hole_i),
        ("IntrinsicCharge", charge_i),
        ("PotentialIntrinsicCharge", pcharge_i)
    ):
        n = i[0]
        e = i[1]
        CreateNodeModel(device, region, n, e)
        CreateNodeModelDerivative(device, region, n, e, "Potential")
        CreateNodeModelDerivative(device, region, n, e, "Temperature")


    ### TODO: Edge Average Model
    for i in (
        ("ElectricField",     "(Potential@n0-Potential@n1)*EdgeInverseLength"),
        ("PotentialEdgeFlux", "Permittivity * ElectricField")
    ):
        n = i[0]
        e = i[1]
        CreateEdgeModel(device, region, n, e)
        CreateEdgeModelDerivatives(device, region, n, e, "Potential")

    equation(device=device, region=region, name="PotentialEquation", variable_name="Potential",
             node_model="PotentialIntrinsicCharge", edge_model="PotentialEdgeFlux", variable_update="log_damp")


def CreatAlGaNPotentialOnly(device, region):
#   The sheet charge is determined by the polarization and 2DEG.
#   For AlGaN, positive charge is accumlated in the bottom due to polarization.
    set_parameter(device=device,region=region,name="surface_charge", value=ns_AlGaN) 

    if not InNodeModelList(device, region, "Potential"):
        print("Creating Node Solution Potential")
        CreateSolution(device, region, "Potential")
    elec_i = "n_i*exp(Potential/V_t_T)"
    hole_i = "n_i^2/IntrinsicElectrons"
    # charge_i = "kahan3(IntrinsicHoles, -IntrinsicElectrons, NetDoping)"
    Sheet_charge="surface_charge * SurfaceArea / NodeVolume"
    charge_i = "kahan4(IntrinsicHoles, -IntrinsicElectrons, NetDoping,Sheet_charge)"
    pcharge_i = "-ElectronCharge * IntrinsicCharge"

    for i in (
        ("IntrinsicElectrons", elec_i),
        ("IntrinsicHoles", hole_i),
        ("Sheet_charge", Sheet_charge),
        ("IntrinsicCharge", charge_i),
        ("PotentialIntrinsicCharge", pcharge_i)
    ):
        n = i[0]
        e = i[1]
        CreateNodeModel(device, region, n, e)
        CreateNodeModelDerivative(device, region, n, e, "Potential")
        CreateNodeModelDerivative(device, region, n, e, "Temperature")

    for i in (
        ("ElectricField",     "(Potential@n0-Potential@n1)*EdgeInverseLength"),
      ("PotentialEdgeFlux", "Permittivity * ElectricField")
    ):
        n = i[0]
        e = i[1]
        CreateEdgeModel(device, region, n, e)
        CreateEdgeModelDerivatives(device, region, n, e, "Potential")

    equation(device=device, region=region, name="PotentialEquation", variable_name="Potential",
             node_model="PotentialIntrinsicCharge", edge_model="PotentialEdgeFlux", variable_update="log_damp")


def CreatGaNPotentialOnly(device, region):
#   The sheet charge is determined by the polarization and 2DEG.
#   For GaN, negative charge is accumlated in the top due to polarization and 2DEG.
    set_parameter(device=device,region=region, name="surface_charge", value=ns_GaN) 

    if not InNodeModelList(device, region, "Potential"):
        print("Creating Node Solution Potential")
        CreateSolution(device, region, "Potential")
    elec_i = "n_i*exp(Potential/V_t_T)"
    hole_i = "n_i^2/IntrinsicElectrons"
    # charge_i = "kahan3(IntrinsicHoles, -IntrinsicElectrons, NetDoping)"
    Sheet_charge="surface_charge * SurfaceArea / NodeVolume"
    charge_i = "kahan4(IntrinsicHoles, -IntrinsicElectrons, NetDoping,Sheet_charge)"
    pcharge_i = "-ElectronCharge * IntrinsicCharge"

    for i in (
        ("IntrinsicElectrons", elec_i),
        ("IntrinsicHoles", hole_i),
        ("Sheet_charge", Sheet_charge),
        ("IntrinsicCharge", charge_i),
        ("PotentialIntrinsicCharge", pcharge_i)
    ):
        n = i[0]
        e = i[1]
        CreateNodeModel(device, region, n, e)
        CreateNodeModelDerivative(device, region, n, e, "Potential")
        CreateNodeModelDerivative(device, region, n, e, "Temperature")

    for i in (
        ("ElectricField",     "(Potential@n0-Potential@n1)*EdgeInverseLength"),
      ("PotentialEdgeFlux", "Permittivity * ElectricField")
    ):
        n = i[0]
        e = i[1]
        CreateEdgeModel(device, region, n, e)
        CreateEdgeModelDerivatives(device, region, n, e, "Potential")

    equation(device=device, region=region, name="PotentialEquation", variable_name="Potential",
             node_model="PotentialIntrinsicCharge", edge_model="PotentialEdgeFlux", variable_update="log_damp")




def CreateSiliconPotentialOnlyContact(device, region, contact, is_circuit=False):
    '''
      Creates the potential equation at the contact
      if is_circuit is true, than use node given by GetContactBiasName
    '''
    # Means of determining contact charge
    # Same for all contacts
    if not InNodeModelList(device, region, "contactcharge_node"):
        CreateNodeModel(device, region, "contactcharge_node", "ElectronCharge*IntrinsicCharge")
    #### TODO: This is the same as D-Field
    if not InEdgeModelList(device, region, "contactcharge_edge"):
        CreateEdgeModel(device, region, "contactcharge_edge", "Permittivity*ElectricField")
        CreateEdgeModelDerivatives(device, region, "contactcharge_edge", "Permittivity*ElectricField", "Potential")

#  set_parameter(device=device, region=region, name=GetContactBiasName(contact), value=0.0)

    contact_model = "Potential -{0} + ifelse(NetDoping > 0, \
    -V_t_T*log({1}/n_i), \
    V_t_T*log({2}/n_i))".format(GetContactBiasName(contact), celec_model, chole_model)

    contact_model_name = GetContactNodeModelName(contact)
    CreateContactNodeModel(device, contact, contact_model_name, contact_model)
    # Simplify it too complicated
    CreateContactNodeModel(device, contact, "{0}:{1}".format(contact_model_name,"Potential"), "1")
    CreateContactNodeModel(device, contact, "{0}:{1}".format(contact_model_name,"Temperature"), "ifelse(NetDoping > 0, \
    -boltzmannConstant/(ElectronCharge)*log({0}/n_i), \
    boltzmannConstant/(ElectronCharge)*log({1}/n_i))".format(celec_model, chole_model))

    if is_circuit:
        CreateContactNodeModel(device, contact, "{0}:{1}".format(contact_model_name,GetContactBiasName(contact)), "-1")

    if is_circuit:
        contact_equation(device=device, contact=contact, name="PotentialEquation",
                         node_model=contact_model_name, edge_model="",
                         node_charge_model="", edge_charge_model="contactcharge_edge",
                         node_current_model="", edge_current_model="", circuit_node=GetContactBiasName(contact))
    else:
        contact_equation(device=device, contact=contact, name="PotentialEquation",
                         node_model=contact_model_name, edge_model="",
                         node_charge_model="", edge_charge_model="contactcharge_edge",
                         node_current_model="", edge_current_model="")

def CreateSchottkyContacts(device,region,contact,workfun):

    if not InEdgeModelList(device, region, "contactcharge_edge"):
        CreateEdgeModel(device, region, "contactcharge_edge", "Permittivity*ElectricField")
        CreateEdgeModelDerivatives(device, region, "contactcharge_edge", "Permittivity*ElectricField", "Potential")
    
    contact_model = "Potential- Affinity - Eg/2 - \
                    boltzmannConstant*Temperature/(2*ElectronCharge)*log(nstate_density/pstate_density) + \
                    ({workfun}) - ({contactbias})".format(workfun=workfun,contactbias=GetContactBiasName(contact))
    #The unit of Affinity, Eg and workfun are eV, assume q=1e to get the V.

    contact_model_name = GetContactNodeModelName(contact)
    CreateContactNodeModel(device, contact, contact_model_name, contact_model)
    CreateContactNodeModel(device, contact, "{0}:{1}".format(contact_model_name,"Potential"), "1")
    CreateContactNodeModelDerivative(device,contact,contact_model_name,contact_model,"Temperature")

    contact_equation(device=device, contact=contact, name="PotentialEquation",
                        node_model=contact_model_name, edge_model="",
                        node_charge_model="", edge_charge_model="contactcharge_edge",
                        node_current_model="", edge_current_model="")


def CreateSRH(device, region):
    RAuger="2.8e-31*(Holes*Electrons^2-Electrons*n1^2)+9.9e-32*(Electrons*Holes^2-Holes*n1^2)" #俄歇复合
    CreateNodeModel(device, region, "RAuger", RAuger)    
    USRH="(Electrons*Holes - n_i^2)/(taup*(Electrons + n1) + taun*(Holes + p1))"  #net recombination rate R-G
    Gn = "-ElectronCharge * USRH"
    Gp = "+ElectronCharge * USRH"
    CreateNodeModel(device, region, "USRH", USRH)
    CreateNodeModel(device, region, "ElectronGeneration", Gn)
    CreateNodeModel(device, region, "HoleGeneration", Gp)
    for i in ("Electrons", "Holes"):
        CreateNodeModelDerivative(device, region, "USRH", USRH, i)
        CreateNodeModelDerivative(device, region, "ElectronGeneration", Gn, i)
        CreateNodeModelDerivative(device, region, "HoleGeneration", Gp, i)

def CreateECE(device, region,mu_n):
    CreateElectronCurrent(device, region,mu_n)
    NCharge = "-ElectronCharge * Electrons"
    CreateNodeModel(device, region, "NCharge", NCharge)
    CreateNodeModelDerivative(device, region, "NCharge", NCharge, "Electrons")

    equation(device=device, region=region, name="ElectronContinuityEquation", variable_name="Electrons",
            time_node_model = "NCharge",
            edge_model="ElectronCurrent", variable_update="positive", node_model="ElectronGeneration")

def CreateHCE(device, region,mu_p):
    CreateHoleCurrent(device, region, mu_p)
    PCharge = "ElectronCharge * Holes"
    CreateNodeModel(device, region, "PCharge", PCharge)
    CreateNodeModelDerivative(device, region, "PCharge", PCharge, "Holes")

    equation(device=device, region=region, name="HoleContinuityEquation", variable_name="Holes",
             time_node_model = "PCharge",
             edge_model="HoleCurrent", variable_update="positive", node_model="HoleGeneration")

def CreatePE(device, region):
    pne = "-ElectronCharge*kahan3(Holes, -Electrons, NetDoping)"
    CreateNodeModel(device, region, "PotentialNodeCharge", pne)
    CreateNodeModelDerivative(device, region, "PotentialNodeCharge", pne, "Electrons")
    CreateNodeModelDerivative(device, region, "PotentialNodeCharge", pne, "Holes")

    equation(device=device, region=region, name="PotentialEquation", variable_name="Potential",
             node_model="PotentialNodeCharge", edge_model="PotentialEdgeFlux",
             time_node_model="", variable_update="log_damp")

def CreateSiliconDriftDiffusion(device, region,mu_n="mu_n", mu_p="mu_p"):
    CreatePE(device, region)
    CreateBernoulli(device, region)
    CreateSRH(device, region)
    CreateECE(device, region, mu_n)
    CreateHCE(device, region, mu_p)


def CreateAlGaNPE(device, region):
    pne = "-ElectronCharge*kahan4(Holes, -Electrons, NetDoping,Sheet_charge)"
    CreateNodeModel(device, region, "PotentialNodeCharge", pne)
    CreateNodeModelDerivative(device, region, "PotentialNodeCharge", pne, "Electrons")
    CreateNodeModelDerivative(device, region, "PotentialNodeCharge", pne, "Holes")

    equation(device=device, region=region, name="PotentialEquation", variable_name="Potential",
             node_model="PotentialNodeCharge", edge_model="PotentialEdgeFlux",
             time_node_model="", variable_update="log_damp")
    
def CreateGaNPE(device, region):
    pne = "-ElectronCharge*kahan4(Holes, -Electrons, NetDoping,Sheet_charge)"
    CreateNodeModel(device, region, "PotentialNodeCharge", pne)
    CreateNodeModelDerivative(device, region, "PotentialNodeCharge", pne, "Electrons")
    CreateNodeModelDerivative(device, region, "PotentialNodeCharge", pne, "Holes")

    equation(device=device, region=region, name="PotentialEquation", variable_name="Potential",
             node_model="PotentialNodeCharge", edge_model="PotentialEdgeFlux",
             time_node_model="", variable_update="log_damp")

def CreateAlGaNDriftDiffusion(device, region,mu_n="mu_n", mu_p="mu_p"):
    CreateEquivalentPotential(device, region)
    CreateAlGaNPE(device, region)
    CreateBernoulli_Ep(device, region)
    CreateBernoulli(device, region)
    CreateSRH(device, region)
    CreateECE(device, region,mu_n)
    CreateHCE(device, region,mu_p)

def CreateGaNDriftDiffusion(device, region,mu_n="mu_n", mu_p="mu_p"):
    CreateEquivalentPotential(device, region)
    CreateGaNPE(device, region)
    CreateBernoulli_Ep(device, region)
    CreateBernoulli(device, region)
    CreateSRH(device, region)
    CreateECE(device, region,mu_n)
    CreateHCE(device, region,mu_p)


def CreateSiliconDriftDiffusionAtContact(device, region, contact, is_circuit=False): 
    '''
      Restrict electrons and holes to their equilibrium values
      Integrates current into circuit
    '''
    contact_electrons_model = "Electrons - ifelse(NetDoping > 0, {0}, n_i^2/{1})".format(celec_model, chole_model)
    contact_holes_model = "Holes - ifelse(NetDoping < 0, +{1}, +n_i^2/{0})".format(celec_model, chole_model)
    contact_electrons_name = "{0}nodeelectrons".format(contact)
    contact_holes_name = "{0}nodeholes".format(contact)

    CreateContactNodeModel(device, contact, contact_electrons_name, contact_electrons_model)
    #TODO: The simplification of the ifelse statement is time consuming
#   CreateContactNodeModelDerivative(device, contact, contact_electrons_name, contact_electrons_model, "Electrons")
    CreateContactNodeModel(device, contact, "{0}:{1}".format(contact_electrons_name, "Electrons"), "1")

    CreateContactNodeModel(device, contact, contact_holes_name, contact_holes_model)
    CreateContactNodeModel(device, contact, "{0}:{1}".format(contact_holes_name, "Holes"), "1")

    #TODO: keyword args
    if is_circuit:
        contact_equation(device=device, contact=contact, name="ElectronContinuityEquation",
                         node_model=contact_electrons_name,
                         edge_current_model="ElectronCurrent", circuit_node=GetContactBiasName(contact))

        contact_equation(device=device, contact=contact, name="HoleContinuityEquation",
                         node_model=contact_holes_name,
                         edge_current_model="HoleCurrent", circuit_node=GetContactBiasName(contact))

    else:
        contact_equation(device=device, contact=contact, name="ElectronContinuityEquation",
                         node_model=contact_electrons_name,
                         edge_current_model="ElectronCurrent")

        contact_equation(device=device, contact=contact, name="HoleContinuityEquation",
                         node_model=contact_holes_name,
                         edge_current_model="HoleCurrent")


def CreateDriftDiffusionAtSchottkyContact(device,region,contact,workfun):
    V_surfe = "110*Temperature^2/ElectronCharge/nstate_density"
    V_surfh = "30*Temperature^2/ElectronCharge/pstate_density"

    CreateContactNodeModel           (device,  contact, "V_surfe", V_surfe)
    CreateContactNodeModel           (device,  contact, "V_surfh", V_surfh)
    CreateContactNodeModelDerivative (device,  contact, "V_surfe", V_surfe, "Temperature")
    CreateContactNodeModelDerivative (device,  contact, "V_surfh", V_surfh, "Temperature")    
    
    deltaPhi="Affinity + Eg/2 + \
                    boltzmannConstant*Temperature/(2*ElectronCharge)*log(nstate_density/pstate_density) - \
                    {0}".format(workfun)
    
    CreateContactNodeModel           (device,  contact, "deltaPhi", deltaPhi)
    CreateContactNodeModelDerivative (device,  contact, "deltaPhi", deltaPhi, "Temperature")

    Electrons_eq = "n_i*exp(ElectronCharge*(deltaPhi)/boltzmannConstant/Temperature)"
    Holes_eq =  "n_i*exp(-ElectronCharge*(deltaPhi)/boltzmannConstant/Temperature)"

    CreateContactNodeModel           (device,  contact, "Electrons_eq", Electrons_eq)
    CreateContactNodeModel           (device,  contact, "Holes_eq", Holes_eq)
    CreateContactNodeModelDerivative (device,  contact, "Electrons_eq", Electrons_eq, "Temperature")
    CreateContactNodeModelDerivative (device,  contact, "Holes_eq", Holes_eq, "Temperature")    

    contact_electrons_model ="ElectronCharge*V_surfe*(Electrons-Electrons_eq) * ContactSurfaceArea/NodeVolume"
    contact_holes_model     ="ElectronCharge*V_surfh*(Holes-Holes_eq) * ContactSurfaceArea/NodeVolume"

    contact_electrons_name="{0}_electrons_schottky".format(contact)
    contact_holes_name="{0}_holes_schottky".format(contact)

    CreateContactNodeModel(device, contact, contact_electrons_name, contact_electrons_model)
    CreateContactNodeModelDerivative(device, contact, contact_electrons_name, contact_electrons_model, "Electrons")
    CreateContactNodeModelDerivative(device, contact, contact_electrons_name, contact_electrons_model, "Temperature")
    CreateContactNodeModel(device, contact, contact_holes_name, contact_holes_model)
    CreateContactNodeModelDerivative(device, contact, contact_holes_name, contact_holes_model, "Holes")
    CreateContactNodeModelDerivative(device, contact, contact_holes_name, contact_holes_model, "Temperature")


    contact_equation(device=device, contact=contact, name="ElectronContinuityEquation", 
                    node_model=contact_electrons_name, edge_model="ElectronCurrent", edge_current_model="ElectronCurrent")
    contact_equation(device=device, contact=contact, name="HoleContinuityEquation",
                    node_model=contact_holes_name, edge_model="HoleCurrent", edge_current_model="HoleCurrent")


def CreateOxidePotentialOnly(device, region, update_type="default"):
    '''
      Create electric field model in oxide
      Creates Potential solution variable if not available
    '''
    if not InNodeModelList(device, region, "Potential"):
        print("Creating Node Solution Potential")
        CreateSolution(device, region, "Potential")

    efield="(Potential@n0 - Potential@n1)*EdgeInverseLength"
    # this needs to remove derivatives w.r.t. independents
    CreateEdgeModel(device, region, "ElectricField", efield)
    CreateEdgeModelDerivatives(device, region, "ElectricField", efield, "Potential")
    dfield="Permittivity*ElectricField"
    CreateEdgeModel(device, region, "PotentialEdgeFlux", dfield)
    CreateEdgeModelDerivatives(device, region, "PotentialEdgeFlux", dfield, "Potential")
    equation(device=device, region=region, name="PotentialEquation", variable_name="Potential",
             edge_model="PotentialEdgeFlux", variable_update=update_type)

def CreateSiliconOxideInterface(device, interface):
    '''
      continuous potential at interface
    '''
    model_name = CreateContinuousInterfaceModel(device, interface, "Potential")
    interface_equation(device=device, interface=interface, name="PotentialEquation", interface_model=model_name, type="continuous")

#
##TODO: similar model for silicon/silicon interface
## should use quasi-fermi potential
def CreateSiliconSiliconInterface(device, interface):
    '''
      Enforces potential, electron, and hole continuity across the interface
    '''
    CreateSiliconOxideInterface(device, interface)
    ename = CreateContinuousInterfaceModel(device, interface, "Electrons")
    interface_equation(device=device, interface=interface, name="ElectronContinuityEquation", interface_model=ename, type="continuous")
    hname = CreateContinuousInterfaceModel(device, interface, "Holes")
    interface_equation(device=device, interface=interface, name="HoleContinuityEquation", interface_model=hname, type="continuous")

def CreateHeterojunctionInterface(device, interface):
    #Make sure the continous of Potential
    model_name = CreateContinuousInterfaceModel(device, interface, "Potential")
    interface_equation(device=device, interface=interface, name="PotentialEquation", interface_model=model_name, type="continuous")


# # Continuous Quasi-Fermi Level Model
# def CreateHeterojunctionInterfaceContinous(device, interface):    
#     model_name = CreateContinuousInterfaceModel(device, interface, "Potential")
#     interface_equation(device=device, interface=interface, name="PotentialEquation", interface_model=model_name, type="continuous")



def CreateThermionicEmission(device,interface): 

    model_name = CreateContinuousInterfaceModel(device, interface, "Potential")
    interface_equation(device=device, interface=interface, name="PotentialEquation", interface_model=model_name, type="continuous")

    #create the thermionic emission in the Heterojunction interface.
    set_parameter(device=device, name="Delta_Ec",   value=0.3766)
    set_parameter(device=device, name="Delta_Ev",   value=0.1614)
    #the alignment of energy bands is determined by the affinity rule 

    Therma_Em_expn="ifelse(  Eg@r0>Eg@r1,\
        -(-An@r0*(Temperature@r0)^2/(nstate_density@r0)*Electrons@r0+\
        An@r1*(Temperature@r1)^2/(nstate_density@r1)*Electrons@r1*exp((-Delta_Ec*ElectronCharge)/(boltzmannConstant*Temperature@r1))),\
        -(-An@r0*(Temperature@r0)^2/(nstate_density@r0)*Electrons@r0*exp((-Delta_Ec*ElectronCharge)/(boltzmannConstant*Temperature@r0))+\
        An@r1*(Temperature@r1)^2/(nstate_density@r1)*Electrons@r1)  )"
    

    Therma_Em_expp="ifelse(  Eg@r0>Eg@r1,\
        -(Ap@r0*(Temperature@r0)^2/(pstate_density@r0)*Holes@r0-\
        Ap@r1*(Temperature@r1)^2/(pstate_density@r1)*Holes@r1*exp((-Delta_Ev*ElectronCharge)/(boltzmannConstant*Temperature@r1))) ,\
        -(Ap@r0*(Temperature@r0)^2/(pstate_density@r0)*Holes@r0*exp((-Delta_Ev*ElectronCharge)/(boltzmannConstant*Temperature@r0))-\
        Ap@r1*(Temperature@r1)^2/(pstate_density@r1)*Holes@r1) )"

    for name, equation in (
        ("Therma_Em_expn", Therma_Em_expn),
        ("Therma_Em_expp", Therma_Em_expp),
        ("Therma_Em_expn:Electrons@r0", "diff(%s,Electrons@r0)" % Therma_Em_expn),
        ("Therma_Em_expn:Electrons@r1", "diff(%s,Electrons@r1)" % Therma_Em_expn),
        ("Therma_Em_expp:Holes@r0", "diff(%s,Holes@r0)" % Therma_Em_expp),
        ("Therma_Em_expp:Holes@r1", "diff(%s,Holes@r1)" % Therma_Em_expp),
        ("Therma_Em_expn:Temperature@r0", "diff(%s,Temperature@r0)" % Therma_Em_expn),
        ("Therma_Em_expn:Temperature@r1", "diff(%s,Temperature@r1)" % Therma_Em_expn),
        ("Therma_Em_expp:Temperature@r0", "diff(%s,Temperature@r0)" % Therma_Em_expp),
        ("Therma_Em_expp:Temperature@r1", "diff(%s,Temperature@r1)" % Therma_Em_expp),
        ("Therma_Em_expn2", Therma_Em_expn),
        ("Therma_Em_expp2", Therma_Em_expp),
        ("Therma_Em_expn2:Electrons@r0", "diff(%s,Electrons@r0)" % Therma_Em_expn),
        ("Therma_Em_expn2:Electrons@r1", "diff(%s,Electrons@r1)" % Therma_Em_expn),
        ("Therma_Em_expp2:Holes@r0", "diff(%s,Holes@r0)" % Therma_Em_expp),
        ("Therma_Em_expp2:Holes@r1", "diff(%s,Holes@r1)" % Therma_Em_expp),
        ("Therma_Em_expn2:Temperature@r0", "diff(%s,Temperature@r0)" % Therma_Em_expn),
        ("Therma_Em_expn2:Temperature@r1", "diff(%s,Temperature@r1)" % Therma_Em_expn),
        ("Therma_Em_expp2:Temperature@r0", "diff(%s,Temperature@r0)" % Therma_Em_expp),
        ("Therma_Em_expp2:Temperature@r1", "diff(%s,Temperature@r1)" % Therma_Em_expp),
        ):
            interface_model(device=device, interface=interface, name=name, equation=equation)
            
    interface_equation(device=device, interface=interface, name="ElectronContinuityEquation", interface_model="Therma_Em_expn2", type="fluxterm")
    interface_equation(device=device, interface=interface, name="HoleContinuityEquation", interface_model="Therma_Em_expp2", type="fluxterm")



def CreateTinitial(device,region):
    initial_Temperature = "initialTemperature"
    for i in (
        ("Init_temperature", initial_Temperature),
    
    ):
        n = i[0]
        e = i[1]
        CreateNodeModel(device, region, n, e)
        # CreateNodeModelDerivative(device, region, n, e, "Temperature")

def CreateMobilityModels(device,region,mu_n,mu_p):
    mu_n_T = "mu_n*pow((edgeTemp)/(T),-2.2)"
    mu_p_T = "mu_p*pow((edgeTemp)/(T),-2.2)"
    for i in (
        ("mu_nT", mu_n_T),
        ("mu_pT", mu_p_T),
    ):
        n = i[0]
        e = i[1]
        CreateEdgeModel(device, region, n, e)
        CreateEdgeModelDerivatives(device, region,n,e,"Temperature")

def CreateSiThermalConductivity(device,region):
    Kappa_Si = "Sithermalconductivity*pow((edgeTemp)/(300),-1.65)"
    for i in (
        ("thermalconductivity", Kappa_Si),
    ):
        n = i[0]
        e = i[1]
        CreateEdgeModel(device, region, n, e)
        CreateEdgeModelDerivatives(device, region,n,e,"Temperature")

def CreateOxideThermalConductivity(device,region):
    Kappa_Oxide = "Oxidethermalconductivity*pow((edgeTemp)/(300),-1.65)"
    for i in (
        ("thermalconductivity", Kappa_Oxide),
    ):
        n = i[0]
        e = i[1]
        CreateEdgeModel(device, region, n, e)
        CreateEdgeModelDerivatives(device, region,n,e,"Temperature")

def CreateGaNThermalConductivity(device,region):
    Kappa_GaN = "GaNthermalconductivity*pow((edgeTemp)/(300),-0.28)"
    for i in (
        ("thermalconductivity", Kappa_GaN),
    ):
        n = i[0]
        e = i[1]
        CreateEdgeModel(device, region, n, e)
        CreateEdgeModelDerivatives(device, region,n,e,"Temperature")

def CreateAlGaNThermalConductivity(device,region):
    Kappa_GaN = "GaNthermalconductivity*pow((edgeTemp)/(300),-0.28)"
    Kappa_AlN = "AlNthermalconductivity*pow((edgeTemp)/(300),-1.64)"
    Kappa_AlGaN="(Kappa_GaN*(1-xcomp)+Kappa_AlN*xcomp)*exp(-phix)"

    for i in (
        ("GaNthermalconductivity", Kappa_GaN),
        ("AlNthermalconductivity", Kappa_AlN),
        ("thermalconductivity", Kappa_AlGaN),
    ):
        n = i[0]
        e = i[1]
        CreateEdgeModel(device, region, n, e)
        CreateEdgeModelDerivatives(device, region,n,e,"Temperature")


def Set_Mobility_Parameters(device, region):
    #As
    set_parameter(device=device, region=region, name="mu_min_e",  value=52.2)       #
    set_parameter(device=device, region=region, name="mu_max_e",  value=1417)       #
    set_parameter(device=device, region=region, name="theta1_e",     value=2.285)   #
    set_parameter(device=device, region=region, name="m_e",      value=1.0)  #
    set_parameter(device=device, region=region, name="f_BH",     value=3.828) #
    set_parameter(device=device, region=region, name="f_CW",     value=2.459)  #
    set_parameter(device=device, region=region, name="c_D",      value=0.21)   #
    set_parameter(device=device, region=region, name="Nref_D",   value=4.0e20)  #
    set_parameter(device=device, region=region, name="Nref_1_e", value=9.68e16)    #
    set_parameter(device=device, region=region, name="alpha_1_e",    value=0.68)  #

    #B
    set_parameter(device=device, region=region, name="mu_min_h",  value=44.9)   #
    set_parameter(device=device, region=region, name="mu_max_h",  value=470.5) #
    set_parameter(device=device, region=region, name="theta1_h",     value=2.247) #
    set_parameter(device=device, region=region, name="m_h",      value=1.258) #
    set_parameter(device=device, region=region, name="c_A",      value=0.50) #
    set_parameter(device=device, region=region, name="Nref_A",   value=7.2e20)  #
    set_parameter(device=device, region=region, name="Nref_1_h", value=2.23e17)    #2.23e17
    set_parameter(device=device, region=region, name="alpha_1_h",    value=0.719)   #


    #F_Pe F_Ph equations
    set_parameter(device=device, region=region, name="r1", value=0.7643)
    set_parameter(device=device, region=region, name="r2", value=2.2999)
    set_parameter(device=device, region=region, name="r3", value=6.5502)
    set_parameter(device=device, region=region, name="r4", value=2.3670)
    set_parameter(device=device, region=region, name="r5", value=-0.8522)  #-0.8522
    set_parameter(device=device, region=region, name="r6", value=0.6478)

    #G_Pe G_Ph equations
    set_parameter(device=device, region=region, name="s1", value=0.892333)
    set_parameter(device=device, region=region, name="s2", value=0.41372)
    set_parameter(device=device, region=region, name="s3", value=0.19778)
    set_parameter(device=device, region=region, name="s4", value=0.28227)
    set_parameter(device=device, region=region, name="s5", value=0.005978)
    set_parameter(device=device, region=region, name="s6", value=1.80618)
    set_parameter(device=device, region=region, name="s7", value=0.72169)

    # velocity saturation
    set_parameter(device=device, region=region, name="vsat_e", value=1.130e7)  #1.13e7
    set_parameter(device=device, region=region, name="vsat_h", value=1.130e7)  #1.13e7

    # #Lucent Mobility
    # set_parameter(device=device, region=region, name="alpha_e", value=6.85e-21)
    # set_parameter(device=device, region=region, name="alpha_h", value=7.82e-21)
    # set_parameter(device=device, region=region, name="A_e", value=2.58)
    # set_parameter(device=device, region=region, name="A_h", value=2.18)
    # set_parameter(device=device, region=region, name="B_e", value=3.61e7)
    # set_parameter(device=device, region=region, name="B_h", value=1.51e7)
    # set_parameter(device=device, region=region, name="C_e", value=1.70e4)
    # set_parameter(device=device, region=region, name="C_h", value=4.18e3)
    # set_parameter(device=device, region=region, name="kappa_e", value=1.7)
    # set_parameter(device=device, region=region, name="kappa_h", value=0.9)
    # set_parameter(device=device, region=region, name="tau_e", value=0.0233)
    # set_parameter(device=device, region=region, name="tau_h", value=0.0119)
    # set_parameter(device=device, region=region, name="eta_e", value=0.0767)
    # set_parameter(device=device, region=region, name="eta_h", value=0.123)
    # set_parameter(device=device, region=region, name="delta_e", value=3.58e18)
    # set_parameter(device=device, region=region, name="delta_h", value=4.10e15)
##create Klaassen model
def Klaassen_Mobility(device, region):
    # require Electrons, Holes, Donors, Acceptors already exist
    mu_L_e="(mu_max_e * (300 / Temperature)^theta1_e)"
    mu_L_h="(mu_max_h * (300 / Temperature)^theta1_h)"
    CreateNodeModel(device, region, "mu_L_e", mu_L_e)
    CreateNodeModel(device, region, "mu_L_h", mu_L_h)
    CreateNodeModelDerivative(device, region, "mu_L_e", mu_L_e, "Temperature")
    CreateNodeModelDerivative(device, region, "mu_L_h", mu_L_h, "Temperature")

    mu_e_N="(mu_max_e * mu_max_e / (mu_max_e - mu_min_e) * (Temperature/300)^(3*alpha_1_e - 1.5))"
    mu_h_N="(mu_max_h * mu_max_h / (mu_max_h - mu_min_h) * (Temperature/300)^(3*alpha_1_h - 1.5))"
    CreateNodeModel(device, region, "mu_e_N", mu_e_N)
    CreateNodeModel(device, region, "mu_h_N", mu_h_N)
    CreateNodeModelDerivative(device, region, "mu_e_N", mu_e_N, "Temperature")
    CreateNodeModelDerivative(device, region, "mu_h_N", mu_h_N, "Temperature")


    mu_e_c="(mu_min_e * mu_max_e / (mu_max_e - mu_min_e)) * (300/Temperature)^(0.5)"
    mu_h_c="(mu_min_h * mu_max_h / (mu_max_h - mu_min_h)) * (300/Temperature)^(0.5)"
    CreateNodeModel(device, region, "mu_e_c", mu_e_c)
    CreateNodeModel(device, region, "mu_h_c", mu_h_c)
    CreateNodeModelDerivative(device, region, "mu_e_c", mu_e_c, "Temperature")
    CreateNodeModelDerivative(device, region, "mu_h_c", mu_h_c, "Temperature")

    # PBH_e="(1.36e20/(Electrons + Holes) * (m_e) * (Temperature/300)^2)"
    # PBH_h="(1.36e20/(Electrons + Holes) * (m_h) * (Temperature/300)^2)"
    PBH_e="(1.36e20/(Electrons) * (m_e) * (Temperature/300)^2)"
    PBH_h="(1.36e20/(Holes) * (m_h) * (Temperature/300)^2)"
    CreateNodeModel (device,  region, "PBH_e", PBH_e)
    CreateNodeModelDerivative(device, region, "PBH_e", PBH_e, "Electrons", "Holes")
    CreateNodeModel (device,  region, "PBH_h", PBH_h)
    CreateNodeModelDerivative(device, region, "PBH_h", PBH_h, "Electrons", "Holes")
    CreateNodeModelDerivative(device, region, "PBH_e", PBH_e, "Temperature")
    CreateNodeModelDerivative(device, region, "PBH_h", PBH_h, "Temperature")

    # Z_D="(1 + 1 / (c_D + (Nref_D / max(Donors,1e-20))^2))"
    # Z_A="(1 + 1 / (c_A + (Nref_A / max(Acceptors,1e-20))^2))"
    Z_D="(1 + 1 / (c_D + (Nref_D / Donors)^2))"
    Z_A="(1 + 1 / (c_A + (Nref_A / Acceptors)^2))"

    CreateNodeModel (device, region, "MZ_D", Z_D)
    CreateNodeModel (device, region, "MZ_A", Z_A)

    ##It found that the name of the nodemodel can't begin with Z, this will make some error in the VISIT visualization
    ##So let's make a deal, when the name is begin with z, add a M to the begining

    N_D="(MZ_D * Donors)"
    N_A="(MZ_A * Acceptors)"
    CreateNodeModel (device,  region, "N_D", N_D)
    CreateNodeModel (device,  region, "N_A", N_A)

    N_e_sc="(N_D + N_A + Holes)"
    N_h_sc="(N_A + N_D + Electrons)"
    CreateNodeModel (device,  region, "N_e_sc", N_e_sc)
    CreateNodeModelDerivative(device, region, "N_e_sc", N_e_sc, "Electrons", "Holes")
    CreateNodeModel (device,  region, "N_h_sc", N_h_sc)
    CreateNodeModelDerivative(device, region, "N_h_sc", N_h_sc, "Electrons", "Holes")

    # PCW_e="(3.97e13 * (1/(MZ_D^3 * N_e_sc) * (Temperature/300)^3)^(2/3))"
    # PCW_h="(3.97e13 * (1/(MZ_A^3 * N_h_sc) * (Temperature/300)^3)^(2/3))"
    PCW_e="(3.97e13 * (1/(MZ_D^3 * Donors) * (Temperature/300)^3)^(2/3))"
    PCW_h="(3.97e13 * (1/(MZ_A^3 * Acceptors) * (Temperature/300)^3)^(2/3))"
    CreateNodeModel (device,  region, "PCW_e", PCW_e)
    CreateNodeModelDerivative(device, region, "PCW_e", PCW_e, "Electrons", "Holes")
    CreateNodeModel (device,  region, "PCW_h", PCW_h)
    CreateNodeModelDerivative(device, region, "PCW_h", PCW_h, "Electrons", "Holes")
    CreateNodeModelDerivative(device, region, "PCW_e", PCW_e, "Temperature")
    CreateNodeModelDerivative(device, region, "PCW_h", PCW_h, "Temperature")


    Pe="(1/(f_CW / PCW_e + f_BH/PBH_e))"
    Ph="(1/(f_CW / PCW_h + f_BH/PBH_h))"
    CreateNodeModel (device,  region, "Pe", Pe)
    CreateNodeModelDerivative(device, region, "Pe", Pe, "Electrons", "Holes")
    CreateNodeModel (device,  region, "Ph", Ph)
    CreateNodeModelDerivative(device, region, "Ph", Ph, "Electrons", "Holes")
    CreateNodeModelDerivative(device, region, "Pe", Pe, "Temperature")
    CreateNodeModelDerivative(device, region, "Ph", Ph, "Temperature")

    G_Pe="(1 - s1 / (s2 + (Temperature/(m_e*300))^s4 * Pe)^s3 + s5/((m_e * 300/Temperature)^s7*Pe)^s6)"
    G_Ph="(1 - s1 / (s2 + (Temperature/(m_h*300))^s4 * Ph)^s3 + s5/((m_h * 300/Temperature)^s7*Ph)^s6)"

    CreateNodeModel (device,  region, "G_Pe", G_Pe)
    CreateNodeModelDerivative(device, region, "G_Pe", G_Pe, "Electrons", "Holes")
    CreateNodeModel (device,  region, "G_Ph", G_Ph)
    CreateNodeModelDerivative(device, region, "G_Ph", G_Ph, "Electrons", "Holes")
    CreateNodeModelDerivative(device, region, "G_Pe", G_Pe, "Temperature")
    CreateNodeModelDerivative(device, region, "G_Ph", G_Ph, "Temperature")

    F_Pe="((r1 * Pe^r6 + r2 + r3 * m_e/m_h)/(Pe^r6 + r4 + r5 * m_e/m_h))"
    F_Ph="((r1 * Ph^r6 + r2 + r3 * m_h/m_e)/(Ph^r6 + r4 + r5 * m_h/m_e))"

    CreateNodeModel (device,  region, "F_Pe", F_Pe)
    CreateNodeModelDerivative(device, region, "F_Pe", F_Pe, "Electrons", "Holes")
    CreateNodeModel (device,  region, "F_Ph", F_Ph)
    CreateNodeModelDerivative(device, region, "F_Ph", F_Ph, "Electrons", "Holes")
    CreateNodeModelDerivative(device, region, "F_Pe", F_Pe, "Temperature")
    CreateNodeModelDerivative(device, region, "F_Ph", F_Ph, "Temperature")

    
    N_e_sc_eff="(N_D + G_Pe * N_A + Holes / F_Pe)"
    N_h_sc_eff="(N_A + G_Ph * N_D + Electrons / F_Ph)"

    CreateNodeModel (device,  region, "N_e_sc_eff", N_e_sc_eff)
    CreateNodeModelDerivative(device, region, "N_e_sc_eff", N_e_sc_eff, "Electrons", "Holes")
    CreateNodeModel (device,  region, "N_h_sc_eff", N_h_sc_eff)
    CreateNodeModelDerivative(device, region, "N_h_sc_eff", N_h_sc_eff, "Electrons", "Holes")
    CreateNodeModelDerivative(device, region, "N_e_sc_eff", N_e_sc_eff, "Temperature")
    CreateNodeModelDerivative(device, region, "N_h_sc_eff", N_h_sc_eff, "Temperature")

    mu_e_D_A_h="mu_e_N * N_e_sc/N_e_sc_eff * (Nref_1_e / N_e_sc)^alpha_1_e + mu_e_c * ((Electrons + Holes)/N_e_sc_eff)"
    mu_h_D_A_e="mu_h_N * N_h_sc/N_h_sc_eff * (Nref_1_h / N_h_sc)^alpha_1_h + mu_h_c * ((Electrons + Holes)/N_h_sc_eff)"

    CreateNodeModel          (device, region, "mu_e_D_A_h", mu_e_D_A_h)
    CreateNodeModelDerivative(device, region, "mu_e_D_A_h", mu_e_D_A_h, "Electrons", "Holes")

    CreateNodeModel          (device, region, "mu_h_D_A_e", mu_h_D_A_e)
    CreateNodeModelDerivative(device, region, "mu_h_D_A_e", mu_h_D_A_e, "Electrons", "Holes")
    CreateNodeModelDerivative(device, region, "mu_e_D_A_h", mu_e_D_A_h, "Temperature")
    CreateNodeModelDerivative(device, region, "mu_h_D_A_e", mu_h_D_A_e, "Temperature")


    mu_bulk_e_Node="mu_e_D_A_h * mu_L_e / (mu_e_D_A_h + mu_L_e)"
    CreateNodeModel          (device, region, "mu_bulk_e_Node", mu_bulk_e_Node)
    CreateNodeModelDerivative(device, region, "mu_bulk_e_Node", mu_bulk_e_Node, "Electrons", "Holes")

    mu_bulk_h_Node="mu_h_D_A_e * mu_L_h / (mu_h_D_A_e + mu_L_h)"
    CreateNodeModel          (device, region, "mu_bulk_h_Node", mu_bulk_h_Node)
    CreateNodeModelDerivative(device, region, "mu_bulk_h_Node", mu_bulk_h_Node, "Electrons", "Holes")
    CreateNodeModelDerivative(device, region, "mu_bulk_e_Node", mu_bulk_e_Node, "Temperature")
    CreateNodeModelDerivative(device, region, "mu_bulk_h_Node", mu_bulk_h_Node, "Temperature")


    CreateGeometricMean          (device, region, "mu_bulk_e_Node", "mu_bulk_e")
    CreateGeometricMeanDerivative(device, region, "mu_bulk_e_Node", "mu_bulk_e", "Electrons", "Holes", "Temperature")

    CreateGeometricMean          (device, region, "mu_bulk_h_Node", "mu_bulk_h")
    CreateGeometricMeanDerivative(device, region, "mu_bulk_h_Node", "mu_bulk_h", "Electrons", "Holes", "Temperature")

##introduce the velocity saturation

def Philips_VelocitySaturation(device, region, mu_sat, mu_bulk, eparallel, vsat):
    '''
    mu_sat  is the saturation velocity model name
    mu_bulk is the original bulk model on the edge
    eparallel is the name of the parallel efield model
    vsat is the name of the parameter for velocity saturation
    '''
    #expr = "2*({mu_bulk}) / (1 + (1 + 4 * (({mu_bulk} * {eparallel})/{vsat})^2)^0.5)"
    #expr = "2*({mu_bulk}) / (1 + (1 + 4 * (max(({mu_bulk}) * {eparallel}, 0)/{vsat})^2)^0.5)"
    #expr_hf = expr.format(mu_bulk=mu_bulk, eparallel=eparallel, vsat=vsat, jmag=jmag, jemin=jemin)
    #expr_lf = expr.format(mu_bulk=mu_bulk, eparallel="ElectricField", vsat=vsat, jmag=jmag, jemin=jemin)
    #mu="ifelse(abs({jmag}) > {jemin}, {expr_hf}, {expr_lf})".format(jmag=jmag, jemin=jemin, expr_lf=expr_lf, expr_hf=expr_hf)
    #mu="ifelse(abs({jmag}) <= {jemin}, {mu_bulk}, 2*({mu_bulk}) / (1 + (1 + 4 * (max(({mu_bulk}) * {eparallel}, 0)/{vsat})^2)^0.5))".format(mu_bulk=mu_bulk, eparallel=eparallel, vsat=vsat, jmag=jmag, jemin=jemin)
    mu="2*({mu_bulk}) / (1 + (1 + 4 * ((({mu_bulk}) * {eparallel})/{vsat})^2)^0.5)".format(mu_bulk=mu_bulk, eparallel=eparallel, vsat=vsat)
    #mu="ifelse({jmag} <= {jemin}, {mu_bulk}, 2*({mu_bulk}) / (1 + (1 + 4 * nodiff(max(({mu_bulk}) * {eparallel}, 0)/{vsat})^2)^0.5))".format(mu_bulk=mu_bulk, eparallel=eparallel, vsat=vsat, jmag=jmag, jemin=jemin)
    CreateEdgeModel           (device,  region, mu_sat, mu)
    CreateEdgeModelDerivatives (device,  region, mu_sat, mu, "Potential")
    CreateEdgeModelDerivatives (device,  region, mu_sat, mu, "Electrons")
    CreateEdgeModelDerivatives (device,  region, mu_sat, mu, "Holes")
    CreateEdgeModelDerivatives (device,  region, mu_sat, mu, "Temperature")

def CaugheyThoms_electronsVelocitySaturation(device, region, mu_sat, mu_bulk, eparallel):
    vsate_Temperaturedependence="2.4e7/(1+0.8*exp( edgeTemp / 600))"
    CreateEdgeModel(device,  region, "vsate_Temperaturedependence", vsate_Temperaturedependence)
    CreateEdgeModelDerivatives (device,  region, "vsate_Temperaturedependence", vsate_Temperaturedependence, "Temperature")

    mu="({mu_bulk})*(1 / (1+((({mu_bulk}) * {eparallel})/vsate_Temperaturedependence)^2))^0.5".format(mu_bulk=mu_bulk, eparallel=eparallel)
    CreateEdgeModel(device,  region, mu_sat, mu)
    CreateEdgeModelDerivatives (device,  region, mu_sat, mu, "Potential")
    CreateEdgeModelDerivatives (device,  region, mu_sat, mu, "Electrons")
    CreateEdgeModelDerivatives (device,  region, mu_sat, mu, "Holes")
    CreateEdgeModelDerivatives (device,  region, mu_sat, mu, "Temperature")

def CaugheyThoms_holesVelocitySaturation(device, region, mu_sat, mu_bulk, eparallel):
    vsath_Temperaturedependence="2.4e7/(1+0.8*exp( edgeTemp / 600))"
    CreateEdgeModel           (device,  region, "vsath_Temperaturedependence", vsath_Temperaturedependence)
    CreateEdgeModelDerivatives (device,  region, "vsath_Temperaturedependence", vsath_Temperaturedependence, "Temperature")
    mu="({mu_bulk})*(1 / (1+((({mu_bulk}) * {eparallel})/vsath_Temperaturedependence)))".format(mu_bulk=mu_bulk, eparallel=eparallel)
    CreateEdgeModel           (device,  region, mu_sat, mu)
    CreateEdgeModelDerivatives (device,  region, mu_sat, mu, "Potential")
    CreateEdgeModelDerivatives (device,  region, mu_sat, mu, "Electrons")
    CreateEdgeModelDerivatives (device,  region, mu_sat, mu, "Holes")
    CreateEdgeModelDerivatives (device,  region, mu_sat, mu, "Temperature")


def CreateAlbrechtMobilityModels(device,region):
    edge_from_node_model(device=device, region=region, node_model="NetDoping")
    edge_average_model(device=device, region=region,edge_model="edgeNetDoping", node_model="NetDoping")
    inverse_mu_n_Albrecht = "2.61e-4*edgeNetDoping/1e17*(edgeTemp/300)^(-3/2)*log(1+3*(edgeTemp/300)^2*(edgeNetDoping/1e17)^(-2/3)) + \
                                2.9e-4*(edgeTemp/300)^(3/2) + \
                                     170e-4/(exp(1065/edgeTemp)-1) "
    # mu_n_Albrecht="1/inverse_mu_n_Albrecht"

    mu_n_Albrecht="1200"
    mu_p_Albrecht = "holes_mobility"
    #contribution of holes is little, set as constant
    for i in (
        ("inverse_mu_n_Albrecht", inverse_mu_n_Albrecht),
        ("mu_n_Albrecht", mu_n_Albrecht),
        ("mu_p_Albrecht", mu_p_Albrecht),
    ):
        n = i[0]
        e = i[1]
        CreateEdgeModel(device, region, n, e)
        CreateEdgeModelDerivatives(device, region,n,e,"Temperature")

#only for electrons,and is not good for convergence
#no need to consider the saturation for holes because its little contribution
def CaugheyGaNMaziarVelocitySaturation(device, region, mu_sat, mu_bulk, eparallel):

    mu="(({mu_bulk})+1.1299e7*(({eparallel})^(5.2217-1)/(3.9552e5)^5.2217))/ \
        (1+3.024*(({eparallel})/3.9552e5)^1.0269+(({eparallel})/3.9552e5)^5.2217)".format(mu_bulk=mu_bulk, eparallel=eparallel)
    CreateEdgeModel           (device,  region, mu_sat, mu)
    CreateEdgeModelDerivatives (device,  region, mu_sat, mu, "Potential")
    CreateEdgeModelDerivatives (device,  region, mu_sat, mu, "Electrons")
    CreateEdgeModelDerivatives (device,  region, mu_sat, mu, "Holes")
    CreateEdgeModelDerivatives (device,  region, mu_sat, mu, "Temperature")



def CreateThermoelectricPowers(device,region):
    Nc = "2.8e19*(Temperature/300)^1.5"
    Nv = "1.04e19*(Temperature/300)^1.5"
    #态密度

    CreateNodeModel           (device,  region, "nstate_density", Nc)
    CreateNodeModel           (device,  region, "pstate_density", Nv)
    CreateNodeModelDerivative (device,  region, "nstate_density", Nc, "Temperature")
    CreateNodeModelDerivative (device,  region, "pstate_density", Nv, "Temperature")


    Pn = "-boltzmannConstant/(ElectronCharge)*(1.5+log(nstate_density/Electrons))"
    Pp = "-boltzmannConstant/(ElectronCharge)*(1.5+log(pstate_density/Holes))" 
    CreateNodeModel           (device,  region, "nThermpower", Pn)
    CreateNodeModel           (device,  region, "pThermpower", Pp)    
    CreateNodeModelDerivative (device,  region, "nThermpower", Pn, "Electrons")
    CreateNodeModelDerivative (device,  region, "pThermpower", Pp, "Holes")
    CreateNodeModelDerivative (device,  region, "nThermpower", Pn, "Temperature")
    CreateNodeModelDerivative (device,  region, "pThermpower", Pp, "Temperature")


def CreateJouleHeatSilicon(device,region):
    edge_average_model(device=device, region=region, edge_model="edgeElectrons", node_model="Electrons")
    edge_average_model(device=device, region=region, edge_model="edgeHoles", node_model="Holes")
    JouleHeat ="ElectronCurrent*ElectronCurrent/(mu_vsat_e*ElectronCharge*edgeElectrons)+ HoleCurrent*HoleCurrent/(mu_bulk_h*ElectronCharge*edgeHoles)"
    #Heatsource ="ElectronCurrent*ElectronCurrent/(mu_nT*ElectronCharge*edgeElectrons)+ HoleCurrent*HoleCurrent/(mu_pT*ElectronCharge*edgeHoles)"
    # Heatsource ="ElectricField * (ElectronCurrent + HoleCurrent)"
    CreateEdgeModel(device, region, "JouleHeat", JouleHeat)
    for v in ("Potential", "Electrons", "Holes", "Temperature"):
        CreateEdgeModelDerivatives(device, region, "JouleHeat", JouleHeat, v)

def CreateJouleHeatOxide(device,region):
    Heatsource ="0"
    CreateEdgeModel(device, region, "Heatsource", Heatsource)


def CreateRecombinationHeatSilicon(device,region):
    if not InEdgeModelList(device, region, "edgeElectrons"):
        edge_average_model(device=device, region=region, edge_model="edgeElectrons", node_model="Electrons")
    if not InEdgeModelList(device, region, "edgeHoles"):
        edge_average_model(device=device, region=region, edge_model="edgeHoles", node_model="Holes")
    if not InEdgeModelList(device, region, "edgePotential"):
        edge_average_model(device=device, region=region, edge_model="edgePotential", node_model="Potential")     
      
    quasiPotential_n= "edgePotential - edgeV_t_T*log(edgeElectrons/n_i)"
    quasiPotential_p= "edgePotential + edgeV_t_T*log(edgeHoles/n_i)"
    CreateEdgeModel(device,  region, "quasiPotential_n", quasiPotential_n)
    CreateEdgeModel(device,  region, "quasiPotential_p", quasiPotential_p) 
    for v in ("Potential", "Electrons", "Holes", "Temperature"):
        CreateEdgeModelDerivatives(device, region, "quasiPotential_n", quasiPotential_n, v)
        CreateEdgeModelDerivatives(device, region, "quasiPotential_p", quasiPotential_p, v)

    edge_average_model(device=device, region=region, edge_model="edgeUSRH", node_model="USRH")  
    edge_average_model(device=device, region=region, edge_model="edgenThermpower", node_model="nThermpower")  
    edge_average_model(device=device, region=region, edge_model="edgepThermpower", node_model="pThermpower")  

    # Recomheat = "ElectronCharge*edgeUSRH*(quasiPotential_p-quasiPotential_n+edgeTemp*(edgepThermpower-edgenThermpower))"
    Recomheat = "edgeUSRH*(1.17*ElectronCharge+3*boltzmannConstant*edgeTemp)"

    CreateEdgeModel(device, region, "Recomheat", Recomheat)    
    for v in ("Potential", "Electrons", "Holes", "Temperature"):
        CreateEdgeModelDerivatives(device, region, "Recomheat", Recomheat, v)

def CreatePJTheatSilicon(device,region):
    edge_average_model(device=device, region=region, edge_model="grad_nThermpower", node_model="nThermpower", average_type="gradient")
    edge_average_model(device=device, region=region, edge_model="grad_pThermpower", node_model="pThermpower", average_type="gradient")


    if not InEdgeModelList(device, region, "edgeTemp"):
        edge_average_model(device=device, region=region, edge_model="edgeTemp", node_model="Temperature")

    PJTheat = "-edgeTemp*(ElectronCurrent*grad_nThermpower+HoleCurrent*grad_pThermpower)"
    CreateEdgeModel(device, region, "PJTheat", PJTheat)
    PJT_pos="(PJTheat+abs(PJTheat))/2"
    PJT_neg="(PJTheat-abs(PJTheat))/2"
    # PJT_pos="ifelse(PJTheat <= 0, PJTheat/1e10, PJTheat)"
    # PJT_neg="ifelse(PJTheat <= 0, PJTheat, PJTheat)"    
    
    CreateEdgeModel(device, region, "PJT_pos", PJT_pos)
    CreateEdgeModel(device, region, "PJT_neg", PJT_neg)



def CreateHeatsourceSilicon(device,region):
    Heatsource = "JouleHeat"
    CreateEdgeModel(device, region, "Heatsource", Heatsource)

    for v in ("Potential", "Electrons", "Holes", "Temperature"):
        CreateEdgeModelDerivatives(device, region, "Heatsource", Heatsource, v)

def CreateHeatsourceSiliconWithPT(device,region):
    Heatsource = "JouleHeat"
    CreateEdgeModel(device, region, "Heatsource", Heatsource)

    for v in ("Potential", "Electrons", "Holes", "Temperature"):
        CreateEdgeModelDerivatives(device, region, "Heatsource", Heatsource, v)


def CreateTemperaturefield(device,region):
    if not InNodeModelList(device, region, "Temperature"):
        print("Creating Node Solution Temperature")
        CreateSolution(device, region, "Temperature")

    edge_average_model(device=device, region=region, edge_model="TemperatureField", node_model="Temperature", average_type="gradient")
    edge_average_model(device=device, region=region, edge_model="TemperatureField", node_model="Temperature", average_type="gradient", derivative="Temperature")
    for i in (
      ("TemperatureEdgeFlux", "thermalconductivity * TemperatureField"),
    ):
        n = i[0]
        e = i[1]
        CreateEdgeModel(device, region, n, e)
        CreateEdgeModelDerivatives(device, region, n, e, "Temperature")

    equation(device=device, region=region, name="HeatconductionEquation", variable_name="Temperature",
            edge_volume_model="Heatsource", edge_model="TemperatureEdgeFlux",
            time_node_model="", variable_update="positive")


def CreateTemperatureboundary(device,region,contact):
    contact_model = "Temperature -{0}".format(GetContactTemperatureName(contact))

    contact_model_name = GetContactTemperatureNodeModelName(contact)
    CreateContactNodeModel(device, contact, contact_model_name, contact_model)
    # Simplify it too complicated
    CreateContactNodeModel(device, contact, "{0}:{1}".format(contact_model_name,"Temperature"), "1")

    contact_equation(device=device, contact=contact, name="HeatconductionEquation",
                        node_model=contact_model_name, edge_model="")


def CreateHeatfluxboundary(device,region,contact,flux):
    ##Createtemperatureboundarycondition
    CreateContactNodeModel(device, contact,"boundaryflux" , -flux)
    CreateContactNodeModel(device, contact, "boundaryflux:Temperature","0")
    contact_equation(device=device, contact=contact, name="HeatconductionEquation",
                    node_model="boundaryflux", edge_model="TemperatureEdgeFlux")
    

# def CreatePotentialboundary(device,region,contact):
#     ##Createtemperatureboundarycondition
#     CreateContactNodeModel(device, contact,"boundaryElectricField" , "0")
#     CreateContactNodeModel(device, contact, "boundaryElectricField:Potential","0")
#     contact_equation(device=device, contact=contact, name="PotentialEquation",
#                     node_model="boundaryElectricField", edge_model="PotentialEdgeFlux")

# def CreateElectronsboundary(device,region,contact):
#     ##Createtemperatureboundarycondition
#     CreateContactNodeModel(device, contact,"boundaryElectrons" , "0")
#     CreateContactNodeModel(device, contact, "boundaryElectrons:Electrons","0")
#     CreateEdgeModel(device, region, "gradElectrons", "(Electrons@n0-Electrons@n1)*EdgeInverseLength")
#     CreateEdgeModelDerivatives(device, region,"gradElectrons", "(Electrons@n0-Electrons@n1)*EdgeInverseLength", "Electrons")
#     contact_equation(device=device, contact=contact, name="PotentialEquation",
#                     node_model="boundaryElectrons", edge_model="gradElectrons")

# def CreateHolesboundary(device,region,contact):
#     ##Createtemperatureboundarycondition
#     CreateContactNodeModel(device, contact,"boundaryHoles" , "0")
#     CreateContactNodeModel(device, contact, "boundaryHoles:Holes","0")
#     CreateEdgeModel(device, region, "gradHoles", "(Holes@n0-Holes@n1)*EdgeInverseLength")
#     CreateEdgeModelDerivatives(device, region,"gradHoles", "(Holes@n0-Holes@n1)*EdgeInverseLength", "Holes")
#     contact_equation(device=device, contact=contact, name="PotentialEquation",
#                     node_model="boundaryHoles", edge_model="gradHoles")
       



def CreateTemperatureinterface(device,interface):

    '''
      continuous Temperature at interface
    '''
    model_name = CreateContinuousInterfaceModel(device, interface, "Temperature")
    interface_equation(device=device, interface=interface, name="HeatconductionEquation", interface_model=model_name, type="continuous")



