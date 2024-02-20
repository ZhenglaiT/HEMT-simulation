from devsim import *
from devsim.python_packages.simple_dd import *
from devsim.python_packages.simple_physics import *
from Band_structure_dd import *

contactcharge_node="contactcharge_node"
contactcharge_edge="contactcharge_edge"
ece_name="ElectronContinuityEquation"
hce_name="HoleContinuityEquation"
celec_model = "(1e-10 + 0.5*abs(NetDoping+(NetDoping^2 + 4 * n_i^2)^(0.5)))"
chole_model = "(1e-10 + 0.5*abs(-NetDoping+(NetDoping^2 + 4 * n_i^2)^(0.5)))"

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
ns_GaN=-1.45e13        #正号代表总电荷为正，负号代表总电荷为负

ns=1.6785e13   #没乘0.8的



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
   


def SetAlGaNParameters(device,region,xcomp):
    set_parameter(device=device, region=region, name="xcomp",   value=xcomp)
    set_parameter(device=device, region=region, name="Permittivity",   value=(8.5*xcomp+8.9*(1-xcomp))*eps_0)
    set_parameter(device=device, region=region, name="Eg",   value=3.97265)
    set_parameter(device=device, region=region, name="Affinity",   value=3.93)
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
    set_parameter(device=device, region=region, name="holes_mobility",     value=100) #The default parameter in Silvaco
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
    set_parameter(device=device, region=region, name="V_surfe", value=2.1723e7) 
    set_parameter(device=device, region=region, name="V_surfh", value=8.9591e5)

    # set_parameter(device=device, region=region, name="V_surfe", value=2.573e6)     
    # set_parameter(device=device, region=region, name="V_surfh", value=1.93e6)

    set_parameter(device=device, region=region, name="Electrons_eq", value=1.2e-8) 
    set_parameter(device=device, region=region, name="Holes_eq", value=8.163e-22)

def SetGaNParameters(device,region):
    set_parameter(device=device, region=region, name="Permittivity",   value=8.9*eps_0)
    set_parameter(device=device, region=region, name="Eg",   value=3.4346)
    set_parameter(device=device, region=region, name="Affinity",   value=4.31) #eV
    set_parameter(device=device, region=region, name="DOS_emasses",   value=0.2)
    set_parameter(device=device, region=region, name="DOS_hmasses",   value=1.0)
    set_parameter(device=device, region=region, name="ElectronCharge", value=q)
    set_parameter(device=device, region=region, name="boltzmannConstant",              value=k)  
    set_parameter(device=device, region=region, name="pi",              value=Pi)    
    set_parameter(device=device, region=region, name="planckconstant",              value=h)    
    set_parameter(device=device, region=region, name="n_i",            value=n_i_GaN)
    set_parameter(device=device, region=region, name="GaNthermalconductivity",     value=GaNthermalconductivity)
    set_parameter(device=device, region=region, name="holes_mobility",     value=8) #The default parameter in Silvaco
    set_parameter(device=device, region=region, name="n1", value=n_i_GaN)
    set_parameter(device=device, region=region, name="p1", value=n_i_GaN)
    set_parameter(device=device, region=region, name="taun", value=1.2e-8) 
    set_parameter(device=device, region=region, name="taup", value=1.2e-8) 
    # set_parameter(device=device, region=region, name="surface_charge", value=ns) 
    set_parameter(device=device, region=region, name="An", value=24.0351) #Richardson constants
    set_parameter(device=device, region=region, name="Ap", value=120.1753) #Richardson constants
    set_parameter(device=device, region=region, name="V_surfe", value=2.7527e7) 
    set_parameter(device=device, region=region, name="V_surfh", value=6.7147e5)
    # set_parameter(device=device, region=region, name="V_surfe", value=2.573e6)     
    # set_parameter(device=device, region=region, name="V_surfh", value=1.93e6)
    set_parameter(device=device, region=region, name="Electrons_eq", value=5.4308e24) 
    set_parameter(device=device, region=region, name="Holes_eq", value=2.0773e-45)


def CreatAlGaNPotentialOnly(device, region):
#   The sheet charge is determined by the polarization and 2DEG.
#   For AlGaN, positive charge is accumlated in the bottom due to polarization.
    set_parameter(device=device,region=region,name="surface_charge", value=0) 
    # set_parameter(device=device,region=region, name="surface_charge", value=0) 

    if not InNodeModelList(device, region, "Potential"):
        print("Creating Node Solution Potential")
        CreateSolution(device, region, "Potential")
    elec_i = "n_i*exp(Potential/V_t_T)"
    hole_i = "n_i^2/IntrinsicElectrons"
 

    # charge_i = "kahan3(IntrinsicHoles, -IntrinsicElectrons, NetDoping)"
    Sheet_charge="surface_charge * SurfaceArea / NodeVolume"
    # Sheet_charge="surface_charge*exp(-abs(y-2.5e-6)/1e-6)"
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

    set_parameter(device=device,region=region, name="surface_charge", value=ns) 
    # set_parameter(device=device,region=region, name="surface_charge", value=ns) 
    
    if not InNodeModelList(device, region, "Potential"):
        print("Creating Node Solution Potential")
        CreateSolution(device, region, "Potential")
    elec_i = "n_i*exp(Potential/V_t_T)"
    hole_i = "n_i^2/IntrinsicElectrons"
    # charge_i = "kahan3(IntrinsicHoles, -IntrinsicElectrons, NetDoping)"
    Sheet_charge="surface_charge * SurfaceArea / NodeVolume"
    # Sheet_charge="surface_charge*exp(-abs(y-2.5e-6)/1e-6)"
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


def CreateECE_BS(device, region,mu_n):
    CreateElectronCurrent_BS(device, region,mu_n)
    NCharge = "-ElectronCharge * Electrons"
    CreateNodeModel(device, region, "NCharge", NCharge)
    CreateNodeModelDerivative(device, region, "NCharge", NCharge, "Electrons")

    equation(device=device, region=region, name="ElectronContinuityEquation", variable_name="Electrons",
            time_node_model = "NCharge",
            edge_model="ElectronCurrent", variable_update="positive", node_model="ElectronGeneration")

def CreateHCE_BS(device, region,mu_p):
    CreateHoleCurrent_BS(device, region, mu_p)
    PCharge = "ElectronCharge * Holes"
    CreateNodeModel(device, region, "PCharge", PCharge)
    CreateNodeModelDerivative(device, region, "PCharge", PCharge, "Holes")

    equation(device=device, region=region, name="HoleContinuityEquation", variable_name="Holes",
             time_node_model = "PCharge",
             edge_model="HoleCurrent", variable_update="positive", node_model="HoleGeneration")


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

def CreateReferencePotential(device,region):
    #use the parameter of one region (Such as AlGaN or GaN) to create the reference potential
    PotentialRefer="3.93+boltzmannConstant*Temperature/ElectronCharge*log(2.8445e18/3.126e-15)"
    CreateNodeModel(device,  region, "Potential_Rerfer", PotentialRefer)

def CreateBandStructure(device,region):
    Ec="Potential_Rerfer-Potential-Affinity"
    Ev="Potential_Rerfer-Potential-Affinity-Eg"

    CreateNodeModel(device,  region, "Conduction_energy", Ec)
    CreateNodeModel(device,  region, "Valence_energy", Ev)

def CreateDriftDiffusionAtSchottkyContact(device,region,contact,workfun):
   
    # V_surfe = "110*300^2/ElectronCharge/(2.8445e18)"
    # V_surfh = "30*300^2/ElectronCharge/(1.8810e19)"

    # CreateContactNodeModel           (device,  contact, "V_surfe", V_surfe)
    # CreateContactNodeModel           (device,  contact, "V_surfh", V_surfh)
    # CreateContactNodeModelDerivative (device,  contact, "V_surfe", V_surfe, "Temperature")
    # CreateContactNodeModelDerivative (device,  contact, "V_surfh", V_surfh, "Temperature")    
    
    # deltaPhi="Affinity + Eg/2 + \
    #                 boltzmannConstant*Temperature/(2*ElectronCharge)*log(nstate_density/pstate_density) - \
    #                 {0}".format(workfun)
    
    # CreateContactNodeModel           (device,  contact, "deltaPhi", deltaPhi)
    # CreateContactNodeModelDerivative (device,  contact, "deltaPhi", deltaPhi, "Temperature")

    # Electrons_eq = "n_i*exp(ElectronCharge*(deltaPhi)/boltzmannConstant/Temperature)"
    # Holes_eq =  "n_i*exp(-ElectronCharge*(deltaPhi)/boltzmannConstant/Temperature)"

    # #refered to the Palankovski. "Analysis and simulation of heterostructure devices"
    # Electrons_eq = "nstate_density*exp(-ElectronCharge*({0}-Affinity)/boltzmannConstant/Temperature)".format(workfun)
    # Holes_eq =  "pstate_density*exp(ElectronCharge*({0}-Affinity-Eg)/boltzmannConstant/Temperature)".format(workfun)


    # CreateContactNodeModel           (device,  contact, "Electrons_eq", Electrons_eq)
    # CreateContactNodeModel           (device,  contact, "Holes_eq", Holes_eq)
    # CreateContactNodeModelDerivative (device,  contact, "Electrons_eq", Electrons_eq, "Temperature")
    # CreateContactNodeModelDerivative (device,  contact, "Holes_eq", Holes_eq, "Temperature")    

    

    contact_electrons_model ="ElectronCharge*V_surfe*(Electrons-Electrons_eq) * ContactSurfaceArea/NodeVolume"
    contact_holes_model     ="-ElectronCharge*V_surfh*(Holes-Holes_eq) * ContactSurfaceArea/NodeVolume"
   


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

    contact_equation(device=device, contact=contact, name="ElectronContinuityEquation", 
                    node_model=contact_electrons_name,  edge_current_model="ElectronCurrent")
    contact_equation(device=device, contact=contact, name="HoleContinuityEquation",
                    node_model=contact_holes_name,  edge_current_model="HoleCurrent")


def CreateDriftDiffusionAtSubstrate(device,region,contact):
   
    contact_electrons_model ="0"
    contact_holes_model     ="0"

    contact_electrons_name="{0}_electrons_Substrate".format(contact)
    contact_holes_name="{0}_holes_Substrate".format(contact)

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


def CreateHeterojunctionInterface(device, interface):
    #Make sure the continous of Potential
    model_name = CreateContinuousInterfaceModel(device, interface, "Potential")
    interface_equation(device=device, interface=interface, name="PotentialEquation", interface_model=model_name, type="continuous")


def CreateThermionicEmission(device,interface): 

    #create the thermionic emission in the Heterojunction interface.

    #the alignment of energy bands is determined by the affinity rule 

    # Therma_Em_expn="ifelse(  Eg@r0>Eg@r1,\
    #     (-An@r0*(Temperature@r0)^2/(nstate_density@r0)*Electrons@r0+\
    #     An@r0*(Temperature@r1)^2/(nstate_density@r1)*Electrons@r1*exp((-Delta_Ec*ElectronCharge)/(boltzmannConstant*Temperature@r1))),\
    #     (-An@r0*(Temperature@r0)^2/(nstate_density@r0)*Electrons@r0*exp((-Delta_Ec*ElectronCharge)/(boltzmannConstant*Temperature@r0))+\
    #     An@r1*(Temperature@r1)^2/(nstate_density@r1)*Electrons@r1)  )"
    # Delta_Ec="abs(Conduction_energy@r0-Conduction_energy@r1)"
    # Delta_Ev="abs(Valence_energy@r0-Valence_energy@r1)"

    # Therma_Em_expp="ifelse(  Eg@r0>Eg@r1,\
    #     (Ap@r0*(Temperature@r0)^2/(pstate_density@r0)*Holes@r0*exp((-Delta_Ev*ElectronCharge)/(boltzmannConstant*Temperature@r0))-\
    #     Ap@r1*(Temperature@r1)^2/(pstate_density@r1)*Holes@r1),\
    #     (Ap@r0*(Temperature@r0)^2/(pstate_density@r0)*Holes@r0-\
    #     Ap@r1*(Temperature@r1)^2/(pstate_density@r1)*Holes@r1*exp((-Delta_Ev*ElectronCharge)/(boltzmannConstant*Temperature@r1))) )"
   
    # set_parameter(device=device, name="Delta_Ec",   value=0.18)
    # set_parameter(device=device, name="Delta_Ev",   value=0.1614)
    # Therma_Em_expn="(-An@r0*(Temperature@r0)^2/(nstate_density@r0)*Electrons@r0+\
    #     An@r1*(Temperature@r1)^2/(nstate_density@r1)*Electrons@r1*exp((-Delta_Ec*ElectronCharge)/(boltzmannConstant*Temperature@r1)))"

    # Therma_Em_expp="(Ap@r0*(Temperature@r0)^2/(pstate_density@r0)*Holes@r0*exp((-Delta_Ev*ElectronCharge)/(boltzmannConstant*Temperature@r0))-\
    #     Ap@r1*(Temperature@r1)^2/(pstate_density@r1)*Holes@r1)"


    Therma_Em_expn="(-An@r0*(Temperature@r0)^2/(nstate_density@r0)*Electrons@r0+\
        An@r1*(Temperature@r1)^2/(nstate_density@r1)*Electrons@r1*exp((-abs(Conduction_energy@r0-Conduction_energy@r1)*ElectronCharge)/(boltzmannConstant*Temperature@r1)))"

    Therma_Em_expp="(Ap@r0*(Temperature@r0)^2/(pstate_density@r0)*Holes@r0*exp((-abs(Valence_energy@r0-Valence_energy@r1)*ElectronCharge)/(boltzmannConstant*Temperature@r0))-\
        Ap@r1*(Temperature@r1)^2/(pstate_density@r1)*Holes@r1)"


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

def CreateVaribleDerivatives(device,region):
    #This function is to reduce the noise in output
    #Not sure its influence on the simulation，consider the independence of varibles, it may be reasonable
    CreateSolution(device, region, "Temperature:Electrons")
    CreateSolution(device, region, "Temperature:Holes")
    CreateSolution(device, region, "nstate_density:Electrons")
    CreateSolution(device, region, "pstate_density:Holes")
    CreateSolution(device, region, "Conduction_energy:Electrons")
    CreateSolution(device, region, "Valence_energy:Holes")
    CreateSolution(device, region, "Conduction_energy:Temperature")
    CreateSolution(device, region, "Valence_energy:Temperature")
    CreateSolution(device, region, "Electrons:Temperature")
    CreateSolution(device, region, "Holes:Temperature")   
    # for i in (
    # ("Temperature:Electrons", 0),
    # ("Temperature:Holes", 0),   
    # ("nstate_density:Electrons", 0),
    # ("pstate_density:Electrons", 0),
    # ("nstate_density:Holes", 0),
    # ("pstate_density:Holes", 0),
    # ):
    #     n = i[0]
    #     e = i[1]
    #     CreateNodeModel(device, region, n, e)




def CreateAlbrechtMobilityModels(device,region):
    edge_from_node_model(device=device, region=region, node_model="NetDoping")
    edge_average_model(device=device, region=region,edge_model="edgeNetDoping", node_model="NetDoping")
    inverse_mu_n_Albrecht = "2.61e-4*edgeNetDoping/1e17*(edgeTemp/300)^(-3/2)*log(1+3*(edgeTemp/300)^2*(edgeNetDoping/1e17)^(-2/3)) + \
                                2.9e-4*(edgeTemp/300)^(3/2) + \
                                     170e-4/(exp(1065/edgeTemp)-1) "
    mu_n_Albrecht="1/inverse_mu_n_Albrecht"
    # mu_n_Albrecht="1200"
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

def CreateGaNThermalConductivity(device,region):
    Kappa_GaN = "GaNthermalconductivity*pow((edgeTemp)/(300),-0.28)"
    for i in (
        ("thermalconductivity", Kappa_GaN),
    ):
        n = i[0]
        e = i[1]
        CreateEdgeModel(device, region, n, e)
        CreateEdgeModelDerivatives(device, region,n,e,"Temperature")



def CreateAlGaNThermalConductivity(device,region,xcomp):
    Kappa_GaN_mid = "GaNthermalconductivity*pow((edgeTemp)/(300),-0.28)"
    Kappa_AlN_mid = "AlNthermalconductivity*pow((edgeTemp)/(300),-1.64)"   
    Kappa_AlGaN = "GaN_thermalconductivity*(1-{xcomp})+AlN_thermalconductivity*{xcomp}".format(xcomp=xcomp)     
    for i in (
        ("GaN_thermalconductivity", Kappa_GaN_mid),
        ("AlN_thermalconductivity", Kappa_AlN_mid),
        ("thermalconductivity", Kappa_AlGaN),       
    ):
        n = i[0]
        e = i[1]
        CreateEdgeModel(device, region, n, e)
        CreateEdgeModelDerivatives(device, region,n,e,"Temperature")
