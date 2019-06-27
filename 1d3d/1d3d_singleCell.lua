-------------------------------------------------------------------------
-- This script sets up a simple test for a 1d/3d hybrid simulation.    --
-- On the 1d domain, it solves the cable equation with HH channels,    --
-- activating specifically set synapses. On the 3d domain, it solves   --
-- a calcium problem (diffusion and buffering) with channels and pumps --
-- in the plasma membrane, where VDCCs are activated according to the  --
-- potential mapped from the 1d domain. Additionally, the 3d domain    --
-- contains an ER on whose membrane pumps and channels cause calcium   --
-- exchange with the cytosol.                                          --
--                                                                     --
-- author: mbreit                                                      --
-- date:   23-02-2017                                                  --
-------------------------------------------------------------------------

-- for profiler output
SetOutputProfileStats(false)

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

AssertPluginsLoaded({"cable_neuron", "neuro_collection"})

dim = 3
InitUG(dim, AlgebraType("CPU", 1));


-- choice of grid and refinement level
gridName1d = util.GetParam("-grid1d", "../apps/calciumDynamics_app/grids/smith_1d.ugx")
gridName3d = util.GetParam("-grid3d", "../apps/calciumDynamics_app/grids/smith_3d_farSyns.ugx")
numRefs = util.GetParamNumber("-numRefs", 0)
gridSyn = string.sub(gridName1d, 1, string.len(gridName1d) - 4) .. "_syns.ugx" 

-- parameters for instationary simulation
dt1d = util.GetParamNumber("-dt1d", 1e-5) -- in s
dt3d = util.GetParamNumber("-dt3d", 1e-2) -- in s
dt3dStart = util.GetParamNumber("-dt3dstart", dt3d)
endTime = util.GetParamNumber("-endTime", 1.0)  -- in s

-- with simulation of single ion concentrations?
withIons = util.HasParamOption("-ions")

-- choice of solver setup
solverID = util.GetParam("-solver", "GMG")
solverID = string.upper(solverID)
validSolverIDs = {}
validSolverIDs["GMG"] = 0;
validSolverIDs["GS"] = 0;
validSolverIDs["ILU"] = 0;
if (validSolverIDs[solverID] == nil) then
    error("Unknown solver identifier " .. solverID)
end

-- specify "-verbose" to output linear solver convergence
verbose1d = util.HasParamOption("-verbose1d")
verbose3d = util.HasParamOption("-verbose3d")

-- vtk output?
generateVTKoutput = util.HasParamOption("-vtk")
pstep = util.GetParamNumber("-pstep", dt3d, "plotting interval")

-- file handling
filename = util.GetParam("-outName", "hybrid_test")
filename = filename.."/"


-- choose length of time step at the beginning
-- if not timeStepStart = 2^(-n)*timeStep, take nearest lower number of that form
function log2(x)
	return math.log(x) / math.log(2)
end
startLv =  math.ceil(log2(dt3d / dt3dStart))
dt3dStartNew = dt3d / math.pow(2, startLv)
if (math.abs(dt3dStartNew - dt3dStart) / dt3dStart > 1e-5) then 
	print("dt3dStart argument ("..dt3dStart..") was not admissible; taking "..dt3dStartNew.." instead.")
end
dt3dStart = dt3dStartNew


print("Chosen parameters:")
print("    grid       = " .. gridName1d)
print("    numRefs    = " .. numRefs)
print("    dt1d       = " .. dt1d)
print("    dt3d       = " .. dt3d)
print("    dt3dStart  = " .. dt3dStart)
print("    endTime    = " .. endTime)
print("    pstep      = " .. pstep)
print("    ions       = " .. tostring(withIons))
print("    solver     = " .. solverID)
print("    verbose1d  = " .. tostring(verbose1d))
print("    verbose3d  = " .. tostring(verbose3d))
print("    vtk        = " .. tostring(generateVTKoutput))
print("    outname    = " .. filename)

------------------
-- set synapses --
------------------
-- firing pattern of the synapse
syn_onset = {0.0, 0.0, 0.0, 0.0, 0.0}
syn_tau = 4e-4
syn_gMax = 1.2e-9
syn_revPot = 0.0

function file_exists(name)
	local f = io.open(name,"r")
	if f ~= nil then io.close(f) return true else return false end
end

-- check if grid version with synapses already exists
-- if so, just use it, otherwise, create it
if not file_exists(gridSyn) then
	-- synapse distributor only works in serial mode
	if NumProcs() > 1 then
		print("Cannot use SynapseDistributor in parallel atm. Please create synapse geometry in serial.")
		exit()
	end
	
	-- use synapse distributor to place synapses on the grid
	synDistr = SynapseDistributor(gridName1d)

	syn1 = AlphaSynapsePair()
	syn1:set_id(0)
	syn1:set_onset(syn_onset[1])
	syn1:set_tau(syn_tau)                   -- default value
	syn1:set_gMax(syn_gMax)                 -- default value
	syn1:set_reversal_potential(syn_revPot) -- default value
	--synDistr:place_synapse_at_coords({-19.9609e-6, 43.9832e-6, 5.82852e-6}, syn1:pre_synapse(), syn1:post_synapse())
	synDistr:place_synapse_at_coords({-42.1714-6, 133.237e-6, 24.8023e-6}, syn1:pre_synapse(), syn1:post_synapse())
	
	syn2 = AlphaSynapsePair()
	syn2:set_id(1)
	syn2:set_onset(syn_onset[2])
	syn2:set_tau(syn_tau)                   -- default value
	syn2:set_gMax(syn_gMax)                 -- default value
	syn2:set_reversal_potential(syn_revPot) -- default value
	--synDistr:place_synapse_at_coords({-18.5196e-6, 40.3637e-6, 5.8769e-6}, syn2:pre_synapse(), syn2:post_synapse())
	synDistr:place_synapse_at_coords({-40.7199e-6, 130.633e-6, 23.3533e-6}, syn2:pre_synapse(), syn2:post_synapse())
	
	syn3 = AlphaSynapsePair()
	syn3:set_id(2)
	syn3:set_onset(syn_onset[3])
	syn3:set_tau(syn_tau)                   -- default value
	syn3:set_gMax(syn_gMax)                 -- default value
	syn3:set_reversal_potential(syn_revPot) -- default value
	--synDistr:place_synapse_at_coords({-18.6295e-6, 35.8769e-6, 5.87629e-6}, syn3:pre_synapse(), syn3:post_synapse())
	synDistr:place_synapse_at_coords({-40.6288e-6, 126.593e-6, 22.5147e-6}, syn3:pre_synapse(), syn3:post_synapse())
	
	syn4 = AlphaSynapsePair()
	syn4:set_id(3)
	syn4:set_onset(syn_onset[4])
	syn4:set_tau(syn_tau)                   -- default value
	syn4:set_gMax(syn_gMax)                 -- default value
	syn4:set_reversal_potential(syn_revPot) -- default value
	--synDistr:place_synapse_at_coords({-19.53e-6, 32.1163e-6, 5.87559e-6}, syn4:pre_synapse(), syn4:post_synapse())
	synDistr:place_synapse_at_coords({-42.0221e-6, 123.26e-6, 21.9985e-6}, syn4:pre_synapse(), syn4:post_synapse())
	
	syn5 = AlphaSynapsePair()
	syn5:set_id(4)
	syn5:set_onset(syn_onset[5])
	syn5:set_tau(syn_tau)                   -- default value
	syn5:set_gMax(syn_gMax)                 -- default value
	syn5:set_reversal_potential(syn_revPot) -- default value
	--synDistr:place_synapse_at_coords({-20.7218e-6, 28.1453e-6, 5.93687e-6}, syn5:pre_synapse(), syn5:post_synapse())
	synDistr:place_synapse_at_coords({-42.8684e-6, 119.663e-6, 21.1337e-6}, syn5:pre_synapse(), syn5:post_synapse())
	
	exportSuccess = synDistr:export_grid(gridSyn)
	if not exportSuccess then
	    print("Synapse distributor failed to export its grid to '" .. gridSyn .. "'.")
	    exit()
	end
end

gridName1d = gridSyn


--------------------------
-- biological settings	--
--------------------------
-- settings are according to T. Branco

-- membrane conductances (in units of S/m^2)
g_k_ax = 400.0	-- axon
g_k_so = 200.0	-- soma
g_k_de = 30		-- dendrite

g_na_ax = 3.0e4
g_na_so = 1.5e3
g_na_de = 40.0

g_l_ax = 200.0
g_l_so = 1.0
g_l_de = 1.0

-- specific capacitance (in units of F/m^2)
spec_cap = 1.0e-2

-- resistivity (in units of Ohm m)
spec_res = 1.5

-- reversal potentials (in units of V)
e_k  = -0.09
e_na = 0.06
e_ca = 0.14

-- equilibrium concentrations (in units of mM)
-- comment: these concentrations will not yield Nernst potentials
-- as given above; pumps will have to be introduced to achieve this
-- in the case where Nernst potentials are calculated from concentrations!
k_out  = 4.0
na_out = 150.0
ca_out = 1.5

k_in   = 140.0
na_in  = 10.0
ca_in  = 5e-5

-- equilibrium potential (in units of V)
v_eq = -0.065

-- diffusion coefficients (in units of m^2/s)
diff_k 	= 1.0e-9
diff_na	= 1.0e-9
diff_ca	= 2.2e-10

-- temperature in units of deg Celsius
temp = 37.0


------------------------------------
-- create 1d domain and approx space --
------------------------------------
neededSubsets1d = {"soma", "dendrite"} --, "axon"}
dom1d = util.CreateDomain(gridName1d, 0, neededSubsets1d)

approxSpace1d = ApproximationSpace(dom1d)
approxSpace1d:add_fct("v", "Lagrange", 1)
if withIons == true then
	approxSpace1d:add_fct("k", "Lagrange", 1)
	approxSpace1d:add_fct("na", "Lagrange", 1)
	approxSpace1d:add_fct("ca", "Lagrange", 1)
end

approxSpace1d:init_levels();
approxSpace1d:init_surfaces();
approxSpace1d:init_top_surface();
approxSpace1d:print_layout_statistic()
approxSpace1d:print_statistic()

OrderCuthillMcKee(approxSpace1d, true);

---------------------------
-- create 1d discretization --
---------------------------
allSubsets = "soma, dendrite" --, axon"

-- cable equation
CE = CableEquation(allSubsets, withIons)
CE:set_spec_cap(spec_cap)
CE:set_spec_res(spec_res)
CE:set_rev_pot_k(e_k)
CE:set_rev_pot_na(e_na)
CE:set_rev_pot_ca(e_ca)
CE:set_k_out(k_out)
CE:set_na_out(na_out)
CE:set_ca_out(ca_out)
CE:set_diff_coeffs({diff_k, diff_na, diff_ca})
CE:set_temperature_celsius(temp)

-- Hodgkin and Huxley channels
--HH = ChannelHHNernst("v, k, na", "axon")
if withIons == true then
	HH = ChannelHHNernst("v", allSubsets)
else
	HH = ChannelHH("v", allSubsets)
end
--HH:set_conductances(g_k_ax, g_na_ax, "axon")
HH:set_conductances(g_k_so, g_na_so, "soma")
HH:set_conductances(g_k_de, g_na_de, "dendrite")

CE:add(HH)


-- leakage
tmp_fct = math.pow(2.3,(temp-23.0)/10.0)

leak = ChannelLeak("v", allSubsets)
--leak:set_cond(g_l_ax*tmp_fct, "axon")
--leak:set_rev_pot(-0.066148458, "axon")
leak:set_cond(g_l_so*tmp_fct, "soma")
leak:set_rev_pot(-0.030654022, "soma")
leak:set_cond(g_l_de*tmp_fct, "dendrite")
leak:set_rev_pot(-0.057803624, "dendrite")

CE:add(leak)


-- synapses
syn_handler = SynapseHandler()
syn_handler:set_ce_object(CE)
CE:set_synapse_handler(syn_handler)


-- domain discretization
domDisc1d = DomainDiscretization(approxSpace1d)
domDisc1d:add(CE)




----------------------------------
-- constants for the 3d problem --
----------------------------------
-- total cytosolic calbindin concentration
-- (four times the real value in order to simulate four binding sites in one)
totalClb = 4*40.0e-6

-- diffusion coefficients
D_cac = 220.0
D_cae = 220.0
D_ip3 = 280.0
D_clb = 20.0

-- calbindin binding rates
k_bind_clb = 	27.0e06
k_unbind_clb = 	19

-- initial concentrations
ca_cyt_init = 5.0e-08
ca_er_init = 2.5e-4
ip3_init = 4.0e-8
clb_init = totalClb / (k_bind_clb/k_unbind_clb*ca_cyt_init + 1)


-- IP3 constants
reactionRateIP3 = 0.11
equilibriumIP3 = 4.0e-08
reactionTermIP3 = -reactionRateIP3 * equilibriumIP3

-- ER densities
IP3Rdensity = 17.3
RYRdensity = 0.86
leakERconstant = 3.8e-17
local v_s = 6.5e-27  -- V_S param of SERCA pump
local k_s = 1.8e-7   -- K_S param of SERCA pump
SERCAfluxDensity =   IP3Rdensity * 3.7606194166520605e-23        -- j_ip3r
			       + RYRdensity * 1.1204582669024472e-21       -- j_ryr
			       + leakERconstant * (ca_er_init-ca_cyt_init) -- j_leak
SERCAdensity = SERCAfluxDensity / (v_s/(k_s/ca_cyt_init+1.0)/ca_er_init)
if (SERCAdensity < 0) then error("SERCA flux density is outward for these density settings!") end

-- PM densities
pmcaDensity = 500.0
ncxDensity  = 15.0
vdccDensity = 1.0
leakPMconstant =  pmcaDensity * 6.9672131147540994e-24	-- single pump PMCA flux (mol/s)
				+ ncxDensity *  6.7567567567567566e-23	-- single pump NCX flux (mol/s)
				+ vdccDensity * (-1.5752042094823713e-25)    -- single channel VGCC flux (mol/s)
				-- *1.5 // * 0.5 for L-type // T-type
if (leakPMconstant < 0) then error("PM leak flux is outward for these density settings!") end


----------------------------------
-- setup 3d approximation space --
----------------------------------
-- create, load, refine and distribute domain
reqSubsets = {"cyt", "er", "pm", "erm"}--, "syn1", "syn2", "syn3", "syn4", "syn5"}
dom3d = util.CreateDomain(gridName3d, 0, reqSubsets)
balancer.partitioner = "parmetis"

-- protect ER membrane from being cut by partitioning
ccw = SubsetCommunicationWeights(dom3d)
ccw:set_weight_on_subset(1000.0, 3) -- mem_er

balancer.communicationWeights = ccw
balancer.staticProcHierarchy = true
balancer.firstDistLvl = -1
balancer.redistSteps = 0

balancer.ParseParameters()
balancer.PrintParameters()

-- in parallel environments: use a load balancer to distribute the grid
-- actual refinement and load balancing after setup of disc.
loadBalancer = balancer.CreateLoadBalancer(dom3d)

-- refining and distributing
-- manual refinement (need to update interface node location in each step)
if loadBalancer ~= nil then
	balancer.Rebalance(dom3d, loadBalancer)
end

if numRefs > 0 then	
	refiner = GlobalDomainRefiner(dom3d)	
	for i = 1, numRefs do
		refiner:refine()
	end
end

if loadBalancer ~= nil then
	print("Edge cut on base level: "..balancer.defaultPartitioner:edge_cut_on_lvl(0))
	loadBalancer:estimate_distribution_quality()
	loadBalancer:print_quality_records()
end
print(dom3d:domain_info():to_string())


--[[
--print("Saving domain grid and hierarchy.")
SaveDomain(dom3d, "refined_grid_p" .. ProcRank() .. ".ugx")
SaveGridHierarchyTransformed(dom3d:grid(), "refined_grid_hierarchy_p" .. ProcRank() .. ".ugx", 2.0)
print("Saving parallel grid layout")
SaveParallelGridLayout(dom3d:grid(), "parallel_grid_layout_p"..ProcRank()..".ugx", 2.0)
--]]

-- create approximation space
approxSpace3d = ApproximationSpace(dom3d)

cytVol = "cyt"
erVol = "er"
plMem = "pm"--, syn1, syn2, syn3, syn4, syn5"
plMem_vec = {"pm"}--, "syn1", "syn2", "syn3", "syn4", "syn5"}
erMem = "erm"

outerDomain = cytVol .. ", " .. plMem .. ", " .. erMem
innerDomain = erVol .. ", " .. erMem 

approxSpace3d:add_fct("ca_cyt", "Lagrange", 1, outerDomain)
approxSpace3d:add_fct("ca_er", "Lagrange", 1, innerDomain)
approxSpace3d:add_fct("clb", "Lagrange", 1, outerDomain)
approxSpace3d:add_fct("ip3", "Lagrange", 1, outerDomain)

approxSpace3d:init_levels();
approxSpace3d:init_surfaces();
approxSpace3d:init_top_surface();
approxSpace3d:print_layout_statistic()
approxSpace3d:print_statistic()

--------------------------
-- setup discretization --
--------------------------
-- diffusion --
diffCaCyt = ConvectionDiffusion("ca_cyt", cytVol, "fv1")
diffCaCyt:set_diffusion(D_cac)

diffCaER = ConvectionDiffusion("ca_er", erVol, "fv1")
diffCaER:set_diffusion(D_cae)

diffClb = ConvectionDiffusion("clb", cytVol, "fv1")
diffClb:set_diffusion(D_clb)

diffIP3 = ConvectionDiffusion("ip3", cytVol, "fv1")
diffIP3:set_diffusion(D_ip3)
diffIP3:set_reaction_rate(reactionRateIP3)
diffIP3:set_reaction(reactionTermIP3)


-- buffering --
discBuffer = BufferFV1(cytVol) -- where buffering occurs
discBuffer:add_reaction(
	"clb",						    -- the buffering substance
	"ca_cyt",						-- the buffered substance
	totalClb,						-- total amount of buffer
	k_bind_clb,					    -- binding rate constant
	k_unbind_clb)				    -- unbinding rate constant


-- er membrane transport systems
ip3r = IP3R({"ca_cyt", "ca_er", "ip3"})
ip3r:set_scale_inputs({1e3,1e3,1e3})
ip3r:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

ryr = RyR({"ca_cyt", "ca_er"})
--ryr = RyRinstat({"ca_cyt", "ca_er"}, erMemVec, approxSpace)
ryr:set_scale_inputs({1e3,1e3})
ryr:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

serca = SERCA({"ca_cyt", "ca_er"})
serca:set_scale_inputs({1e3,1e3})
serca:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

leakER = Leak({"ca_er", "ca_cyt"})
leakER:set_scale_inputs({1e3,1e3})
leakER:set_scale_fluxes({1e3}) -- from mol/(m^2 s) to (mol um)/(dm^3 s)


discIP3R = MembraneTransportFV1(erMem, ip3r)
discIP3R:set_density_function(IP3Rdensity)

discRyR = MembraneTransportFV1(erMem, ryr)
discRyR:set_density_function(RYRdensity)

discSERCA = MembraneTransportFV1(erMem, serca)
discSERCA:set_density_function(SERCAdensity)

discERLeak = MembraneTransportFV1(erMem, leakER)
discERLeak:set_density_function(1e12*leakERconstant/(1e3)) -- from mol/(um^2 s M) to m/s


-- plasma membrane transport systems
pmca = PMCA({"ca_cyt", ""})
pmca:set_constant(1, 1.0)
pmca:set_scale_inputs({1e3,1.0})
pmca:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

ncx = NCX({"ca_cyt", ""})
ncx:set_constant(1, 1.0)
ncx:set_scale_inputs({1e3,1.0})
ncx:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

leakPM = Leak({"", "ca_cyt"})
leakPM:set_constant(0, 1.0)
leakPM:set_scale_inputs({1.0,1e3})
leakPM:set_scale_fluxes({1e3}) -- from mol/(m^2 s) to (mol um)/(dm^3 s)

vdcc = VDCC_BG_CN({"ca_cyt", ""}, plMem_vec, approxSpace1d, approxSpace3d, "v")
vdcc:set_domain_disc_1d(domDisc1d)
vdcc:set_cable_disc(CE)
vdcc:set_coordinate_scale_factor_3d_to_1d(1e-6)
if withIons then
	vdcc:set_initial_values({v_eq, k_in, na_in, ca_in})
else
	vdcc:set_initial_values({v_eq})
end
vdcc:set_time_steps_for_simulation_and_potential_update(dt1d, dt1d)
vdcc:set_solver_output_verbose(verbose1d)
vdcc:set_vtk_output(filename.."vtk/solution1d", pstep)
vdcc:set_constant(1, 1.5)
vdcc:set_scale_inputs({1e3,1.0})
vdcc:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)
vdcc:set_channel_type_L() --default, but to be sure
vdcc:init(0.0)

discPMCA = MembraneTransportFV1(plMem, pmca)
discPMCA:set_density_function(pmcaDensity)

discNCX = MembraneTransportFV1(plMem, ncx)
discNCX:set_density_function(ncxDensity)

discPMLeak = MembraneTransportFV1(plMem, leakPM)
discPMLeak:set_density_function(1e12*leakPMconstant / (1.0-1e3*ca_cyt_init))

discVDCC = MembraneTransportFV1(plMem, vdcc)
discVDCC:set_density_function(vdccDensity)


synapseInflux = HybridSynapseCurrentAssembler(approxSpace3d, approxSpace1d, syn_handler, {"pm"}, "ca_cyt", "ip3")
synapseInflux:set_current_percentage(0.01)
synapseInflux:set_3d_neuron_ids({0})
synapseInflux:set_scaling_factors(1e-15, 1e-6, 1.0, 1e-15)
synapseInflux:set_valency(2)
synapseInflux:set_ip3_production_params(6e-20, 1.188)



-- domain discretization --
domDisc3d = DomainDiscretization(approxSpace3d)

domDisc3d:add(diffCaCyt)
domDisc3d:add(diffCaER)
domDisc3d:add(diffClb)
domDisc3d:add(diffIP3)

domDisc3d:add(discBuffer)

domDisc3d:add(discIP3R)
domDisc3d:add(discRyR)
domDisc3d:add(discSERCA)
domDisc3d:add(discERLeak)

domDisc3d:add(discPMCA)
domDisc3d:add(discNCX)
domDisc3d:add(discPMLeak)
domDisc3d:add(discVDCC)

domDisc3d:add(synapseInflux)

-- setup time discretization --
timeDisc = ThetaTimeStep(domDisc3d)
timeDisc:set_theta(1.0) -- 1.0 is implicit Euler

-- create operator from discretization
op = AssembledOperator()
op:set_discretization(timeDisc)
op:init()


------------------
-- solver setup --
------------------
-- debug writer
dbgWriter = GridFunctionDebugWriter(approxSpace3d)
dbgWriter:set_base_dir(filename)
dbgWriter:set_vtk_output(false)

-- biCGstab --
convCheck = ConvCheck()
convCheck:set_minimum_defect(1e-50)
convCheck:set_reduction(1e-8)
convCheck:set_verbose(verbose3d)

if (solverID == "ILU") then
    bcgs_steps = 10000
    ilu = ILU()
    ilu:set_sort(true)
    bcgs_precond = ilu
elseif (solverID == "GS") then
    bcgs_steps = 10000
    bcgs_precond = GaussSeidel()
else -- (solverID == "GMG")
	gmg = GeometricMultiGrid(approxSpace3d)
	gmg:set_discretization(timeDisc)
	gmg:set_base_level(0)
	gmg:set_gathered_base_solver_if_ambiguous(true)
	
	-- treat SuperLU problems with Dirichlet constraints by using constrained version
	gmg:set_base_solver(SuperLU())
	
	ilu_gmg = ILU()
	ilu_gmg:set_sort(true)		-- <-- SUPER-important!
	gmg:set_smoother(ilu_gmg)
	gmg:set_smooth_on_surface_rim(true)
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(3)
	gmg:set_num_postsmooth(3)
	--gmg:set_rap(true) -- causes error in base solver!!
	--gmg:set_debug(GridFunctionDebugWriter(approxSpace))
	
    bcgs_steps = 1000
	bcgs_precond = gmg
end

convCheck:set_maximum_steps(bcgs_steps)

bicgstabSolver = BiCGStab()
bicgstabSolver:set_preconditioner(bcgs_precond)
bicgstabSolver:set_convergence_check(convCheck)
--bicgstabSolver:set_debug(dbgWriter)

--- non-linear solver ---
-- convergence check
newtonConvCheck = CompositeConvCheck(approxSpace3d, 10, 1e-18, 1e-12)
--newtonConvCheck:set_component_check("ca_cyt, ca_er, clb, ip3", 1e-18, 1e-12)
newtonConvCheck:set_verbose(true)
newtonConvCheck:set_time_measurement(true)
--newtonConvCheck:set_adaptive(true)

-- Newton solver
newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(bicgstabSolver)
newtonSolver:set_convergence_check(newtonConvCheck)
--newtonSolver:set_debug(dbgWriter)

newtonSolver:init(op)


-------------
-- solving --
-------------
-- get grid function
u = GridFunction(approxSpace3d)

-- set initial value
InterpolateInner(ca_cyt_init, u, "ca_cyt", 0.0)
InterpolateInner(ca_er_init, u, "ca_er", 0.0)
InterpolateInner(clb_init, u, "clb", 0.0)
InterpolateInner(ip3_init, u, "ip3", 0.0)

-- timestep in seconds
dt = dt3dStart
time = 0.0
step = 0

-- initial vtk output
if (generateVTKoutput) then
	out = VTKOutput()
	out:print(filename .. "vtk/solution3d", u, step, time)
end

-- create new grid function for old value
uOld = u:clone()

-- store grid function in vector of old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)

min_dt = dt3d / math.pow(2,15)
cb_interval = 4
lv = startLv
levelUpDelay = math.max(6*syn_tau, 5e-3)
cb_counter = {}
for i=0,startLv do cb_counter[i]=0 end
while endTime-time > 0.001*dt do
	print("++++++ POINT IN TIME  " .. math.floor((time+dt)/dt+0.5)*dt .. "s  BEGIN ++++++")
	
	-- setup time Disc for old solutions and timestep
	timeDisc:prepare_step(solTimeSeries, dt)
	
	-- apply newton solver
	if newtonSolver:apply(u) == false
	then
		-- in case of failure:
		print ("Newton solver failed at point in time " .. time .. " with time step " .. dt)
		
		dt = dt/2
		lv = lv + 1
		VecScaleAssign(u, 1.0, solTimeSeries:latest())
		
		-- halve time step and try again unless time step below minimum
		if dt < min_dt
		then 
			print ("Time step below minimum. Aborting. Failed at point in time " .. time .. ".")
			time = endTime
		else
			print ("Trying with half the time step...")
			cb_counter[lv] = 0
		end
	else
		-- update new time
		time = solTimeSeries:time(0) + dt
		
		-- update check-back counter and, if applicable, reset dt
		cb_counter[lv] = cb_counter[lv] + 1
		while cb_counter[lv] % (2*cb_interval) == 0 and lv > 0 and (time >= levelUpDelay or lv > startLv) do
			print ("Doubling time due to continuing convergence; now: " .. 2*dt)
			dt = 2*dt;
			lv = lv - 1
			cb_counter[lv] = cb_counter[lv] + cb_counter[lv+1] / 2
			cb_counter[lv+1] = 0
		end
		
		-- plot solution every pstep seconds
		if (generateVTKoutput) then
			if math.abs(time/pstep - math.floor(time/pstep+0.5)) < 1e-5 then
				out:print(filename .. "vtk/solution3d", u, math.floor(time/pstep+0.5), time)
			end
		end
		
		-- get oldest solution
		oldestSol = solTimeSeries:oldest()
		
		-- copy values into oldest solution (we reuse the memory here)
		VecScaleAssign(oldestSol, 1.0, u)
		
		-- push oldest solutions with new values to front, oldest sol pointer is popped from end
		solTimeSeries:push_discard_oldest(oldestSol, time)
		
		print("++++++ POINT IN TIME  " .. math.floor(time/dt+0.5)*dt .. "s  END ++++++++");
	end

end

-- end timeseries, produce gathering file
if (generateVTKoutput) then out:write_time_pvd(filename .. "vtk/solution3d", u) end

