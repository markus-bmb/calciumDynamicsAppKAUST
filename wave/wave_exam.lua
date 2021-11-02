--------------------------------------------------------------------------------
-- Examination of calcium wave prerequisites                                  --
--                                                                            --
-- This script is intended to be used for simulations on 2d representations   --
-- of perfectly rotationally symmetric model dendrites.                       --
-- It is supposed to be called by the script 'wave_exam_batch.lua' during a   --
-- bisection to find thresholds for ER radius and RyR channel density above   --
-- which a calcium wave is elicited.                                          --
-- The LIMEX method is used for time stepping.                                --
--                                                                            --
-- This script has been used to achieve the results published in:             --
-- Breit & Queisser: "What is required for neuronal calcium waves?            --
-- A numerical parameter study.", Journal of Mathematical Neuroscience (2018) --
--                                                                            --
-- Author: Markus Breit                                                       --
-- Date:   2017-08-29                                                         --
--------------------------------------------------------------------------------

-- for profiler output
--SetOutputProfileStats(true)

-- load pre-implemented lua functions
ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")
ug_load_script("plugins/Limex/limex_util.lua")

AssertPluginsLoaded({"neuro_collection", "Limex", "Parmetis"})

EnableLUA2C(true)  -- speed up evaluation of lua functions by c program
--SetDebugLevel(debugID.LUACompiler, 0) 


-------------------------------------
-- parse command line parameters  ---
-------------------------------------

-- choice of grid name
gridName = util.GetParam("-grid", "modelDendrite.ugx")

-- grid parameters
dendLength = util.GetParamNumber("-dendLength", 50.0)
dendRadius = util.GetParamNumber("-dendRadius", dendRadius or 0.5)
erRadius = util.GetParamNumber("-erRadius", erRadius or 0.158)
nSeg = util.GetParamNumber("-nSeg", 96)

-- refinements (global and at ERM)
numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numGlobRefs = util.GetParamNumber("-numGlobRefs", 0)
numERMRefs = util.GetParamNumber("-numERMRefs", 0)

-- which ER mechanisms are to be activated?
setting = util.GetParam("-setting", "ryr")
setting = string.lower(setting)
validSettings = {}
validSettings["all"] = 0;
validSettings["none"] = 0;
validSettings["ip3r"] = 0;
validSettings["ryr"] = 0;
if (validSettings[setting] == nil) then
    error("Unknown setting " .. setting)
end

-- densities
ryrDens = util.GetParamNumber("-ryrDens", ryrDens or 0.86)

-- buffer
totalBuffer = util.GetParamNumber("-totBuf", 4*40.0e-6)

-- whether to scale synaptic influx with dendritic radius
scaledInflux = util.HasParamOption("-scaledInflux")

-- choice of algebra
useBlockAlgebra = util.HasParamOption("-block")

-- choice of solver setup
solverID = util.GetParam("-solver", "GMG")
solverID = string.upper(solverID)
validSolverIDs = {}
validSolverIDs["GMG"] = 0
validSolverIDs["GS"] = 0
validSolverIDs["ILU"] = 0
if (validSolverIDs[solverID] == nil) then
    error("Unknown solver identifier " .. solverID)
end

-- error tolerance for Limex iteration
toleratedError = util.GetParamNumber("-tol", 0.01)
nstages = util.GetParamNumber("-nst", 3)

-- specify "-verbose" to output linear solver convergence
verbose = util.HasParamOption("-verbose")

-- parameters for instationary simulation
dt = util.GetParamNumber("-dt", 1e-2)
endTime = util.GetParamNumber("-endTime", 1.0)

-- choose outfile directory
outDir = util.GetParam("-outName", "caWaveExploration")
outDir = outDir .. "/"

-- specify -vtk to generate vtk output
generateVTKoutput = util.HasParamOption("-vtk")
pstep = util.GetParamNumber("-pstep", dt, "plotting interval")



-- init with dimension and algebra
if useBlockAlgebra then
	InitUG(2, AlgebraType("CPU", 7))
else
	InitUG(2, AlgebraType("CPU", 1))
end

-----------------------
-- geometry creation --
-----------------------
if ProcRank() == 0 then
	gen = DendriteGenerator()
	gen:set_dendrite_length(dendLength)
	gen:set_dendrite_radius(dendRadius)
	gen:set_er_radius(erRadius)
	gen:set_num_segments(nSeg)
	
	gridName = outDir .. "grid/" .. gridName
	gen:create_dendrite(gridName)
end

PclDebugBarrierAll()

-------------------------
--  problem constants  --
-------------------------
-- setting-dependent variables
withIP3R = true
withRyR = true
withSERCAandLeak = true

if setting == "none" then 
	withIP3R = false
	withRyR = false
	withSERCAandLeak = false
end

if setting == "ip3r" then
	withRyR = false
end

if setting == "ryr" then
	withIP3R = false
end


-- total cytosolic calbindin concentration
-- (four times the real value in order to simulate four binding sites in one)
totalClb = totalBuffer

-- diffusion coefficients
D_cac = 220.0
D_cae = 220.0
D_ip3 = 280.0
D_clb = 20.0

-- calbindin binding rates
k_bind_clb = 27.0e06
k_unbind_clb = 19

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
RYRdensity = ryrDens --0.86
local v_s = 6.5e-27  -- V_S param of SERCA pump
local k_s = 1.8e-7   -- K_S param of SERCA pump
local j_ip3r = 3.7606194166520605e-23   -- single channel IP3R flux (mol/s) - to be determined via gdb
local j_ryr = 1.1201015633466695e-21    -- single channel RyR flux (mol/s) - to be determined via gdb
				  						-- ryr1: 1.1204582669024472e-21	
---[[
-- equilibration using SERCA
leakERconstant = 3.8e-17
local j_leak = ca_er_init-ca_cyt_init	-- leak proportionality factor


SERCAfluxDensity = leakERconstant * j_leak
if withIP3R then 
	SERCAfluxDensity = SERCAfluxDensity + IP3Rdensity * j_ip3r
end
if withRyR then
	SERCAfluxDensity = SERCAfluxDensity + RYRdensity * j_ryr
end
SERCAdensity = SERCAfluxDensity / (v_s/(k_s/ca_cyt_init+1.0)/ca_er_init)
if (SERCAdensity < 0) then error("SERCA flux density is outward for these density settings!") end
--]]
--[[
-- equilibration using leakage
SERCAdensity = 1973.0
SERCAflux = v_s / (k_s / ca_cyt_init + 1.0) / ca_er_init

netEquilFlux = SERCAdensity*SERCAflux
if withIP3R then 
	netEquilFlux = netEquilFlux - IP3Rdensity * j_ip3r
end
if withRyR then
	netEquilFlux = netEquilFlux - RYRdensity * j_ryr
end

leakERconstant = netEquilFlux / (ca_er_init - ca_cyt_init)
if (leakERconstant < 0) then
	error("ER leakage flux density is outward for these density settings!")
end
--]]

-- PM densities
pmcaDensity = 500.0
ncxDensity  = 15.0
vdccDensity = 0.0  -- 1.0
leakPMconstant =  pmcaDensity * 6.9672131147540994e-24	-- single pump PMCA flux (mol/s)
				+ ncxDensity *  6.7567567567567566e-23	-- single pump NCX flux (mol/s)
				--+ vdccDensity * (-1.5752042094823713e-25)    -- single channel VGCC flux (mol/s)
				-- *1.5 // * 0.5 for L-type // T-type
if (leakPMconstant < 0) then error("PM leak flux is outward for these density settings!") end


-- activation pattern
caEntryDuration = 0.001
function synCurrentDensityCa(z, r, t, si)
	-- single spike (~1200 ions)
	local influx = 0.0
	if t <= caEntryDuration	then
		influx = 2.5e-3 * (1.0 - t/caEntryDuration)
	end
	
    return 2.0*math.pi*r * influx
end

ip3EntryDelay = 0.000
ip3EntryDuration = 0.2
function synCurrentDensityIP3(z, r, t, si)
	local influx = 0.0
	if t > ip3EntryDelay and t <= ip3EntryDelay+ip3EntryDuration then
		influx = 7.5e-5 * (1.0 - t/ip3EntryDuration)
	end
	
    return 2.0*math.pi*r * influx
end

-------------------------------
-- setup approximation space --
-------------------------------

-- load domain
reqSubsets = {"cyt", "er", "pm", "erm", "act", "meas"}
dom = util.CreateDomain(gridName, numPreRefs, reqSubsets)

-- create approximation space
approxSpace = ApproximationSpace(dom)

cytVol = "cyt"
erVol = "er"
plMem = "pm"
plMem_vec = {"pm"}
erMem = "erm"
erMemVec = {"erm"}

outerDomain = cytVol .. ", " .. plMem .. ", " .. erMem .. ", act, meas"
innerDomain = erVol .. ", " .. erMem

if useBlockAlgebra then
	approxSpace:add_fct("ca_cyt", "Lagrange", 1)
	approxSpace:add_fct("ca_er", "Lagrange", 1)
	approxSpace:add_fct("clb", "Lagrange", 1)
	if withIP3R then
		approxSpace:add_fct("ip3", "Lagrange", 1)
	end
	if withRyR then
		approxSpace:add_fct("o2", "Lagrange", 1)
		approxSpace:add_fct("c1", "Lagrange", 1)
		approxSpace:add_fct("c2", "Lagrange", 1)
	end
else
	approxSpace:add_fct("ca_cyt", "Lagrange", 1, outerDomain)
	approxSpace:add_fct("ca_er", "Lagrange", 1, innerDomain)
	approxSpace:add_fct("clb", "Lagrange", 1, outerDomain)
	if withIP3R then
		approxSpace:add_fct("ip3", "Lagrange", 1, outerDomain)
	end
	if withRyR then
		approxSpace:add_fct("o2", "Lagrange", 1, erMem)
		approxSpace:add_fct("c1", "Lagrange", 1, erMem)
		approxSpace:add_fct("c2", "Lagrange", 1, erMem)
	end
end
approxSpace:init_levels()
approxSpace:init_surfaces()
approxSpace:init_top_surface()
approxSpace:print_layout_statistic()
approxSpace:print_statistic()

if useBlockAlgebra then
	OrderCuthillMcKee(approxSpace, true)
end


-- in parallel environments: domain distribution
balancer.partitioner = "parmetis"
balancer.staticProcHierarchy = true
balancer.firstDistLvl = numPreRefs
balancer.firstDistProcs = 96
balancer.redistSteps = 4
balancer.parallelElementThreshold = 4

balancer.ParseParameters()
balancer.PrintParameters()

loadBalancer = balancer.CreateLoadBalancer(dom)
if loadBalancer ~= nil then
	loadBalancer:enable_vertical_interface_creation(solverID == "GMG")
	if balancer.partitioner == "parmetis" then
		au = AnisotropyUnificator(dom)
		au:set_threshold_ratio(0.1)
		mu = ManifoldUnificator(dom)
		mu:add_protectable_subsets("erm")
		cdgm = ClusteredDualGraphManager()
		cdgm:add_unificator(SiblingUnificator())
		cdgm:add_unificator(au)
		cdgm:add_unificator(mu)
		balancer.defaultPartitioner:set_dual_graph_manager(cdgm)
	end
	balancer.Rebalance(dom, loadBalancer)
	loadBalancer:estimate_distribution_quality()
	loadBalancer:print_quality_records()
	if balancer.partitioner == "parmetis" then
		print("Edge cut on base level: "..balancer.defaultPartitioner:edge_cut_on_lvl(0))
	end
end


-- ERM refinements
refiner = HangingNodeDomainRefiner(dom)

for i = 1, numERMRefs do
	MarkAlongSurface(refiner, dom, {"erm"}, {"cyt"})
	refiner:refine()
end

-- global refinements
for i = 1, numGlobRefs do
	MarkGlobal(refiner, dom)
	refiner:refine()
end

print(dom:domain_info():to_string())
--SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), outDir .. "grid/refined_grid_hierarchy_p" .. ProcRank() .. ".ugx", 1.0)
--SaveParallelGridLayout(dom:grid(), outDir .. "grid/parallel_grid_layout_p"..ProcRank()..".ugx", 1.0)


--------------------------
-- setup discretization --
--------------------------
-- scaled coefficient functions
function scaled_diff_cac(z,r,t)
	local c = 2*math.pi*r*D_cac
	return  c, 0,
			0, c
end
function scaled_diff_cae(z,r,t)
	local c = 2*math.pi*r*D_cae
	return  c, 0,
			0, c
end
function scaled_diff_clb(z,r,t)
	local c = 2*math.pi*r*D_clb
	return  c, 0,
			0, c
end
function scaled_diff_ip3(z,r,t)
	local c = 2*math.pi*r*D_ip3
	return  c, 0,
			0, c
end
function scaled_reactionRate_ip3(z,r,t)
	return 2*math.pi*r*reactionRateIP3
end
function scaled_reactionTerm_ip3(z,r,t)
	return 2*math.pi*r*reactionTermIP3
end
function buffer_kBind(z,r,t)
	return 2*math.pi*r*k_bind_clb
end
function buffer_kUnbind(z,r,t)
	return 2*math.pi*r*k_unbind_clb
end
function rotSym_scale(z,r,t)
	return 2*math.pi*r
end


-- diffusion --
diffCaCyt = ConvectionDiffusion("ca_cyt", cytVol, "fv1")
diffCaCyt:set_mass_scale("rotSym_scale")
diffCaCyt:set_diffusion("scaled_diff_cac")

diffCaER = ConvectionDiffusion("ca_er", erVol, "fv1")
diffCaER:set_mass_scale("rotSym_scale")
diffCaER:set_diffusion("scaled_diff_cae")

diffClb = ConvectionDiffusion("clb", cytVol, "fv1")
diffClb:set_mass_scale("rotSym_scale")
diffClb:set_diffusion("scaled_diff_clb")

if withIP3R then
	diffIP3 = ConvectionDiffusion("ip3", cytVol, "fv1")
	diffIP3:set_mass_scale("rotSym_scale")
	diffIP3:set_diffusion("scaled_diff_ip3")
	diffIP3:set_reaction_rate("scaled_reactionRate_ip3")
	diffIP3:set_reaction("scaled_reactionTerm_ip3")
end

-- buffering --
discBuffer = BufferFV1(cytVol) -- where buffering occurs
discBuffer:add_reaction(
	"clb",                     -- the buffering substance
	"ca_cyt",                  -- the buffered substance
	totalClb,                  -- total amount of buffer
	"buffer_kBind",            -- binding rate constant
	"buffer_kUnbind"           -- unbinding rate constant
)


-- er membrane transport systems
if withIP3R then
	ip3r = IP3R({"ca_cyt", "ca_er", "ip3"})
	ip3r:set_scale_inputs({1e3,1e3,1e3})
	ip3r:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

	discIP3R = MembraneTransportFV1(erMem, ip3r)
	discIP3R:set_density_function(IP3Rdensity)
	discIP3R:set_flux_scale("rotSym_scale")  -- to achieve 3d rot. symm. simulation in 2d
end

if withRyR then
	ryr = RyRImplicit({"ca_cyt", "ca_er", "o2", "c1", "c2"}, erMemVec)
	ryr:set_scale_inputs({1e3, 1e3, 1.0, 1.0, 1.0})
	ryr:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

	discRyR = MembraneTransportFV1(erMem, ryr)
	discRyR:set_density_function(RYRdensity)
	discRyR:set_flux_scale("rotSym_scale")  -- to achieve 3d rot. symm. simulation in 2d
end

if withSERCAandLeak then
	serca = SERCA({"ca_cyt", "ca_er"})
	serca:set_scale_inputs({1e3,1e3})
	serca:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)
	
	leakER = Leak({"ca_er", "ca_cyt"})
	leakER:set_scale_inputs({1e3,1e3})
	leakER:set_scale_fluxes({1e3}) -- from mol/(m^2 s) to (mol um)/(dm^3 s)

	discSERCA = MembraneTransportFV1(erMem, serca)
	discSERCA:set_density_function(SERCAdensity)
	discSERCA:set_flux_scale("rotSym_scale")  -- to achieve 3d rot. symm. simulation in 2d
	
	discERLeak = MembraneTransportFV1(erMem, leakER)
	discERLeak:set_density_function(1e12*leakERconstant/(1e3)) -- from mol/(um^2 s M) to m/s
	discERLeak:set_flux_scale("rotSym_scale")  -- to achieve 3d rot. symm. simulation in 2d
end

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

--[[
vdcc = VDCC_BG_CN({"ca_cyt", ""}, plMem_vec, approxSpace1d, approxSpace, "v")
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
if generateVTKoutput then
	vdcc:set_vtk_output(outDir.."vtk/solution1d", pstep)
end
vdcc:set_constant(1, 1.5)
vdcc:set_scale_inputs({1e3,1.0})
vdcc:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)
vdcc:set_channel_type_L() --default, but to be sure
vdcc:init(0.0)
--]]

discPMCA = MembraneTransportFV1(plMem, pmca)
discPMCA:set_density_function(pmcaDensity)
discPMCA:set_flux_scale("rotSym_scale")  -- to achieve 3d rot. symm. simulation in 2d

discNCX = MembraneTransportFV1(plMem, ncx)
discNCX:set_density_function(ncxDensity)
discNCX:set_flux_scale("rotSym_scale")  -- to achieve 3d rot. symm. simulation in 2d

discPMLeak = MembraneTransportFV1(plMem, leakPM)
discPMLeak:set_density_function(1e12*leakPMconstant / (1.0-1e3*ca_cyt_init))
discPMLeak:set_flux_scale("rotSym_scale")  -- to achieve 3d rot. symm. simulation in 2d

--[[
discVDCC = MembraneTransportFV1(plMem, vdcc)
discVDCC:set_density_function(vdccDensity)
discVDCC:set_flux_scale("rotSym_scale")  -- to achieve 3d rot. symm. simulation in 2d
--]]

-- synaptic activity
synapseInfluxCa = UserFluxBoundaryFV1("ca_cyt", "act")
synapseInfluxCa:set_flux_function("synCurrentDensityCa")
if withIP3R then
	synapseInfluxIP3 = UserFluxBoundaryFV1("ip3", "act")
	synapseInfluxIP3:set_flux_function("synCurrentDensityIP3")
end

-- domain discretization --
domDisc = DomainDiscretization(approxSpace)

domDisc:add(diffCaCyt)
domDisc:add(diffCaER)
domDisc:add(diffClb)

domDisc:add(discBuffer)

if withIP3R then
	domDisc:add(diffIP3)
	domDisc:add(discIP3R)
	domDisc:add(synapseInfluxIP3)
end
if withRyR then
	domDisc:add(discRyR)
	domDisc:add(ryr) -- also add ryr as elem disc (for state variables)
end
if withSERCAandLeak then
	domDisc:add(discSERCA)
	domDisc:add(discERLeak)
end

domDisc:add(discPMCA)
domDisc:add(discNCX)
domDisc:add(discPMLeak)

domDisc:add(synapseInfluxCa)


-- Dirichlet for superfluous dofs
if useBlockAlgebra then
	uselessDofDiri = DirichletBoundary()
	uselessDofDiri:add(ca_cyt_init, "ca_cyt", "er")
	if withIP3R then
		uselessDofDiri:add(ip3_init, "ip3", "er")
	end
	uselessDofDiri:add(clb_init, "clb", "er")
	uselessDofDiri:add(ca_er_init, "ca_er", "cyt, pm, act, meas")
	if withRyR then
		uselessDofDiri:add(0, "o2", "cyt, er, pm, act, meas")
		uselessDofDiri:add(0, "c1", "cyt, er, pm, act, meas")
		uselessDofDiri:add(0, "c2", "cyt, er, pm, act, meas")
	end
	domDisc:add(uselessDofDiri)
end

-- constraints for adaptivity
if numERMRefs > 0 then
	hangingConstraint = SymP1Constraints()--OneSideP1Constraints()
	domDisc:add(hangingConstraint)
end

-- setup time discretization --
timeDisc = ThetaTimeStep(domDisc)
timeDisc:set_theta(1.0) -- 1.0 is implicit Euler

-- create operator from discretization
op = AssembledOperator()
op:set_discretization(timeDisc)
op:init()


------------------
-- solver setup --
------------------
-- debug writer
dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_base_dir(outDir)
dbgWriter:set_vtk_output(false)

-- biCGstab --
convCheck = ConvCheck()
convCheck:set_minimum_defect(1e-50)
convCheck:set_reduction(1e-8)
convCheck:set_verbose(verbose)

if (solverID == "ILU") then
    bcgs_steps = 10000
    ilu = ILU()
    ilu:set_sort(true)
    bcgs_precond = ilu
elseif (solverID == "GS") then
    bcgs_steps = 1000
    bcgs_precond = GaussSeidel()
else -- (solverID == "GMG")
	gmg = GeometricMultiGrid(approxSpace)
	gmg:set_discretization(timeDisc)
	gmg:set_base_level(numPreRefs)
	gmg:set_gathered_base_solver_if_ambiguous(true)
	
	-- treat SuperLU problems with Dirichlet constraints by using constrained version
	gmg:set_base_solver(SuperLU())
	
	smoother = GaussSeidel()
	gmg:set_smoother(smoother)
	gmg:set_smooth_on_surface_rim(true)
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(3)
	gmg:set_num_postsmooth(3)
	--gmg:set_rap(true)
	--gmg:set_debug(GridFunctionDebugWriter(approxSpace))
	
    bcgs_steps = 100
	bcgs_precond = gmg
end

convCheck:set_maximum_steps(bcgs_steps)

bicgstabSolver = BiCGStab()
bicgstabSolver:set_preconditioner(bcgs_precond)
bicgstabSolver:set_convergence_check(convCheck)
--bicgstabSolver:set_debug(dbgWriter)

-- Newton solver
newtonSolver = LimexNewtonSolver()
newtonSolver:set_linear_solver(bicgstabSolver)
--newtonSolver:set_convergence_check(newtonConvCheck)
--newtonSolver:set_debug(dbgWriter)

newtonSolver:init(op)


-------------
-- solving --
-------------
-- get grid function
u = GridFunction(approxSpace)

-- set initial value
InterpolateInner(ca_cyt_init, u, "ca_cyt", 0.0)
InterpolateInner(ca_er_init, u, "ca_er", 0.0)
InterpolateInner(clb_init, u, "clb", 0.0)
if withIP3R then
	InterpolateInner(ip3_init, u, "ip3", 0.0)
end
if withRyR then
	ryr:calculate_steady_state(u)
end

-- timestep in seconds
dtmin = 1e-9
dtmax = 1e-2
time = 0.0
step = 0

-- initial vtk output
--if generateVTKoutput then
--	out = VTKOutput()
--	out:print(outDir .. "vtk/solution", u, step, time)
--end


------------------
--  LIMEX setup --
------------------
stageNSteps = {}    -- number of time steps for each stage
for i = 1, nstages do stageNSteps[i] = i end

limex = LimexTimeIntegrator(nstages)
for i = 1, nstages do
	limex:add_stage(stageNSteps[i], newtonSolver, domDisc)
end

limex:set_tolerance(toleratedError)
limex:set_time_step(dt)
limex:set_dt_min(dtmin)
limex:set_dt_max(dtmax)
limex:set_increase_factor(2.0)
limex:set_reduction_factor(0.1)
limex:set_stepsize_greedy_order_factor(1)
limex:set_stepsize_safety_factor(0.25)

-- GridFunction error estimator (relative norm)
limexEstimator = ScaledGridFunctionEstimator()

errorEvalCaCyt = H1ComponentSpace("ca_cyt", "cyt", 3)  -- fct names, subset names, order
errorEvalCaER = H1ComponentSpace("ca_er", "er", 3)
errorEvalClb = H1ComponentSpace("clb", "cyt", 3)
limexEstimator:add(errorEvalCaCyt)
limexEstimator:add(errorEvalCaER)
limexEstimator:add(errorEvalClb)

if withIP3R then
	errorEvalIP3 = H1ComponentSpace("ip3", "cyt", 3)
	limexEstimator:add(errorEvalIP3)
end
if withRyR then
	errorEvalO2 = L2ComponentSpace("o2", 3, 1.0, "erm")  -- fct names, order, weight, subset names
	errorEvalC1 = L2ComponentSpace("c1", 3, 1.0, "erm")
	errorEvalC2 = L2ComponentSpace("c2", 3, 1.0, "erm")
	limexEstimator:add(errorEvalO2)
	limexEstimator:add(errorEvalC1)
	limexEstimator:add(errorEvalC2)
end

limex:add_error_estimator(limexEstimator)


-- for vtk output
if generateVTKoutput then
	local vtkObserver = VTKOutputObserver(outDir .."vtk/solution", out, pstep)
	limex:attach_observer(vtkObserver)
end

-- for measurements
waveHitEnd = false
waveGotStuck = false
waveFrontXPos = -0.5*dendLength
maxRyRFluxDens = 0.0
stuckWaveXPos = 0.0
interruptTime = 0.0


-- prepare output
if ProcRank() == 0 then
	waveFrontPosFile = outDir.."meas/waveFrontX.dat"
	waveFrontPosFH = assert(io.open(waveFrontPosFile, "a"))
	waveFrontPosFH:write(0.0, "\t", waveFrontXPos + 0.5*dendLength, "\n")
	waveFrontPosFH:flush()
end
measInterval = 1e-4
lastMeasPt = 0
--wpe = WaveProfileExporter(approxSpace, "ca_cyt", "erm", outDir .. "meas/waveProfile")


function measWaveActivity(step, time, dt)
	curSol = measObserver:get_current_solution()

	-- measure free Ca in dendrite
	take_measurement(curSol, time, "cyt", "ca_cyt", outDir.."meas/meanCaCyt.dat")
	
	-- measure concentration at right end
	local measConc = take_measurement(curSol, time, "meas", "ca_cyt", outDir.."meas/caAtRightEnd_erRad"..erRadius.."_ryrDens"..ryrDens)
	if measConc > 4*ca_cyt_init then
		waveHitEnd = true
		interruptTime = time
		limex:interrupt()
	end
	
	-- measure wave front x position
	if withRyR then
		lastWaveFrontXPos = waveFrontXPos
		waveFrontXPos = wave_front_x(curSol, "c1, c2", "erm", 0.1)
		--waveFrontXPos = wave_front_x(curSol, "ca_cyt", "erm", 6e-8)
		if waveFrontXPos < -1e10 then
			waveFrontXPos = -0.5*dendLength
		end		

		local measPt = math.floor(time/measInterval)
		if measPt > lastMeasPt then
			-- write current wave front location to file
			if ProcRank() == 0 then
				waveFrontPosFH:write(time, "\t", waveFrontXPos + 0.5*dendLength, "\n")
				waveFrontPosFH:flush()
			end
			
			-- write complete wave profile to file
			--wpe:exportWaveProfileX(curSol, time)
			
			lastMeasPt = measPt
		end
		
		
		-- in case the wave front becomes stuck: terminate
		---[[
		local ryrFluxDens = max_ryr_flux_density(curSol, "ca_cyt, ca_er, c1, c2", "erm", ryr)
		if ryrFluxDens > maxRyRFluxDens then
			maxRyRFluxDens = ryrFluxDens
		elseif ryrFluxDens < 0.25 * maxRyRFluxDens then -- wave got stuck
			waveGotStuck = true
			stuckWaveXPos = waveFrontXPos
			interruptTime = time
			limex:interrupt()
		end
		--]]
		
		--[[
		if waveFrontXPos < lastWaveFrontXPos then 
			waveGotStuck = true
			stuckWaveXPos = lastWaveFrontXPos
			interruptTime = time
			limex:interrupt()
		end
		--]]
	end
			
	print("Current (real) time: " .. time .. ",   last dt: " .. dt)
		
	return 0.0
end

measObserver = LuaCallbackObserver()
measObserver:set_callback("measWaveActivity")
limex:attach_observer(measObserver)

	
-- solve problem
limex:apply(u, endTime, u, time)

-- output outcome
print("")
if waveHitEnd then
	avgVelocity = dendLength / interruptTime / 1000.0
	print("#######################################")
	print("## A calcium wave has been detected. ##")
	print(string.format("## Average velocity was %5.3f um/ms. ##", avgVelocity))
	print("#######################################")
elseif waveGotStuck then
	interruptTimeMs = interruptTime * 1000.0
	print("##################################")
	print("## A calcium wave terminated at ##")
	print(string.format("## t = %5.2f ms, x = %5.2f um.  ##", interruptTimeMs, stuckWaveXPos+0.5*dendLength))
	print("##################################")
else
	print("###########################################")
	print("## A calcium wave has been elicited, but ##")
	print("## neither terminated nor hit the end.   ##")
	print("###########################################")
end


-- close output files
if ProcRank() == 0 then
	waveFrontPosFH:close()
end


--if (generateVTKoutput) then
--	out:write_time_pvd(outDir .. "vtk/solution", u)
--end


if doProfiling then
	WriteProfileData(outDir .."pd.pdxml")
end
