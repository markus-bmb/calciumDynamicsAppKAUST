--------------------------------------------------------------------------------
-- Example script for calcium simulations on a 3d reconstructed spine using   --
-- LIMEX for time stepping and treatment of non-linearities.                  --
--                                                                            --
--  Author: Markus Breit                                                      --
--  Date: 2019-07-11                                                          --
--------------------------------------------------------------------------------

-- for profiler output
SetOutputProfileStats(false)

-- load pre-implemented lua functions
ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

-- choose dimension and algebra
InitUG(3, AlgebraType("CPU", 1))

-- choice of grid
gridName = util.GetParam("-grid", "calciumDynamics_app/grids/reconstructed_spine.ugx")

-- total refinements
numRefs = util.GetParamNumber("-numRefs", 0)

-- choose length of maximal time step during the whole simulation
dt = util.GetParamNumber("-dt", 1e-05)

-- choose end time
endTime = util.GetParamNumber("-endTime")
if (endTime == nil)
then
	-- choose number of time steps
	nTimeSteps = util.GetParamNumber("-nTimeSteps", 1)
	endTime = nTimeSteps*timeStep
end

-- choose solver setup
solverID = util.GetParam("-solver", "GS")
solverID = string.upper(solverID)
validSolverIDs = {}
validSolverIDs["GMG"] = 0
validSolverIDs["GS"] = 0
validSolverIDs["ILU"] = 0
if (validSolverIDs[solverID] == nil) then
    error("Unknown solver identifier " .. solverID)
end

-- error tolerance for Limex iteration
tol = util.GetParamNumber("-tol", 0.01)

-- specify -verbose to display linear solver output
verbose = util.HasParamOption("-verbose")
 
-- choose outfile directory
outDir = util.GetParam("-outName", "CD/reconstructed")
outDir = outDir.."/"

-- specify -vtk to generate vtk output
generateVTKoutput = util.HasParamOption("-vtk")

-- choose plotting interval
pstep = util.GetParamNumber("-pstep", dt, "plotting interval")


---------------
-- constants --
---------------
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
ca_cyt_init = 5.0e-08 --4.0e-8
ca_er_init = 2.5e-4
ca_ext = 1.0e-3
ip3_init = 4.0e-8
clb_init = totalClb / (k_bind_clb/k_unbind_clb*ca_cyt_init + 1)


-- reaction reate IP3
reactionRateIP3 = 0.11

-- equilibrium concentration IP3
equilibriumIP3 = 4.0e-08

-- reation term IP3
reactionTermIP3 = -reactionRateIP3 * equilibriumIP3

IP3Rdensity = 17.3
RYRdensity = 0.86
leakERconstant = 3.8e-17

-- this is a little bit more complicated, since it must be ensured that
-- the net flux for equilibrium concentrations is zero
-- MUST be adapted whenever any parameterization of ER flux mechanisms is changed!
local v_s = 6.5e-27						-- V_S param of SERCA pump
local k_s = 1.8e-7						-- K_S param of SERCA pump
local j_ip3r = 3.7606194166520605e-23   -- single channel IP3R flux (mol/s) - to be determined via gdb
local j_ryr = 1.1201015633466695e-21    -- single channel RyR flux (mol/s) - to be determined via gdb
				  							-- ryr1: 1.1204582669024472e-21	
local j_leak = ca_er_init-ca_cyt_init	-- leak proportionality factor
	
SERCAdensity = leakERconstant * j_leak
if withIP3R then 
	SERCAdensity = SERCAdensity + IP3Rdensity * j_ip3r
end
if withRyR then
	SERCAdensity = SERCAdensity + RYRdensity * j_ryr
end
SERCAdensity = SERCAdensity / (v_s/(k_s/ca_cyt_init+1.0)/ca_er_init)

pmcaDensity = 500.0
ncxDensity  = 15.0
vdccDensity = 1.0

leakPMconstant =  pmcaDensity * 6.9672131147540994e-24	-- single pump PMCA flux (mol/s)
				+ ncxDensity *  6.7567567567567566e-23	-- single pump NCX flux (mol/s)
				- vdccDensity * 1.435114751437757757e-28  -- single channel VGCC flux (mol/s) at 1mM ext. Ca2+ and V_m = -0.07V
if leakPMconstant < 0 then error("PM leak flux is outward for these density settings!") end


-- firing pattern of the synapse
synSubset = 4
synStartTime = 0.0
caEntryDuration = 0.01

-- burst of calcium influx for active synapses (~1200 ions)
freq = 50      -- spike train frequency (Hz) (the ineq. 1/freq > caEntryDuration must hold)
nSpikes = 10   -- number of spikes	
function neumannBndCa(x, y, z, t, si)	
	-- spike train
	if (si == synSubset and t <= synStartTime + caEntryDuration + (nSpikes - 1) * 1.0/freq) then
        t = t % (1.0/freq)
	end --> now, treat like single spike
	
	-- single spike
	if 	(si == synSubset and t <= caEntryDuration)
	then influx = 0.01 * (1.0 - t/caEntryDuration)
	else influx = 0.0
	end
	
    return influx
end


-- burst of ip3 at active synapse (triangular, immediate)
ip3EntryDelay = 0.000
ip3EntryDuration = 0.2
function neumannBndIP3(x, y, z, t, si)
	if 	(si == synSubset and synStartTime+ip3EntryDelay < t and t <= synStartTime+ip3EntryDelay+ip3EntryDuration)
	then influx = 2e-3 * (1.0 - t/ip3EntryDuration)
	else influx = 0.0
	end
    return influx
end

-------------------------------
-- setup approximation space --
-------------------------------
-- create, load, refine and distribute domain
reqSubsets = {"cyt", "er", "pm", "erm", "syn", "meas_dend", "meas_neck", "meas_head"}
dom = util.CreateDomain(gridName, 0, reqSubsets)

balancer.partitioner = "parmetis"
balancer.staticProcHierarchy = true
balancer.firstDistLvl = -1
balancer.redistSteps = 0

-- in parallel environments: use a load balancer to distribute the grid
-- actual refinement and load balancing after setup of disc.
balancer.ParseParameters()
balancer.PrintParameters()
loadBalancer = balancer.CreateLoadBalancer(dom)

-- add distribution protection for ER membrane, then distribute
if loadBalancer ~= nil then
	if balancer.partitioner == "parmetis" then
		mu = ManifoldUnificator(dom)
		mu:add_protectable_subsets("erm")
		cdgm = ClusteredDualGraphManager()
		cdgm:add_unificator(SiblingUnificator())
		cdgm:add_unificator(mu)
		balancer.defaultPartitioner:set_dual_graph_manager(cdgm)
	end
	balancer.Rebalance(dom, loadBalancer)
end

if numRefs > 0 then	
	refiner = GlobalDomainRefiner(dom)
	for i = 1, numRefs do
		refiner:refine()
	end
end

if loadBalancer ~= nil then
	print("Edge cut on base level: "..balancer.defaultPartitioner:edge_cut_on_lvl(0))
	loadBalancer:estimate_distribution_quality()
	loadBalancer:print_quality_records()
end
print(dom:domain_info():to_string())


--[[
SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), "refined_grid_hierarchy_p" .. ProcRank() .. ".ugx", 2.0)
SaveParallelGridLayout(dom:grid(), "parallel_grid_layout_p"..ProcRank()..".ugx", 2.0)
--]]

-- create approximation space
approxSpace = ApproximationSpace(dom)

cytVol = "cyt"
measZones = "meas_head, meas_neck, meas_dend"
cytVol = cytVol .. ", " .. measZones
erVol = "er"
plMem = "pm, syn"
plMem_vec = {"pm", "syn"}
erMem = "erm"
erMemVec = {"erm"}

outerDomain = cytVol .. ", " .. plMem .. ", " .. erMem
innerDomain = erVol .. ", " .. erMem

approxSpace:add_fct("ca_er", "Lagrange", 1, innerDomain)
approxSpace:add_fct("ca_cyt", "Lagrange", 1, outerDomain)
approxSpace:add_fct("ip3", "Lagrange", 1, outerDomain)
approxSpace:add_fct("clb", "Lagrange", 1, outerDomain)
approxSpace:add_fct("o2", "Lagrange", 1, erMem)
approxSpace:add_fct("c1", "Lagrange", 1, erMem)
approxSpace:add_fct("c2", "Lagrange", 1, erMem)
approxSpace:add_fct("m", "Lagrange", 1, plMem)

approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_layout_statistic()
approxSpace:print_statistic()


--------------------------
-- setup discretization --
--------------------------
-- diffusion --
elemDiscCYT = ConvectionDiffusion("ca_cyt", cytVol, "fv1")
elemDiscCYT:set_diffusion(D_cac)

elemDiscER = ConvectionDiffusion("ca_er", erVol, "fv1") 
elemDiscER:set_diffusion(D_cae)

elemDiscIP3 = ConvectionDiffusion("ip3", cytVol, "fv1")
elemDiscIP3:set_diffusion(D_ip3)
elemDiscIP3:set_reaction_rate(reactionRateIP3)
elemDiscIP3:set_reaction(reactionTermIP3)

elemDiscClb = ConvectionDiffusion("clb", cytVol, "fv1")
elemDiscClb:set_diffusion(D_clb)

-- buffering --
elemDiscBuffering = BufferFV1(cytVol)	-- where buffering occurs
elemDiscBuffering:add_reaction(
	"clb",						    -- the buffering substance
	"ca_cyt",						-- the buffered substance
	totalClb,						-- total amount of buffer
	k_bind_clb,					    -- binding rate constant
	k_unbind_clb)				    -- unbinding rate constant


-- er membrane transport systems

-- We pass the function needed to evaluate the flux function here.
-- The order, in which the discrete fcts are passed, is crucial!
ip3r = IP3R({"ca_cyt", "ca_er", "ip3"})
ip3r:set_scale_inputs({1e3,1e3,1e3})
ip3r:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

ryr = RyRImplicit({"ca_cyt", "ca_er", "o2", "c1", "c2"}, erMemVec)
ryr:set_scale_inputs({1e3, 1e3, 1.0, 1.0, 1.0})
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
pmca:set_constant(1, ca_ext)
pmca:set_scale_inputs({1e3,1e3})
pmca:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

ncx = NCX({"ca_cyt", ""})
ncx:set_constant(1, ca_ext)
ncx:set_scale_inputs({1e3,1e3})
ncx:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

leakPM = Leak({"", "ca_cyt"})
leakPM:set_constant(0, ca_ext)
leakPM:set_scale_inputs({1e3,1e3})
leakPM:set_scale_fluxes({1e15}) -- from mol/(m^2 s) to (mol um)/(dm^3 s)

vdccMode = 0  -- 0 for script function voltage, 1 for data file voltage
if vdccMode == 0 then
	apSignal = ActionPotentialTrain(0.0, 0.02, 50, -70.0)
	function membranePotential(x, y, z, t, si)
		-- 80% of the intensity of the real AP
		return 1e-3*(-70.0 + 0.8*(apSignal:membrane_potential(t) + 70.0))
	end
	vdcc = VDCC_BG_UserData({"ca_cyt", "", "m"}, plMem_vec, approxSpace)
	vdcc:set_potential_function("membranePotential")
elseif vdccMode == 1 then
	vdcc = VDCC_BG_VM2UG({"ca_cyt", "", "m"}, plMem_vec, approxSpace,
			"voltageData/vm_", "%.5f", ".dat", false)
	vdcc:set_file_times(1e-5, 0.0) -- file interval is 0.01ms, starting at 0.0
end
vdcc:set_constant(1, ca_ext)
vdcc:set_scale_inputs({1e3, 1e3, 1.0})
vdcc:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)
vdcc:set_channel_type_L()
vdcc:init(0.0)


discPMCA = MembraneTransportFV1(plMem, pmca)
discPMCA:set_density_function(pmcaDensity)

discNCX = MembraneTransportFV1(plMem, ncx)
discNCX:set_density_function(ncxDensity)

discPMLeak = MembraneTransportFV1(plMem, leakPM)
discPMLeak:set_density_function(leakPMconstant / (1e3*(ca_ext - ca_cyt_init)))

discVDCC = MembraneTransportFV1(plMem, vdcc)
discVDCC:set_density_function(vdccDensity)


-- synaptic activity
synapseInfluxCa = UserFluxBoundaryFV1("ca_cyt", plMem)
synapseInfluxCa:set_flux_function("neumannBndCa")

synapseInfluxIP3 = UserFluxBoundaryFV1("ip3", plMem)
synapseInfluxIP3:set_flux_function("neumannBndIP3")


------------------------------------------
-- setup complete domain discretization --
------------------------------------------
domainDisc = DomainDiscretization(approxSpace)

domainDisc:add(elemDiscER)
domainDisc:add(elemDiscCYT)
domainDisc:add(elemDiscIP3)
domainDisc:add(elemDiscClb)

domainDisc:add(elemDiscBuffering)

domainDisc:add(discPMCA)
domainDisc:add(discNCX)
domainDisc:add(discPMLeak)
domainDisc:add(discVDCC)
domainDisc:add(vdcc)  -- for gating param discretization

domainDisc:add(discIP3R)
domainDisc:add(discRyR)
domainDisc:add(ryr) -- also add ryr as elem disc (for state variables)
domainDisc:add(discSERCA)
domainDisc:add(discERLeak)

domainDisc:add(synapseInfluxCa)
domainDisc:add(synapseInfluxIP3)

-- setup time discretization --
timeDisc = ThetaTimeStep(domainDisc)
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
    bcgs_steps = 1000
    ilu = ILU()
    ilu:set_sort(true)
    bcgs_precond = ilu
elseif (solverID == "GS") then
    bcgs_steps = 1000
    bcgs_precond = GaussSeidel()
else -- (solverID == "GMG")
	gmg = GeometricMultiGrid(approxSpace)
	gmg:set_discretization(timeDisc)
	gmg:set_base_level(0)
	gmg:set_gathered_base_solver_if_ambiguous(true)
	
	-- treat SuperLU problems with Dirichlet constraints by using constrained version
	gmg:set_base_solver(SuperLU())
	
	ilu_gmg = ILU()
	gmg:set_smoother(ilu_gmg)
	gmg:set_smooth_on_surface_rim(true)
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(3)
	gmg:set_num_postsmooth(3)
	gmg:set_rap(true)
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
newtonSolver = LimexNewtonSolver()
newtonSolver:set_linear_solver(bicgstabSolver)
--newtonSolver:set_debug(dbgWriter)

newtonSolver:init(op)

-------------
-- solving --
-------------
-- get grid function
u = GridFunction(approxSpace)

-- set initial value
Interpolate(ca_cyt_init, u, "ca_cyt", 0.0)
Interpolate(ca_er_init, u, "ca_er", 0.0)
Interpolate(ip3_init, u, "ip3", 0.0)
Interpolate(clb_init, u, "clb", 0.0)
vdcc:calculate_steady_state(u, -0.07)
ryr:calculate_steady_state(u)

-- timestep in seconds
dtmin = 1e-9
dtmax = 1e-2
time = 0.0
step = 0

-- initial vtk output
--if (generateVTKoutput) then
--	out = VTKOutput()
--	out:print(outDir .. "vtk/solution", u, step, time)
--end

-- initial measurements
take_measurement(u, time, measZones, "ca_cyt, ip3, clb", outDir .. "meas/data")

compute_volume(approxSpace, "meas_head, meas_neck, syn")


------------------
--  LIMEX setup --
------------------
nstages = 3              -- number of stages
stageNSteps = {1,2,3,4}  -- number of time steps for each stage

limex = LimexTimeIntegrator(nstages)
for i = 1, nstages do
	limex:add_stage(stageNSteps[i], newtonSolver, domainDisc)
end

limex:set_tolerance(tol)
limex:set_time_step(dt)
limex:set_dt_min(dtmin)
limex:set_dt_max(dtmax)
limex:set_increase_factor(2.0)
limex:set_reduction_factor(0.1)
limex:set_stepsize_greedy_order_factor(1)
limex:set_stepsize_safety_factor(0.25)

-- GridFunction error estimator (relative norm)
errorEvalCa = H1ComponentSpace("ca_cyt", "cyt", 3) -- function name, subset names
errorEvalCaE = H1ComponentSpace("ca_er", "er", 3) -- function name, subset names
errorEvalClb = H1ComponentSpace("clb", "cyt", 3) -- function name, subset names
errorEvalIP3 = H1ComponentSpace("ip3", "cyt", 3) -- function name, subset names
errorEvalC1 = L2ComponentSpace("c1", 3, 1.0, "erm") -- function name, subset names
errorEvalC2 = L2ComponentSpace("c2", 3, 1.0, "erm") -- function name, subset names
errorEvalO2 = L2ComponentSpace("o2", 3, 1.0, "erm") -- function name, subset names
limexEstimator = ScaledGridFunctionEstimator()
limexEstimator:add(errorEvalCa)
limexEstimator:add(errorEvalCaE)
limexEstimator:add(errorEvalClb)
limexEstimator:add(errorEvalIP3)
limexEstimator:add(errorEvalC1)
limexEstimator:add(errorEvalC2)
limexEstimator:add(errorEvalO2)
limex:add_error_estimator(limexEstimator)

-- for vtk output
if generateVTKoutput then 
	local vtkObserver = VTKOutputObserver(outDir .."vtk/solution", out, pstep)
	limex:attach_observer(vtkObserver)
end

-- for concentration measurements
function measurementCallback(step, time, dt)
	local curSol = measObserver:get_current_solution()

	take_measurement(curSol, time, measZones, "ca_cyt, ip3, clb", outDir .. "meas/data")

	return 0.0
end

measObserver = LuaCallbackObserver()
measObserver:set_callback("measurementCallback")
limex:attach_observer(measObserver)


-- solve problem
limex:apply(u, endTime, u, time)


--if generateVTKoutput then 
--	out:write_time_pvd(outDir .. "vtk/solution", u)
--end




