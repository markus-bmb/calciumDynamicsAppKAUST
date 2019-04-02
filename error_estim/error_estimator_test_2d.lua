----------------------------------------------------------------
--  Example script for error estimation	in 2D				  --
--                                                            --
--  Author: Markus Breit                                      --
----------------------------------------------------------------

-- for profiler output
SetOutputProfileStats(false)

-- load pre-implemented lua functions
ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

-- dimension
dim = 2

-- choose dimension and algebra
InitUG(dim, AlgebraType("CPU", 1));

-- choice of grid
gridName = util.GetParam("-grid", "calciumDynamics_app/grids/error_estimator_test_2d_long.ugx")
--gridName = "paper_test_2d.ugx"

-- total refinements
numRefs = util.GetParamNumber("-numRefs",    0)

-- choose length of maximal time step during the whole simulation
timeStep = util.GetParamNumber("-tstep", 0.01)

-- choose length of time step at the beginning
-- if not timeStepStart = 2^(-n)*timeStep, take nearest lower number of that form
timeStepStart = util.GetParamNumber("-tstepStart", timeStep)
function log2(x)
	return math.log(x)/math.log(2)
end
startLv =  math.ceil(log2(timeStep/timeStepStart))
timeStepStartNew = timeStep / math.pow(2, startLv)
if (math.abs(timeStepStartNew-timeStepStart)/timeStepStart > 1e-5) then 
	print("timeStepStart argument ("..timeStepStart..") was not admissible; taking "..timeStepStartNew.." instead.")
end
timeStepStart = timeStepStartNew
	
-- choose end time
nTimeSteps = util.GetParamNumber("-nTimeSteps", 1)
def_endTime = nTimeSteps*timeStep
endTime = util.GetParamNumber("-endTime", def_endTime)

-- choose plotting interval
plotStep = util.GetParamNumber("-pstep", timeStep)

-- choose solver setup
solverID = util.GetParam("-solver", "GS")
solverID = string.upper(solverID)
validSolverIDs = {}
validSolverIDs["GMG-GS"] = 0;
validSolverIDs["GMG-ILU"] = 0;
validSolverIDs["GMG-LU"] = 0;
validSolverIDs["GS"] = 0;
validSolverIDs["ILU"] = 0;
if (validSolverIDs[solverID] == nil) then
    error("Unknown solver identifier " .. solverID)
end
 
-- choose outfile directory
fileName = util.GetParam("-outName", "error_estimator")
fileName = fileName.."/"

-- specify -vtk to generate vtk output
generateVTKoutput = util.HasParamOption("-vtk")

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
ip3_init = 4.0e-8
clb_init = totalClb / (k_bind_clb/k_unbind_clb*ca_cyt_init + 1)

-- reaction reate IP3
reactionRateIP3 = 0.11

-- equilibrium concentration IP3
equilibriumIP3 = 4.0e-08

-- reation term IP3
reactionTermIP3 = -reactionRateIP3 * equilibriumIP3

---------------------------------------------------------------------
-- functions steering tempo-spatial parameterization of simulation --
---------------------------------------------------------------------
IP3Rdensity = 17.3
RYRdensity = 0.4
LEAKERconstant = 3.8e-17

-- this is a little bit more complicated, since it must be ensured that
-- the net flux for equilibrium concentrations is zero
-- MUST be adapted whenever any parameterization of ER flux mechanisms is changed!
local v_s = 6.5e-27						-- V_S param of SERCA pump
local k_s = 1.8e-7						-- K_S param of SERCA pump
local j_ip3r = 3.7606194166520605e-23 -- 2.7817352713488838e-23	-- single channel IP3R flux (mol/s) - to be determined via gdb
local j_ryr = 1.1204582669024472e-21 -- 4.6047720062808216e-22	-- single channel RyR flux (mol/s) - to be determined via gdb
local j_leak = ca_er_init-ca_cyt_init	-- leak proportionality factor

SERCAdensity =  IP3Rdensity * j_ip3r
			+ RYRdensity * j_ryr
			+ LEAKERconstant * j_leak
SERCAdensity = SERCAdensity / (v_s/(k_s/ca_cyt_init+1.0)/ca_er_init)


pmcaDensity = 500.0
ncxDensity  = 15.0
vgccDensity = 1.0

leakPMconstant =  pmcaDensity * 6.9672131147540994e-24	-- single pump PMCA flux (mol/s)
				+ ncxDensity *  6.7567567567567566e-23	-- single pump NCX flux (mol/s)
				+ vgccDensity * (-7.475181253921754074e-28)  -- single channel VGCC flux (mol/s)
				-- *1.5 // * 0.5 for L-type // T-type
if (leakPMconstant < 0) then error("PM leak flux is outward for these density settings!") end




-- firing pattern of the synapse
caEntryDuration = 0.01
syn_start = 0

-- burst of calcium influx for active synapses (~1200 ions)
freq = 50      -- spike train frequency (Hz) (the ineq. 1/freq > caEntryDuration must hold)
nSpikes = 10   -- number of spikes	
function ourNeumannBndCA(x, y, t, si)	
	-- spike train
	if t <= syn_start + caEntryDuration + nSpikes * 1.0/freq then
        t = t % (1.0/freq)
	end --> now, treat like single spike
	
	-- single spike
	if 	t>syn_start and t<=syn_start+caEntryDuration
	then efflux = -2e-3
	else efflux = 0.0
	end	
    
    --[[
 	-- more smoothness for the input signal
    if y < -0.5 then
		efflux = efflux*(2*(y+1)*(y+1))
	elseif y < 0.5 then 
		efflux = efflux*(1-2*y*y)
	else efflux = efflux*(2*(y-1)*(y-1))
	end
    --]]
    
    return -efflux
end


-- burst of ip3 at active synapse (triangular, immediate)
ip3EntryDelay = 0.000
ip3EntryDuration = 2.0
corrFact = -10.4
function ourNeumannBndIP3(x, y, t, si)
	if 	t>syn_start+ip3EntryDelay and t<=syn_start+ip3EntryDelay+ip3EntryDuration
	--then efflux = - math.exp(corrFact*t) * 2.1e-5/1.188 * (1.0 - (t-syn_start)/ip3EntryDuration)
	then efflux = - 2.1e-5 * (1.0 - (t-syn_start)/ip3EntryDuration)
	else efflux = 0.0
	end
	
	--[[
	-- more smoothness for the input signal
    if y < -0.5 then
		efflux = efflux*(2*(y+1)*(y+1))
	elseif y < 0.5 then 
		efflux = efflux*(1-2*y*y)
	else efflux = efflux*(2*(y-1)*(y-1))
	end
	--]]
	
    return -efflux
end

-------------------------------
-- setup approximation space --
-------------------------------
-- create, load, refine and distribute domain
dom = util.CreateDomain(gridName, numRefs)

balancer.partitioner = "parmetis"
balancer.staticProcHierarchy = -1
balancer.firstDistLvl = 0
balancer.firstDistProcs = 120
balancer.ParseParameters()

loadBalancer = balancer.CreateLoadBalancer(dom)
if loadBalancer ~= nil then
	cdgm = ClusteredDualGraphManager()
	mu = ManifoldUnificator(dom)
    cdgm:add_unificator(mu)
	balancer.defaultPartitioner:set_dual_graph_manager(cdgm)
	balancer.qualityRecordName = "coarse"
	balancer.Rebalance(dom, loadBalancer)
	
	edgeCut = balancer.defaultPartitioner:edge_cut_on_lvl(0)
	print("Edge cut on base level: " .. edgeCut)

	loadBalancer:estimate_distribution_quality()
	loadBalancer:print_quality_records()
end

--SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), fileName.."grid/refined_grid_hierarchy_p" .. ProcRank() .. ".ugx", 2.0)
--SaveParallelGridLayout(dom:grid(), "parallel_grid_layout_p"..ProcRank()..".ugx", 2.0)

-- create approximation space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)

cytVol = "cyt"

erVol = "er"

plMem = "mem_cyt, syn"
plMem_vec = {"mem_cyt", "syn"}

erMem = "mem_er"
measZonesERM = "measZoneERM"..1
for i=2,6 do
	measZonesERM = measZonesERM .. ", measZoneERM" .. i
end
erMem = erMem .. ", " .. measZonesERM

outerDomain = cytVol .. ", " .. plMem .. ", " .. erMem
innerDomain = erVol .. ", " .. erMem

approxSpace:add_fct("ca_er", "Lagrange", 1, innerDomain)
approxSpace:add_fct("ca_cyt", "Lagrange", 1, outerDomain)
approxSpace:add_fct("ip3", "Lagrange", 1, outerDomain)
approxSpace:add_fct("clb", "Lagrange", 1, outerDomain)

approxSpace:init_levels()
approxSpace:print_layout_statistic()
approxSpace:print_statistic()


----------------------------
-- setup error estimators --
----------------------------
eeCaCyt = SideAndElemErrEstData(4,4,cytVol)
eeCaER 	= SideAndElemErrEstData(4,4,erVol)
eeIP3 	= SideAndElemErrEstData(4,4,cytVol)
eeClb 	= SideAndElemErrEstData(4,4,cytVol)

eeMult = MultipleSideAndElemErrEstData(approxSpace)
eeMult:add(eeCaCyt, "ca_cyt")
eeMult:add(eeCaER, "ca_er")
eeMult:add(eeIP3, "ip3")
eeMult:add(eeClb, "clb")
eeMult:set_consider_me(false) -- not necessary (default)

----------------------------------------------------------
-- setup FV convection-diffusion element discretization --
----------------------------------------------------------
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


elemDiscCYT:set_error_estimator(eeCaCyt)
elemDiscER:set_error_estimator(eeCaER)
elemDiscIP3:set_error_estimator(eeIP3)
elemDiscClb:set_error_estimator(eeClb)

---------------------------------------
-- setup reaction terms of buffering --
---------------------------------------
elemDiscBuffering = BufferFV1(cytVol)	-- where buffering occurs
elemDiscBuffering:add_reaction(
	"clb",						    -- the buffering substance
	"ca_cyt",						-- the buffered substance
	totalClb,						-- total amount of buffer
	k_bind_clb,					    -- binding rate constant
	k_unbind_clb)				    -- unbinding rate constant

elemDiscBuffering:set_error_estimator(eeMult)

----------------------------------------------------
-- setup inner boundary (channels on ER membrane) --
----------------------------------------------------

-- We pass the function needed to evaluate the flux function here.
-- The order, in which the discrete fcts are passed, is crucial!
ip3r = IP3R({"ca_cyt", "ca_er", "ip3"})
ip3r:set_scale_inputs({1e3,1e3,1e3})
ip3r:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

ryr = RyR({"ca_cyt", "ca_er"})
ryr:set_scale_inputs({1e3,1e3})
ryr:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

serca = SERCA({"ca_cyt", "ca_er"})
serca:set_scale_inputs({1e3,1e3})
serca:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

leakER = Leak({"ca_er", "ca_cyt"})
leakER:set_scale_inputs({1e3,1e3})
leakER:set_scale_fluxes({1e3}) -- from mol/(m^2 s) to (mol um)/(dm^3 s)


innerDiscIP3R = MembraneTransportFV1(erMem, ip3r)
innerDiscIP3R:set_density_function(IP3Rdensity)

innerDiscRyR = MembraneTransportFV1(erMem, ryr)
innerDiscRyR:set_density_function(RYRdensity)

innerDiscSERCA = MembraneTransportFV1(erMem, serca)
innerDiscSERCA:set_density_function(SERCAdensity)

innerDiscLeak = MembraneTransportFV1(erMem, leakER)
innerDiscLeak:set_density_function(1e12*LEAKERconstant/(1e3)) -- from mol/(um^2 s M) to m/s


innerDiscIP3R:set_error_estimator(eeMult)
innerDiscRyR:set_error_estimator(eeMult)
innerDiscSERCA:set_error_estimator(eeMult)
innerDiscLeak:set_error_estimator(eeMult)

------------------------------
-- setup Neumann boundaries --
------------------------------
-- synaptic activity
---[[
neumannDiscCA = UserFluxBoundaryFV1("ca_cyt", "syn")
neumannDiscCA:set_flux_function("ourNeumannBndCA")
neumannDiscIP3 = UserFluxBoundaryFV1("ip3", "syn")
neumannDiscIP3:set_flux_function("ourNeumannBndIP3")

neumannDiscCA:set_error_estimator(eeMult)
neumannDiscIP3:set_error_estimator(eeMult)
--]]
--[[
neumannDiscCA = NeumannBoundary("ca_cyt")
neumannDiscCA:add("ourNeumannBndCA", plMem, cytVol)
neumannDiscIP3 = NeumannBoundary("ip3")
neumannDiscIP3:add("ourNeumannBndIP3", plMem, cytVol)

neumannDiscCA:set_error_estimator(eeCaCyt)
neumannDiscIP3:set_error_estimator(eeIP3)
--]]


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

vdcc = VDCC_BG_UserData({"ca_cyt", ""}, plMem_vec, approxSpace)
vdcc:set_constant(1, 1.0)
vdcc:set_scale_inputs({1e3,1.0})
vdcc:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)
vdcc:set_channel_type_L() --default, but to be sure
vdcc:set_potential_function(-0.07)
vdcc:init(0.0)

neumannDiscPMCA = MembraneTransportFV1(plMem, pmca)
neumannDiscPMCA:set_density_function(pmcaDensity)

neumannDiscNCX = MembraneTransportFV1(plMem, ncx)
neumannDiscNCX:set_density_function(ncxDensity)

neumannDiscLeak = MembraneTransportFV1(plMem, leakPM)
neumannDiscLeak:set_density_function(1e12*leakPMconstant / (1.0-1e3*ca_cyt_init))

neumannDiscVGCC = MembraneTransportFV1(plMem, vdcc)
neumannDiscVGCC:set_density_function(vgccDensity)


neumannDiscPMCA:set_error_estimator(eeMult)
neumannDiscNCX:set_error_estimator(eeMult)
neumannDiscLeak:set_error_estimator(eeMult)
neumannDiscVGCC:set_error_estimator(eeMult)

------------------------------------------
-- setup complete domain discretization --
------------------------------------------
domainDisc = DomainDiscretization(approxSpace)

-- diffusion discretizations
domainDisc:add(elemDiscER)
domainDisc:add(elemDiscCYT)
domainDisc:add(elemDiscIP3)
domainDisc:add(elemDiscClb)
--domainDisc:add(elemDiscClmC)
--domainDisc:add(elemDiscClmN)

-- buffering disc
domainDisc:add(elemDiscBuffering)
--domainDisc:add(elemDiscBuffering_clm)

-- (outer) boundary conditions
domainDisc:add(neumannDiscCA)
domainDisc:add(neumannDiscIP3)
domainDisc:add(neumannDiscPMCA)
domainDisc:add(neumannDiscNCX)
domainDisc:add(neumannDiscLeak)
domainDisc:add(neumannDiscVGCC)

-- ER flux
domainDisc:add(innerDiscIP3R)
domainDisc:add(innerDiscRyR)
domainDisc:add(innerDiscSERCA)
domainDisc:add(innerDiscLeak)


-- constraints for adatptivity
domainDisc:add(OneSideP1Constraints())

-------------------------------
-- setup time discretization --
-------------------------------
timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(1.0) -- 1.0 is implicit Euler

-- create operator from discretization
op = AssembledOperator()
op:set_discretization(timeDisc)
op:init()

------------------
-- solver setup --
------------------
-- create algebraic preconditioner
jac = Jacobi()
jac:set_damp(0.8)
gs = GaussSeidel()
sgs = SymmetricGaussSeidel()
bgs = BackwardGaussSeidel()
ilu = ILU()
ilut = ILUT()

-- exact solver
exactSolver = LU()


-- geometric multi-grid --
-- base solver
baseConvCheck = ConvCheck()
baseConvCheck:set_maximum_steps(1000)
baseConvCheck:set_minimum_defect(1e-28)
baseConvCheck:set_reduction(1e-1)
baseConvCheck:set_verbose(false)

-- debug writer
dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_vtk_output(false)

if (solverID == "GMG-LU") then
    base = exactSolver
else
    base = LinearSolver()
    base:set_convergence_check(baseConvCheck)
    if (solverID == "GMG-ILU") then
        base:set_preconditioner(ilu)
    else
        base:set_preconditioner(gs)
    end
end

-- gmg
gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(timeDisc)
gmg:set_base_level(0)
if (solverID == "GMG-LU") then
    gmg:set_gathered_base_solver_if_ambiguous(true)
end
gmg:set_base_solver(base)
if (solverID == "GMG-ILU") then
    gmg:set_smoother(ilu)
else
    gmg:set_smoother(ilu)
end 
gmg:set_cycle_type(1)
gmg:set_num_presmooth(5)
gmg:set_num_postsmooth(3)
--gmg:set_debug(dbgWriter)

-- biCGstab --
convCheck = ConvCheck()
convCheck:set_minimum_defect(1e-24)
convCheck:set_reduction(1e-06)
convCheck:set_verbose(false)
bicgstabSolver = BiCGStab()
if (solverID == "ILU") then
    convCheck:set_maximum_steps(2000)
    bicgstabSolver:set_preconditioner(ilu)
elseif (solverID == "GS") then
    convCheck:set_maximum_steps(2000)
    bicgstabSolver:set_preconditioner(gs)
else
    convCheck:set_maximum_steps(100)
    bicgstabSolver:set_preconditioner(gmg)
end
bicgstabSolver:set_convergence_check(convCheck)
--print(bicgstabSolver:config_string())

-----------------------
-- non linear solver --
-----------------------
-- convergence check
newtonConvCheck = CompositeConvCheck(approxSpace, 10, 1e-21, 1e-08)
newtonConvCheck:set_verbose(true)
newtonConvCheck:set_time_measurement(true)
newtonConvCheck:set_adaptive(true)

newtonLineSearch = StandardLineSearch()
newtonLineSearch:set_maximum_steps(5)
newtonLineSearch:set_accept_best(true)
newtonLineSearch:set_verbose(false)

-- Newton solver
newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(bicgstabSolver)
newtonSolver:set_convergence_check(newtonConvCheck)
--newtonSolver:set_line_search(newtonLineSearch)

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

-- timestep in seconds
dt = timeStepStart
time = 0.0
step = 0

if (generateVTKoutput) then
	out = VTKOutput()
	out:print(fileName .. "vtk/result", u, step, time)
end


-- refiner setup
refiner = HangingNodeDomainRefiner(dom)
TOL = 1e-15
maxLevel = 6
maxElem = 1000000
refStrat = StdRefinementMarking(TOL, maxLevel)
coarsStrat = StdCoarseningMarking(TOL)


-- set indicators for refinement in space and time to 0
space_refine_ind = 0.0
time_refine_ind = 0.0


outRefinement = VTKOutput()
approxSpace_vtk = ApproximationSpace(dom)
approxSpace_vtk:add_fct("eta_squared", "piecewise-constant");
u_vtk = GridFunction(approxSpace_vtk)
out_error = VTKOutput()
out_error:clear_selection()
out_error:select_all(false)
out_error:select_element("eta_squared", "error")

take_measurement(u, time, measZonesERM, "ca_cyt, ca_er, ip3, clb", fileName .. "meas/data")


-- create new grid function for old value
uOld = u:clone()

-- store grid function in vector of  old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)

--computeVolume(approxSpace, "cyt, er, mem_cyt, syn, mem_er")

newTime = true
min_dt = timeStep / math.pow(2,15)
cb_interval = 10
lv = startLv
levelUpDelay = caEntryDuration + (nSpikes - 1) * 1.0/freq;
cb_counter = {}
for i=0,startLv do cb_counter[i]=0 end
while endTime-time > 0.001*dt do
	if newTime then
		print("++++++ POINT IN TIME  " .. math.floor((time+dt)/dt+0.5)*dt .. "s  BEGIN ++++++")
		newTime = false
		numCoarsenOld = -1.0;
		n=0
	end
	
	-- setup time disc for old solutions and timestep
	timeDisc:prepare_step(solTimeSeries, dt)
	
	-- prepare Newton solver
	if newtonSolver:prepare(u) == false then print ("Newton solver failed at step "..step.."."); exit(); end 
	
	-- prepare BG channel state
	--[[ never update, so that always -72mV
	if (time+dt<0.2) then
		vm_time = math.floor((time+dt)/voltageFilesInterval)*voltageFilesInterval	-- truncate to last time that data exists for
		neumannDiscVGCC:update_potential(vm_time)
	end
	neumannDiscVGCC:update_gating(time+dt)
	--]]
	
	-- apply Newton solver
	newton_fail = false
	error_fail = false
	
	if newtonSolver:apply(u) == false then
		-- in case of Newton convergence failure:
		newton_fail = true
		print ("Newton solver failed at point in time " .. time .. " with time step " .. dt)
	else
		-- check error
		-- timeDisc:mark_error(new solution, refiner, tolerance for total squared l2 error,
		--					   refinement fraction of max error, coarsening fraction (-1) of min error,
		--					   max #refinements)
		timeDisc:calc_error(u, u_vtk)
		
		changedGrid = false
		
		-- refining
		domainDisc:mark_with_strategy(refiner, refStrat)
		if refiner:num_marked_elements() > 0 then
			error_fail = true
			print ("Error estimator is above required error.")
		end
	end
	
	if (error_fail or newton_fail)
	then
		print ("refinement score (time/space): " .. time_refine_ind .." / " .. space_refine_ind)
		if (time_refine_ind > space_refine_ind or not timeDisc:is_error_valid())
		then
			-- TIME STEP ADJUSTMENT ------------------------------------------------
			
			-- clear marks in the case where they might have been set (no newton_fail)
			refiner:clear_marks()
			
			-- correction for Borg-Graham channels: have to set back time
			--neumannDiscVGCC:update_gating(time)
			
			dt = dt/2
			lv = lv + 1
			VecScaleAssign(u, 1.0, solTimeSeries:latest())
			
			-- error is invalid, since time step has changed
			timeDisc:invalidate_error()
			
			-- halve time step and try again unless time step below minimum
			if dt < min_dt
			then 
				print ("Time step below minimum. Aborting. Failed at point in time " .. time .. ".")
				time = endTime
			else
				print ("Retrying with half the time step...")
				cb_counter[lv] = 0
			end
			-------------------------------------------------------------------
			time_refine_ind = time_refine_ind / 2
		else
			-- GRID REFINEMENT ------------------------------------------------
		
			-- if too many elements: 
			numElemBeforeRefinement = dom:domain_info():num_elements()
			if (numElemBeforeRefinement > maxElem) then
				print ("Adaptive refinement failed - too many elements. Aborting.")
				print ("Failed at point in time " .. time .. ".")
				time = endTime
			else
				if newton_fail then
					domainDisc:mark_with_strategy(refiner, refStrat)
				end
				adaptTolerance = 1.0;
				while (refiner:num_marked_elements() == 0) do
					adaptTolerance = adaptTolerance*0.1
					refStrat:set_tolerance(TOL*adaptTolerance)
					domainDisc:mark_with_strategy(refiner, refStrat)
				end
				refiner:refine()
				refiner:clear_marks()
				
				if (adaptTolerance < 1.0) then
					print ("Adaptive refinement tolerance has temporarily been reduced to "
						   .. TOL*adaptTolerance .. " in order to mark elements for refinement.")
				end
				
				numElemAfterRefinement = dom:domain_info():num_elements()
				space_refine_ind = space_refine_ind * numElemBeforeRefinement / numElemAfterRefinement;
				
				-- error is invalid, since grid has changed
				timeDisc:invalidate_error()
				
				--SaveDomain(dom, "refined_grid_".. math.floor((time+dt)/dt+0.5)*dt .. "_" .. n ..".ugx")
				n = n+1
				
				-- solve again on adapted grid
				VecScaleAssign(u, 1.0, solTimeSeries:latest())
				
				print ("Retrying with refined grid...")
			end
			-------------------------------------------------------------------
		end
	else
		-- GRID COARSENING ------------------------------------------------
		
		domainDisc:mark_with_strategy(refiner, coarsStrat)
		numElemBeforeCoarsening = dom:domain_info():num_elements()
		numCoarsenNew = refiner:num_marked_elements()
		if (numCoarsenNew >= numElemBeforeCoarsening/5
		   and (numCoarsenNew < 0.8*numCoarsenOld or numCoarsenOld < 0)) then
			refiner:coarsen()
			numCoarsenOld = numCoarsenNew
		end
		numElemAfterCoarsening = dom:domain_info():num_elements()
		
		refiner:clear_marks()
		
		-- grid is changed iff number of elements is different than before
		if (numElemAfterCoarsening ~= numElemBeforeCoarsening) then
			space_refine_ind = space_refine_ind * numElemBeforeCoarsening / numElemAfterCoarsening;
			
			-- error is invalid, since grid has changed
			timeDisc:invalidate_error()
			
			--SaveDomain(dom, "refined_grid_".. math.floor((time+dt)/dt+0.5)*dt .. "_" .. n ..".ugx")
			n = n+1
			
			-- solve again on adapted grid
			VecScaleAssign(u, 1.0, solTimeSeries:latest())
			
			print ("Error estimator is below required error.")
			print ("Retrying with coarsened grid...")
		-------------------------------------------------------------------
		else
			-- NORMAL PROCEEDING ----------------------------------------------
			
			-- on first success, init the refinement indicators
			if (time_refine_ind == 0) then
				time_refine_ind = 2.0	-- first refine will be in time, then in space
				space_refine_ind = 1.0			
			end
			
			numCoarsenOld = -1.0;
			-- update new time
			time = solTimeSeries:time(0) + dt
			newTime = true
			
			-- update check-back counter and if applicable, reset dt
			cb_counter[lv] = cb_counter[lv] + 1
			while cb_counter[lv] % (2*cb_interval) == 0 and lv > 0 and (time >= levelUpDelay or lv > startLv) do
				print ("Doubling time due to continuing convergence; now: " .. 2*dt)
				dt = 2*dt;
				lv = lv - 1
				cb_counter[lv] = cb_counter[lv] + cb_counter[lv+1] / 2
				cb_counter[lv+1] = 0
				time_refine_ind = time_refine_ind * 2
			end
			
			-- plot solution (& error) every plotStep seconds
			if (generateVTKoutput) then
				if math.abs(time/plotStep - math.floor(time/plotStep+0.5)) < 1e-5 then
					out:print(fileName .. "vtk/result", u, math.floor(time/plotStep+0.5), time)
					out_error:print(fileName .. "vtk/error_estimator_"..n, u_vtk, math.floor(time/plotStep+0.5), time)
				end
			end
			
			-- take measurement in nucleus every timeStep seconds 
			take_measurement(u, time, measZonesERM, "ca_cyt, ca_er, ip3, clb", fileName .. "meas/data")
			
			-- export solution of ca on mem_er
			--exportSolution(u, approxSpace, time, "mem_cyt", "ca_cyt", fileName .. "sol/sol");
			
			print("++++++ POINT IN TIME  " .. math.floor(time/dt+0.5)*dt .. "s  END ++++++++");
		
			-- get oldest solution
			oldestSol = solTimeSeries:oldest()
				
			-- copy values into oldest solution (we reuse the memory here)
			VecScaleAssign(oldestSol, 1.0, u)
			
			-- push oldest solutions with new values to front, oldest sol pointer is popped from end
			solTimeSeries:push_discard_oldest(oldestSol, time)
			-------------------------------------------------------------------
		end
	end
end

out_error:write_time_pvd(fileName .. "vtk/error_estimator", u_vtk)

-- end timeseries, produce gathering file
if (generateVTKoutput) then out:write_time_pvd(fileName .. "vtk/result", u)end
