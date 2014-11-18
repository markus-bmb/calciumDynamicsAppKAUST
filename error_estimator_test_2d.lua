----------------------------------------------------------------
--  Example script for error estimation	in 2D				  --
--                                                            --
--  Author: Markus Breit                                      --
----------------------------------------------------------------

-- for profiler output
SetOutputProfileStats(false)

-- load pre-implemented lua functions
ug_load_script("ug_util.lua")

-- dimension
dim = 2

-- choose dimension and algebra
InitUG(dim, AlgebraType("CPU", 1));

-- choice of grid
gridName = "error_estimator_test_2d.ugx"
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
function CaCytStart(x, y, t)
    return ca_cyt_init
end

function CaERStart(x, y, t)
    return ca_er_init
end

function IP3Start(x, y, t)
    return ip3_init
end

function clbStart(x, y, t)
    return clb_init
end

function ourDiffTensorCAcyt(x, y, t)
    return	D_cac, 0,
            0, D_cac
end

function ourDiffTensorCAer(x, y, t)
    return	D_cae, 0,
            0, D_cae
end

function ourDiffTensorIP3(x, y, t)
    return	D_ip3, 0,
            0, D_ip3
end

function ourDiffTensorClb(x, y, t)
    return	D_clb, 0,
            0, D_clb
end

function ourRhs(x, y, t)
    return 0;
end


function IP3Rdensity(x,y,t,si)
	return 17.3
end

function RYRdensity(x,y,t,si)
	return 0.4--0.86
end

-- this is a little bit more complicated, since it must be ensured that
-- the net flux for equilibrium concentrations is zero
-- MUST be adapted whenever any parameterization of ER flux mechanisms is changed!
function SERCAdensity(x,y,t,si)
	local v_s = 6.5e-27						-- V_S param of SERCA pump
	local k_s = 1.8e-7						-- K_S param of SERCA pump
	local j_ip3r = 3.7606194166520605e-23 -- 2.7817352713488838e-23	-- single channel IP3R flux (mol/s) - to be determined via gdb
	local j_ryr = 1.1204582669024472e-21 -- 4.6047720062808216e-22	-- single channel RyR flux (mol/s) - to be determined via gdb
	local j_leak = ca_er_init-ca_cyt_init	-- leak proportionality factor
	
	local dens =  IP3Rdensity(x,y,z,t,si) * j_ip3r
				+ RYRdensity(x,y,z,t,si) * j_ryr
				+ LEAKERconstant(x,y,z,t,si) * j_leak
	dens = dens / (v_s/(k_s/ca_cyt_init+1.0)/ca_er_init)
	
	return dens
end

function LEAKERconstant(x,y,t,si)
	return 3.8e-17 --3.4e-17
end

function PMCAdensity(x,y,t,si)
	return 500.0
end

function NCXdensity(x,y,t,si)
	return 15.0
end

function VGCCdensity(x,y,t,si)
	return 1.0
end

function LEAKPMconstant(x,y,t,si)
	local j_pmca = - 6.9672131147540994e-24 -- - 5.230769230769231e-24	-- single pump PMCA flux (mol/s) - to be determined via gdb
	local j_ncx = - 6.7567567567567566e-23 -- - 5.4347826086956515e-23	-- single pump NCX flux (mol/s) - to be determined via gdb
	local j_vgcc = 1.5752042094823713e-25	-- single channel VGCC flux (mol/s) - to be determined via gdb
				-- *1.5 // * 0.5 for L-type // T-type
	local flux =  PMCAdensity(x,y,z,t,si) * j_pmca
				+ NCXdensity(x,y,z,t,si) * j_ncx
				+ VGCCdensity(x,y,z,t,si) * j_vgcc
	
	if (-flux < 0) then error("PM leak flux is outward for these density settings!") end
	return -flux -- 6.85e-22
end



-- firing pattern of the synapse
caEntryDuration = 0.01
syn_start = 0

-- burst of calcium influx for active synapses (~1200 ions)
freq = 50      -- spike train frequency (Hz) (the ineq. 1/freq > caEntryDuration must hold)
nSpikes = 10   -- number of spikes	
function ourNeumannBndCA(x, y, t, si)	
	-- spike train
	if (si==4 and t <= syn_start + caEntryDuration + nSpikes * 1.0/freq) then
        t = t % (1.0/freq)
	end --> now, treat like single spike
	
	-- single spike
	if 	(si==4 and t>syn_start and t<=syn_start+caEntryDuration)
	then efflux = -2e-4
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
	if 	(si==4 and t>syn_start+ip3EntryDelay and t<=syn_start+ip3EntryDelay+ip3EntryDuration)
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
print("create, refine and distribute domain")
neededSubsets = {}
distributionMethod = "metisReweigh"
weightingFct = InterSubsetPartitionWeighting()
weightingFct:set_default_weights(1,1)
weightingFct:set_inter_subset_weight(0, 1, 1000)
dom = util.CreateAndDistributeDomain(gridName, numRefs, 0, neededSubsets, distributionMethod, nil, nil, nil, weightingFct)

---[[
--print("Saving domain grid and hierarchy.")
SaveDomain(dom, "refined_grid_p" .. ProcRank() .. ".ugx")
SaveGridHierarchyTransformed(dom:grid(), "refined_grid_hierarchy_p" .. ProcRank() .. ".ugx", 2.0)
--print("Saving parallel grid layout")
SaveParallelGridLayout(dom:grid(), "parallel_grid_layout_p"..ProcRank()..".ugx", 2.0)
--]]

-- create approximation space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)

cytVol = "cyt"

erVol = "er"

plMem = "mem_cyt, syn"

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

--------------------------
-- setup user functions --
--------------------------
print ("Setting up Assembling")

-- start value function setup
CaCytStartValue = LuaUserNumber2d("CaCytStart")
CaERStartValue = LuaUserNumber2d("CaERStart")
IP3StartValue = LuaUserNumber2d("IP3Start")
ClbStartValue = LuaUserNumber2d("clbStart")

-- diffusion Tensor setup
diffusionMatrixCAcyt = LuaUserMatrix2d("ourDiffTensorCAcyt")
diffusionMatrixCAer = LuaUserMatrix2d("ourDiffTensorCAer")
diffusionMatrixIP3 = LuaUserMatrix2d("ourDiffTensorIP3")
diffusionMatrixClb = LuaUserMatrix2d("ourDiffTensorClb")

-- rhs setup
rhs = LuaUserNumber2d("ourRhs")


----------------------------------------------------------
-- setup FV convection-diffusion element discretization --
----------------------------------------------------------
-- Note: No VelocityField and Reaction is set. The assembling assumes default
--       zero values for them

if dim == 2 then 
    upwind = NoUpwind2d()
elseif dim == 3 then 
    upwind = NoUpwind3d()
end

elemDiscCYT = ConvectionDiffusion("ca_cyt", cytVol, "fv1")
elemDiscCYT:set_diffusion(diffusionMatrixCAcyt)
elemDiscCYT:set_upwind(upwind)

elemDiscER = ConvectionDiffusion("ca_er", erVol, "fv1") 
elemDiscER:set_diffusion(diffusionMatrixCAer)
elemDiscER:set_upwind(upwind)

elemDiscIP3 = ConvectionDiffusion("ip3", cytVol, "fv1")
elemDiscIP3:set_diffusion(diffusionMatrixIP3)
elemDiscIP3:set_reaction_rate(reactionRateIP3)
elemDiscIP3:set_reaction(reactionTermIP3)
elemDiscIP3:set_upwind(upwind)

elemDiscClb = ConvectionDiffusion("clb", cytVol, "fv1")
elemDiscClb:set_diffusion(diffusionMatrixClb)
elemDiscClb:set_upwind(upwind)

-- error estimators
eeCaCyt = SideAndElemErrEstData(4, 4, cytVol)
eeCaER 	= SideAndElemErrEstData(4, 4, erVol)
eeIP3 	= SideAndElemErrEstData(4, 4, cytVol)
eeClb 	= SideAndElemErrEstData(4, 4, cytVol)

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

-- error estimator
eeBuffering = MultipleSideAndElemErrEstData()
eeBuffering:add(eeClb)
eeBuffering:add(eeCaCyt)
eeBuffering:set_consider_me(false)

elemDiscBuffering:set_error_estimator(eeBuffering)

----------------------------------------------------
-- setup inner boundary (channels on ER membrane) --
----------------------------------------------------

-- We pass the function needed to evaluate the flux function here.
-- The order, in which the discrete fcts are passed, is crucial!
innerDiscIP3R = TwoSidedIP3RFV1("ca_cyt, ca_er, ip3", erMem)
innerDiscIP3R:set_density_function("IP3Rdensity")

innerDiscRyR = TwoSidedRyRFV1("ca_cyt, ca_er", erMem)
innerDiscRyR:set_density_function("RYRdensity")

innerDiscSERCA = TwoSidedSERCAFV1("ca_cyt, ca_er", erMem)
innerDiscSERCA:set_density_function("SERCAdensity")

innerDiscLeak = TwoSidedERCalciumLeakFV1("ca_cyt, ca_er", erMem)
innerDiscLeak:set_density_function("LEAKERconstant")

-- error estimators
eeERM = MultipleSideAndElemErrEstData()
eeERM:add(eeCaCyt)
eeERM:add(eeCaER)
eeERM:add(eeIP3)
eeERM:set_consider_me(false)

innerDiscIP3R:set_error_estimator(eeERM)
innerDiscRyR:set_error_estimator(eeERM)
innerDiscSERCA:set_error_estimator(eeERM)
innerDiscLeak:set_error_estimator(eeERM)

------------------------------
-- setup Neumann boundaries --
------------------------------
-- synaptic activity
---[[
neumannDiscCA = UserFluxBoundaryFV1("ca_cyt", plMem)
neumannDiscCA:set_flux_function("ourNeumannBndCA")
neumannDiscIP3 = UserFluxBoundaryFV1("ip3", plMem)
neumannDiscIP3:set_flux_function("ourNeumannBndIP3")
--]]
--[[
neumannDiscCA = NeumannBoundary("ca_cyt")
neumannDiscCA:add("ourNeumannBndCA", plMem, cytVol)
neumannDiscIP3 = NeumannBoundary("ip3")
neumannDiscIP3:add("ourNeumannBndIP3", plMem, cytVol)
--]]
---[[
eeNeumannCA = MultipleSideAndElemErrEstData()
eeNeumannCA:add(eeCaCyt)
eeNeumannCA:set_consider_me(false)
neumannDiscCA:set_error_estimator(eeNeumannCA)
eeNeumannIP3 = MultipleSideAndElemErrEstData()
eeNeumannIP3:add(eeIP3)
eeNeumannIP3:set_consider_me(false)
neumannDiscIP3:set_error_estimator(eeNeumannIP3)
--]]
--[[
neumannDiscCA:set_error_estimator(eeCaCyt)
neumannDiscIP3:set_error_estimator(eeIP3)
--]]


-- plasme membrane transport systems
neumannDiscPMCA = OneSidedPMCAFV1("ca_cyt", plMem)
neumannDiscPMCA:set_density_function("PMCAdensity")

neumannDiscNCX = OneSidedNCXFV1("ca_cyt", plMem)
neumannDiscNCX:set_density_function("NCXdensity")

neumannDiscLeak = OneSidedPMCalciumLeakFV1("ca_cyt", plMem)
neumannDiscLeak:set_density_function("LEAKPMconstant")

neumannDiscVGCC = OneSidedBorgGrahamFV1WithVM2UG("ca_cyt", plMem, approxSpace,
		"neuronRes/timestep".."_order".. 0 .."_jump"..string.format("%1.1f", 5.0).."_", "%.3f", ".dat", false)
neumannDiscVGCC:set_channel_type_L() --default, but to be sure
neumannDiscVGCC:set_density_function("VGCCdensity")
neumannDiscVGCC:init(0.0)

voltageFilesInterval = 0.001;

-- error estimators
eePM = MultipleSideAndElemErrEstData()
eePM:add(eeCaCyt)
eePM:set_consider_me(false)


neumannDiscPMCA:set_error_estimator(eePM)
neumannDiscNCX:set_error_estimator(eePM)
neumannDiscLeak:set_error_estimator(eePM)
neumannDiscVGCC:set_error_estimator(eePM)

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
Interpolate(CaCytStartValue, u, "ca_cyt", 0.0)
Interpolate(CaERStartValue, u, "ca_er", 0.0)
Interpolate(IP3StartValue, u, "ip3", 0.0)
Interpolate(ClbStartValue, u, "clb", 0.0)

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
refineFrac = 0.01
coarseFrac = 0.9
maxLevel = 6
maxElem = 50000

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

takeMeasurement(u, time, measZonesERM, "ca_cyt, ca_er, ip3, clb", fileName .. "meas/data")


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
		timeDisc:mark_for_refinement(refiner, TOL, refineFrac, maxLevel)
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
					timeDisc:mark_for_refinement(refiner, TOL, refineFrac, maxLevel)
				end
				adaptTolerance = 1.0;
				while (refiner:num_marked_elements() == 0) do
					adaptTolerance = adaptTolerance*0.1
					timeDisc:mark_for_refinement(refiner, TOL*adaptTolerance, refineFrac, maxLevel)
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
		
		timeDisc:mark_for_coarsening(refiner, TOL, coarseFrac, maxLevel)
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
			takeMeasurement(u, time, measZonesERM, "ca_cyt, ca_er, ip3, clb", fileName .. "meas/data")
			
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
