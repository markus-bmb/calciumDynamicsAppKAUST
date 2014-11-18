----------------------------------------------------------------
--  This script calculates fluorescence values for a specific --
--	dye and activation pattern.								  --
--	It is intended to reproduce experimental results.		  --
--                                                            --
--  Author: Markus Breit                                      --
--    Date:	16-09-2014	                                      --
----------------------------------------------------------------

-- load pre-implemented lua functions
ug_load_script("ug_util.lua")

-- dimension
dim = 3

-- choose dimension and algebra
InitUG(dim, AlgebraType("CPU", 1));

-- speed up lua functions
--EnableLUA2C(true)
--SetDebugLevel(debugID.LUACompiler, 0) 

-- choice of grid
gridName = "camh36.ugx"

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
startTime = util.GetParamNumber("-startTime", 0.0)
nTimeSteps = util.GetParamNumber("-nTimeSteps", 1)
def_endTime = nTimeSteps*timeStep
endTime = util.GetParamNumber("-endTime", def_endTime)

-- choose plotting interval
plotStep = util.GetParamNumber("-pstep", 0.01)

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
fileName = util.GetParam("-outName", "fluorescence_data")
fileName = fileName.."/"

-- specify -vtk to generate vtk output
generateVTKoutput = util.HasParamOption("-vtk")

---------------
-- constants --
---------------
-- total cytosolic calbindin concentration
-- (four times the real value in order to simulate four binding sites in one)
totalf2ff = 200.0e-6

-- diffusion coefficients
D_cac = 220.0
D_cae = 220.0
D_ip3 = 280.0
D_f2ff = 20.0

-- Fura-2ff binding rates (Bollmann et al., 2000)
k_bind_f2ff = 	5.5e08
k_unbind_f2ff = 4.9e03

-- initial concentrations
ca_cyt_init = 5.0e-08 --4.0e-8
ca_er_init = 5.0e-4
ip3_init = 4.0e-8
f2ff_init = totalf2ff / (k_bind_f2ff/k_unbind_f2ff*ca_cyt_init + 1)

-- reaction reate IP3
reactionRateIP3 = 0.11

-- equilibrium concentration IP3
equilibriumIP3 = 4.0e-08

-- reation term IP3
reactionTermIP3 = -reactionRateIP3 * equilibriumIP3

---------------------------------------------------------------------
-- functions steering tempo-spatial parameterization of simulation --
---------------------------------------------------------------------

function IP3Rdensity(x,y,z,t,si)
	local d = z/120.0
	-- fourth order polynomial, distance to soma
	local dens = 1.4 -2.8*d +6.6*math.pow(d,2) -7.0*math.pow(d,3) +2.8*math.pow(d,4)
	dens = dens * 500.0--17.3

	-- extra ip3rs in fat part of dendrite
	if z < 30 then dens = dens*2 end
	return dens
end

function RYRdensity(x,y,z,t,si)
	local d = z/120.0
	-- fourth order polynomial, distance to soma
	local dens = 1.5 -3.5*d +9.1*math.pow(d,2) -10.5*math.pow(d,3) +4.3*math.pow(d,4)
	dens = 0.86
	return dens
end

-- this is a little bit more complicated, since it must be ensured that
-- the net flux for equilibrium concentrations is zero
-- MUST be adapted whenever any parameterization of ER flux mechanisms is changed!
function SERCAdensity(x,y,z,t,si)
	local v_s = 6.5e-27						-- V_S param of SERCA pump
	local k_s = 1.8e-7						-- K_S param of SERCA pump
	local j_ip3r = 7.52199e-23--3.7606194166520605e-23 -- 2.7817352713488838e-23	-- single channel IP3R flux (mol/s) - to be determined via gdb
	local j_ryr = 2.24114e-21--1.1204582669024472e-21 -- 4.6047720062808216e-22	-- single channel RyR flux (mol/s) - to be determined via gdb
	local j_leak = ca_er_init-ca_cyt_init	-- leak proportionality factor
	
	local dens =  IP3Rdensity(x,y,z,t,si) * j_ip3r
				+ RYRdensity(x,y,z,t,si) * j_ryr
				+ LEAKERconstant(x,y,z,t,si) * j_leak
	dens = dens / (v_s/(k_s/ca_cyt_init+1.0)/ca_er_init)
	
	return dens
end

function LEAKERconstant(x,y,z,t,si)
	return 3.3e-17
end

function PMCAdensity(x,y,z,t,si)
	return 1000.0
end

function NCXdensity(x,y,z,t,si)
	return 30.0
end

function VGCCdensity(x,y,z,t,si)
	return 1.0
end

function LEAKPMconstant(x,y,z,t,si)
	local j_pmca = - 6.9672131147540994e-24 -- - 5.230769230769231e-24	-- single pump PMCA flux (mol/s) - to be determined via gdb
	local j_ncx = - 6.7567567567567566e-23 -- - 5.4347826086956515e-23	-- single pump NCX flux (mol/s) - to be determined via gdb
	local j_vgcc = 1.5752042094823713e-25	-- single channel VGCC flux (mol/s) - to be determined via gdb
				-- *1.5 // * 0.5 for L-type // T-type
	local flux =  PMCAdensity(x,y,z,t,si) * j_pmca
				+ NCXdensity(x,y,z,t,si) * j_ncx
				+ VGCCdensity(x,y,z,t,si) * j_vgcc
	
	if (-flux < 0) then error("PM leak flux is outward for these density settings!") end
	return -flux
end


-- firing pattern of the synapse
caEntryDuration = 0.01
syn_start = {}
syn_start[4] = 1.2
syn_start[5] = 1.3
syn_start[6] = 1.4

-- burst of calcium influx for active synapses (~1200 ions)
freq = 50      -- spike train frequency (Hz) (the ineq. 1/freq > caEntryDuration must hold)
nSpikes = 10   -- number of spikes

function ourNeumannBndCA(x, y, z, t, si)	
	-- spike train
	if (si==4 and t>syn_start[4] and t <= syn_start[4] + caEntryDuration + (nSpikes-1) * 1.0/freq) then
        t = syn_start[4] + (t-syn_start[4]) % (1.0/freq)
    elseif (si==5 and t>syn_start[5] and t <= syn_start[5] + caEntryDuration + (nSpikes-1) * 1.0/freq) then
        t = syn_start[5] + (t-syn_start[5]) % (1.0/freq)
    elseif (si==6 and t>syn_start[6] and t <= syn_start[6] + caEntryDuration + (nSpikes-1) * 1.0/freq) then
        t = syn_start[6] + (t-syn_start[6]) % (1.0/freq)
	end --> now, treat like single spike
	
	-- single spike
	if 	(si==4 and t>syn_start[4] and t<=syn_start[4]+caEntryDuration)
		or (si==5 and t>syn_start[5] and t<=syn_start[5]+caEntryDuration)
		or (si==6 and t>syn_start[6] and t<=syn_start[6]+caEntryDuration)
	then efflux = 2e-4
	else efflux = 0.0
	end	
    return efflux
end


-- burst of ip3 at active synapse (triangular, immediate)
ip3EntryDelay = 0.000
ip3EntryDuration = 2.0
function ourNeumannBndIP3(x, y, z, t, si)
	---[[
	if 	(si==4 and t>syn_start[4]+ip3EntryDelay and t<=syn_start[4]+ip3EntryDelay+ip3EntryDuration)
	then efflux = 2.1e-5/1.36 * (1.0 - (t-syn_start[4])/ip3EntryDuration)
	elseif 	(si==5 and t>syn_start[5]+ip3EntryDelay and t<=syn_start[5]+ip3EntryDelay+ip3EntryDuration)
	then efflux = 2.1e-5/1.36 * (1.0 - (t-syn_start[5])/ip3EntryDuration)
	elseif 	(si==6 and t>syn_start[6]+ip3EntryDelay and t<=syn_start[6]+ip3EntryDelay+ip3EntryDuration)
	then efflux = 2.1e-5/1.36 * (1.0 - (t-syn_start[6])/ip3EntryDuration)
	else efflux = 0.0
	end
    return efflux
	--]]
	--return true, 0.0
end

-------------------------------
-- setup approximation space --
-------------------------------
-- create, load, refine and distribute domain
print("create, refine and distribute domain")
neededSubsets = {"cyt","er","mem_cyt","mem_er"}
--[[
distributionMethod = "bisection"
dom = util.CreateAndDistributeDomain(gridName, numRefs, 0, neededSubsets, distributionMethod)
--]]
---[[
distributionMethod = "metisReweigh"
weightingFct = InterSubsetPartitionWeighting()
weightingFct:set_default_weights(1,1)
weightingFct:set_inter_subset_weight(0, 1, 1000)
dom = util.CreateAndDistributeDomain(gridName, numRefs, 0, neededSubsets, distributionMethod, nil, nil, nil, weightingFct)
--]]
--[[
--print("Saving domain grid and hierarchy.")
--SaveDomain(dom, "refined_grid_p" .. ProcRank() .. ".ugx")
--SaveGridHierarchyTransformed(dom:grid(), "refined_grid_hierarchy_p" .. ProcRank() .. ".ugx", 20.0)
print("Saving parallel grid layout")
SaveParallelGridLayout(dom:grid(), "parallel_grid_layout_p"..ProcRank()..".ugx", 20.0)
--]]

-- create approximation space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)

cytVol = "cyt"

erVol = "er"

plMem = "mem_cyt, syn1, syn2, syn3"

erMem = "mem_er"
--[[
measZonesERM = "measZoneERM"..1
for i=2,6 do
	measZonesERM = measZonesERM .. ", measZoneERM" .. i
end
erMem = erMem .. ", " .. measZonesERM
--]]

outerDomain = cytVol .. ", " .. plMem .. ", " .. erMem
innerDomain = erVol .. ", " .. erMem

approxSpace:add_fct("ca_er", "Lagrange", 1, innerDomain)
approxSpace:add_fct("ca_cyt", "Lagrange", 1, outerDomain)
approxSpace:add_fct("ip3", "Lagrange", 1, outerDomain)
approxSpace:add_fct("f2ff", "Lagrange", 1, outerDomain)

approxSpace:init_levels()
approxSpace:print_layout_statistic()
approxSpace:print_statistic()

----------------------------------------------------------
-- setup FV convection-diffusion element discretization --
----------------------------------------------------------
print ("Setting up Assembling")

if dim == 2 then 
    upwind = NoUpwind2d()
elseif dim == 3 then 
    upwind = NoUpwind3d()
end

elemDiscER = ConvectionDiffusion("ca_er", erVol, "fv1") 
elemDiscER:set_diffusion(D_cae)
elemDiscER:set_upwind(upwind)

elemDiscCYT = ConvectionDiffusion("ca_cyt", cytVol, "fv1")
elemDiscCYT:set_diffusion(D_cac)
elemDiscCYT:set_upwind(upwind)

elemDiscIP3 = ConvectionDiffusion("ip3", cytVol, "fv1")
elemDiscIP3:set_diffusion(D_ip3)
elemDiscIP3:set_reaction_rate(reactionRateIP3)
elemDiscIP3:set_reaction(reactionTermIP3)
elemDiscIP3:set_upwind(upwind)

elemDiscf2ff = ConvectionDiffusion("f2ff", cytVol, "fv1")
elemDiscf2ff:set_diffusion(D_f2ff)
elemDiscf2ff:set_upwind(upwind)

---------------------------------------
-- setup reaction terms of buffering --
---------------------------------------
elemDiscBuffering = BufferFV1(cytVol)	-- where buffering occurs
elemDiscBuffering:add_reaction(
	"f2ff",						    -- the buffering substance
	"ca_cyt",						-- the buffered substance
	totalf2ff,						-- total amount of buffer
	k_bind_f2ff,					    -- binding rate constant
	k_unbind_f2ff)				    -- unbinding rate constant

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

------------------------------
-- setup Neumann boundaries --
------------------------------
-- synaptic activity
neumannDiscCA = UserFluxBoundaryFV1("ca_cyt", plMem)
neumannDiscCA:set_flux_function("ourNeumannBndCA")
neumannDiscIP3 = UserFluxBoundaryFV1("ip3", plMem)
neumannDiscIP3:set_flux_function("ourNeumannBndIP3")

-- plasme membrane transport systems
neumannDiscPMCA = OneSidedPMCAFV1("ca_cyt", plMem)
neumannDiscPMCA:set_density_function("PMCAdensity")

neumannDiscNCX = OneSidedNCXFV1("ca_cyt", plMem)
neumannDiscNCX:set_density_function("NCXdensity")

neumannDiscLeak = OneSidedPMCalciumLeakFV1("ca_cyt", plMem)
neumannDiscLeak:set_density_function("LEAKPMconstant")

neumannDiscVGCC = OneSidedBorgGrahamFV1WithVM2UG("ca_cyt", plMem, approxSpace,
		"fluorescence_data/camh36/Vm/Vm_", "%.4f", ".dat", false)
neumannDiscVGCC:set_channel_type_L() --default, but to be sure
neumannDiscVGCC:set_density_function("VGCCdensity")

firstVoltageFile = 0.9999
lastVoltageFile = 1.202
voltageFilesInterval = 0.0001;

neumannDiscVGCC:init(firstVoltageFile)

------------------------------------------
-- setup complete domain discretization --
------------------------------------------
domainDisc = DomainDiscretization(approxSpace)

-- diffusion discretizations
domainDisc:add(elemDiscER)
domainDisc:add(elemDiscCYT)
domainDisc:add(elemDiscIP3)
domainDisc:add(elemDiscf2ff)

-- buffering disc
domainDisc:add(elemDiscBuffering)


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
--domainDisc:add(OneSideP1Constraints())

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

-- gmg --
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
gmg:set_rap(false)
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
newtonConvCheck = CompositeConvCheck3dCPU1(approxSpace, 10, 1e-16, 1e-08)
--newtonConvCheck:set_all_component_check(1e-16, 1e-08)
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
--newtonSolver:set_debug(GridFunctionDebugWriter(approxSpace))

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
Interpolate(f2ff_init, u, "f2ff", 0.0)

-- timestep in seconds
dt = timeStepStart
time = startTime
step = 0

if (generateVTKoutput) then
	out = VTKOutput()
	out:print(fileName .. "vtk/fluorescence", u, step, time)
end

--takeMeasurement(u, time, measZonesERM, "ca_cyt, ca_er, ip3, f2ff", fileName .. "meas/data")

-- create new grid function for old value
uOld = u:clone()

-- store grid function in vector of  old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)

min_dt = timeStep / math.pow(2,15)
cb_interval = 10
lv = startLv
levelUpDelay = caEntryDuration + (nSpikes - 1) * 1.0/freq;
cb_counter = {}
for i=0,startLv do cb_counter[i]=0 end
while endTime-time > 0.001*dt do
	print("++++++ POINT IN TIME  " .. math.floor((time+dt)/dt+0.5)*dt .. "s  BEGIN ++++++")
	
	-- setup time disc for old solutions and timestep
	timeDisc:prepare_step(solTimeSeries, dt)
	
	-- prepare Newton solver
	if newtonSolver:prepare(u) == false then print ("Newton solver failed at step "..step.."."); exit(); end 
	
	-- prepare BG channel state
	if time+dt >= firstVoltageFile and time+dt <= lastVoltageFile then
		vm_time = math.floor((time+dt)/voltageFilesInterval)*voltageFilesInterval	-- truncate to last time that data exists for
		neumannDiscVGCC:update_potential(vm_time)	
	end
	neumannDiscVGCC:update_gating(time+dt)
	
	-- apply Newton solver
	newton_fail = false
	error_fail = false
	
	if newtonSolver:apply(u) == false then
		-- in case of Newton convergence failure:
		print ("Newton solver failed at point in time " .. time .. " with time step " .. dt)
		
		-- correction for Borg-Graham channels: have to set back time
		neumannDiscVGCC:update_gating(time)
		
		dt = dt/2
		lv = lv + 1
		VecScaleAssign(u, 1.0, solTimeSeries:latest())
		
		-- halve time step and try again unless time step below minimum
		if dt < min_dt
		then 
			print ("Time step below minimum. Aborting. Failed at point in time " .. time .. ".")
			time = endTime
		else
			print ("Retrying with half the time step...")
			cb_counter[lv] = 0
		end
	else
		-- update new time
		time = solTimeSeries:time(0) + dt
		
		-- update check-back counter and if applicable, reset dt
		cb_counter[lv] = cb_counter[lv] + 1
		while cb_counter[lv] % (2*cb_interval) == 0 and lv > 0 and (time >= levelUpDelay or lv > startLv) do
			print ("Doubling time due to continuing convergence; now: " .. 2*dt)
			dt = 2*dt;
			lv = lv - 1
			cb_counter[lv] = cb_counter[lv] + cb_counter[lv+1] / 2
			cb_counter[lv+1] = 0
		end
		
		-- ensure that time step is not bigger than voltage file resolution when stimulation occurs
		if time >= firstVoltageFile and time <= lastVoltageFile then
			while dt > voltageFilesInterval do
				dt = dt/2
				lv = lv + 1
				cb_counter[lv] = 0
			end
		end
		
		-- plot solution every plotStep seconds
		if (generateVTKoutput) then
			if math.abs((time-startTime)/plotStep - math.floor((time-startTime)/plotStep+0.5)) < 1e-5 then
				out:print(fileName .. "vtk/fluorescence", u, math.floor((time-startTime)/plotStep+0.5), time)
			end
		end
		
		-- take measurement in nucleus every timeStep seconds 
		--takeMeasurement(u, time, measZonesERM, "ca_cyt, ca_er, ip3, f2ff", fileName .. "meas/data")
		
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

-- end timeseries, produce gathering file
if (generateVTKoutput) then out:write_time_pvd(fileName .. "vtk/fluorescence", u) end
