--------------------------------------------------------------
--  Example script for simulation on 3d spine model			--
--                                                          --
--  Author: Markus Breit / Marcus Kessler                                   --
--------------------------------------------------------------

-- for profiler output
SetOutputProfileStats(false)

-- load pre-implemented lua functions
ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

-- dimension
dim = 3

-- choose dimension and algebra
InitUG(dim, AlgebraType("CPU", 1));

-- choice of grid
gridName = util.GetParam("-grid", "dend.ugx")

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
endTime = util.GetParamNumber("-endTime")
if (endTime == nil)
then
	-- choose number of time steps
	nTimeSteps = util.GetParamNumber("-nTimeSteps", 1)
	endTime = nTimeSteps*timeStep
end

-- choose plotting interval
plotStep = util.GetParamNumber("-pstep", 0.01)

-- choose outfile directory
fileName = util.GetParam("-outName", "spine/reconstructed")
fileName = fileName.."/"

-- specify -vtk to generate vtk output
generateVTKoutput = util.HasParamOption("-vtk")

---------------
-- constants --
---------------

-- diffusion coefficients
D_chl = 2000.0

-- initial concentrations
chl_cyt_init = 5.0e-03 --4.0e-8
chl_influx = 10.0e-03

---------------------------------------------------------------------
-- functions steering tempo-spatial parameterization of simulation --
---------------------------------------------------------------------

-- firing pattern of the synapse
synSubset = 2
synStartTime = 0.0
chlEntryDuration = 0.01

-- burst of chloride influx for active synapses (~1200 ions)
freq = 50      -- spike train frequency (Hz) (the ineq. 1/freq > chlEntryDuration must hold)
nSpikes = 1   -- number of spikes	
function ourNeumannBndCA(x, y, z, t, si)    
    if (synStartTime < t and t <= synStartTime+chlEntryDuration)
    then efflux = 0.0 -- -2e6
    else efflux = 0.0
    end
    
    return -efflux
end

-------------------------------
-- setup approximation space --
-------------------------------
-- create, load, refine and distribute domain
print("create, refine and distribute domain")

dom = util.CreateDomain(gridName, 0, neededSubsets)

balancer.partitioner = "parmetis"
balancer.staticProcHierarchy = true
balancer.redistSteps		= 0
balancer.firstDistLvl = -1

balancer.ParseParameters()
balancer.PrintParameters()

write(">> distributing and refining domain ...\n")
-- in parallel environments: use a load balancer to distribute the grid
loadBalancer = balancer.CreateLoadBalancer(dom)
balancer.RefineAndRebalanceDomain(dom, numRefs, loadBalancer)
write(">> distributing done\n")

print(dom:domain_info():to_string())

if load_balancer ~= nil then
	loadBalancer:print_quality_records()
end


--[[
--print("Saving domain grid and hierarchy.")
--SaveDomain(dom, "refined_grid_p" .. ProcRank() .. ".ugx")
--SaveGridHierarchyTransformed(dom:grid(), "refined_grid_hierarchy_p" .. ProcRank() .. ".ugx", 2.0)
print("Saving parallel grid layout")
SaveParallelGridLayout(dom:grid(), "parallel_grid_layout_p"..ProcRank()..".ugx", 2.0)
--]]

-- create approximation space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)

--cytVol = "cyt"
plMem = "mem_cyt, syn"
plMem_vec = {"mem_cyt", "syn"}

-- collect several subset names in subdomain variables
measZones = "measZone"..1
for i=2,201 do
	measZones = measZones .. ", measZone" .. i
end
--cytVol = cytVol .. ", " .. measZones
cytVol = measZones

outerDomain = cytVol .. ", " .. plMem

approxSpace:add_fct("chl_cyt", "Lagrange", 1, outerDomain)

approxSpace:init_levels()
approxSpace:print_layout_statistic()
approxSpace:print_statistic()

OrderCuthillMcKee(approxSpace, true);

--------------------------------------------------
-- setup FV convection-diffusion element discretization --
----------------------------------------------------------
print ("Setting up Assembling")

-- Note: No VelocityField and Reaction is set. The assembling assumes default
--       zero values for them
elemDiscCYT = ConvectionDiffusion("chl_cyt", cytVol, "fv1")
elemDiscCYT:set_diffusion(D_chl)
elemDiscCYT:set_upwind(NoUpwind())

------------------------------
-- setup Neumann boundaries --
------------------------------
-- synaptic activity
neumannDiscCA = UserFluxBoundaryFV1("chl_cyt", "syn")
neumannDiscCA:set_flux_function("ourNeumannBndCA")

------------------------------------------
-- setup complete domain discretization --
------------------------------------------
domainDisc = DomainDiscretization(approxSpace)

-- diffusion discretizations
domainDisc:add(elemDiscCYT)

-- (outer) boundary conditions
domainDisc:add(neumannDiscCA)

assTuner = domainDisc:ass_tuner()

-------------------------------
-- setup time discretization --
-------------------------------
timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(1.0) -- 1.0 is implicit Euler

-- create operator from discretization
op = AssembledLinearOperator()
op:set_discretization(timeDisc)

------------------
-- solver setup --
------------------

convCheck = CompositeConvCheck(approxSpace, 20, 1e-16, 1e-08)
convCheck:set_verbose(true)
convCheck:set_time_measurement(true)

cgSolver = CG()
cgSolver:set_preconditioner(ilu)
cgSolver:set_convergence_check(convCheck)

-------------
-- solving --
-------------

-- get grid function
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)

u:set(0.0)
b:set(0.0)

-- set initial value
Interpolate(chl_influx, u, "chl_cyt", "measZone101" , 0.0)
Interpolate(chl_cyt_init, u, "chl_cyt", 0.0)
Interpolate(chl_influx, u, "chl_cyt", "measZone101" , 0.0)


-- timestep in seconds
dt = timeStepStart
time = 0.0
step = 0

if (generateVTKoutput) then
	out = VTKOutput()
	out:print(fileName .. "vtk/result", u, step, time)
end

-- taking an initial measurement of all unknwons in all measurement zones on the ER membrane
-- the folder "meas" must exist in your file output directory specified in fileName
takeMeasurement(u, time, measZones, "chl_cyt", fileName .. "meas/data")

-- create new grid function for old value
uOld = u:clone()

-- store grid function in vector of  old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)


dtChanged = true
min_dt = timeStep / math.pow(2,15)
cb_interval = 10
lv = startLv
levelUpDelay = chlEntryDuration;
cb_counter = {}
for i=0,startLv do cb_counter[i]=0 end
while endTime-time > 0.001*dt do
	-- setup time Disc for old solutions and timestep
	timeDisc:prepare_step(solTimeSeries, dt)
	
	-- assemble linear problem
	matrixIsConst = dtChanged == false
	assTuner:set_matrix_is_const(matrixIsConst)
	if AssembleLinearOperatorRhsAndSolution(op, u, b) == false then 
		print("Could not assemble operator"); exit(); 
	end
	
	-- apply linear solver
	ilu:set_disable_preprocessing(matrixIsConst)
	if ApplyLinearSolver(op, u, b, cgSolver) == false then
		dtChanged = true
		
		curr_dt = curr_dt/dtred
		lv = lv + 1
		VecScaleAssign(u, 1.0, solTimeSeries:latest())

		-- halve time step and try again unless time step below minimum
		if curr_dt < min_dt
		then 
			print ("Time step below minimum. Aborting. Failed at point in time " .. time .. ".")
			time = endTime
		else
			print ("Retrying with half the time step...")
			cb_counter[lv] = 0
		end
	else
		dtChanged = false
		
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
			dtChanged = false
		end
		
		-- plot solution every plotStep seconds
		if (generateVTKoutput) then
			if math.abs(time/plotStep - math.floor(time/plotStep+0.5)) < 1e-5 then
				out:print(fileName .. "vtk/result", u, math.floor(time/plotStep+0.5), time)
			end
		end
	
		-- take measurements in measurement zones
		takeMeasurement(u, time, measZones, "chl_cyt", fileName .. "meas/data")	
	
		-- get oldest solution
		oldestSol = solTimeSeries:oldest()
		
		-- copy values into oldest solution (we reuse the memory here)
		VecScaleAssign(oldestSol, 1.0, u)
		
		-- push oldest solutions with new values to front, oldest sol pointer is popped from end
		solTimeSeries:push_discard_oldest(oldestSol, time)
		
		print("++++++ POINT IN TIME  " .. math.floor(time/dt+0.5)*dt .. "s  END ++++++++");
	end
end

-- output of load balancing quality statistics
if loadBalancer ~= nil then
	loadBalancer:print_quality_records()
end

--out_error:write_time_pvd(fileName .. "vtk/error_estimator", u_vtk)

-- end timeseries, produce gathering file
if (generateVTKoutput) then out:write_time_pvd(fileName .. "vtk/result", u) end
