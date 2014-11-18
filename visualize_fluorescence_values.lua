----------------------------------------------------------------
--  This script produces a fluorescence movie from data files --
--  created using the script "data_import_from_anna.lua" and. --
--															  --
--  Author: Markus Breit                                      --
----------------------------------------------------------------

-- load pre-implemented lua functions
ug_load_script("ug_util.lua")

-- dimension
dim = 3

-- choose dimension and algebra
InitUG(dim, AlgebraType("CPU", 1));

-- choice of grid
gridName = "camh36.ugx"

-- choose length of maximal time step during the whole simulation
timeStep = util.GetParamNumber("-tstep", 0.01)
	
-- choose start and end time
startTime = util.GetParamNumber("-startTime", 0.0)
nTimeSteps = util.GetParamNumber("-nTimeSteps", 1)
def_endTime = nTimeSteps*timeStep
endTime = util.GetParamNumber("-endTime", def_endTime)

-- chose plotting interval
plotStep = util.GetParamNumber("-pstep", timeStep)

-- choose outfile directory
fileName = util.GetParam("-outName", "fluorescence_data")
fileName = fileName.."/"


---------------
-- constants --
---------------

-------------------------------
-- setup approximation space --
-------------------------------
-- create, load, refine and distribute domain
print("create, refine and distribute domain")
neededSubsets = {}
distributionMethod = "bisection"
dom = util.CreateAndDistributeDomain(gridName, 0, 0, neededSubsets, distributionMethod)

-- create approximation space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)

all_subsets = "cyt, er, mem_cyt, mem_er"

approxSpace:add_fct("dfof", "Lagrange", 1, all_subsets)

approxSpace:init_levels()
approxSpace:print_layout_statistic()
approxSpace:print_statistic()

------------------------------------------
-- setup complete domain discretization --
------------------------------------------
domainDisc = DomainDiscretization(approxSpace)

-------------------------------
-- setup time discretization --
-------------------------------
timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(1.0) -- 1.0 is implicit Euler

-------------------
-- time stepping --
-------------------

-- create operator from discretization
op = AssembledOperator()
op:set_discretization(timeDisc)
op:init()

-- get grid function
u = GridFunction(approxSpace)
u:enforce_consistent_type()

-- vtk output interface
out = VTKOutput()


-- time loop
time = startTime
dt = timeStep
while endTime-time > 0.001*dt do
	-- construct file name
	file = string.format("/Users/mbreit/GCSC/Calcium_data_Anna/flourescence_data_camh36/dFoF_%.4f.dat", time)
	
	-- import values from file
	import_solution(u, all_subsets, "dfof", file)
	
	-- plot
	if math.abs((time-startTime)/plotStep - math.floor((time-startTime)/plotStep+0.5)) < 1e-5 then
		out:print(fileName .. "vtk/fluorescence", u, math.floor((time-startTime)/plotStep+0.5), time)
	end

	time = time + dt
end

-- end timeseries, produce gathering file
out:write_time_pvd(fileName .. "vtk/fluorescence", u)



