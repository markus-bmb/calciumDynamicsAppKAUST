--------------------------------------------------------------------------------
-- This script produces a fluorescence movie from data files created using    --
-- the script "data_import_from_anna.lua".                                    --
--                                                                            --
-- Author: Markus Breit                                                       --
-- Date:   2014-09-16	                                                      --
--------------------------------------------------------------------------------

-- load pre-implemented lua functions
ug_load_script("ug_util.lua")


-- choose dimension and algebra
InitUG(3, AlgebraType("CPU", 1))

-- choice of grid
gridName = util.GetParam("-grid", "calciumDynamics_app/grids/camh36.ugx")

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


-------------------------------
-- setup approximation space --
-------------------------------
-- create, load, refine and distribute domain
dom = util.CreateDomain(gridName, 0, {"cyt", "er", "mem_cyt", "mem_er"})

-- create approximation space
all_subsets = "cyt, er, mem_cyt, mem_er"
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("dfof", "Lagrange", 1, all_subsets)


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
	file = string.format("/Users/mbreit/GCSC/Calcium_data_Anna/flourescence_data_camh36/dFoF_%.4f", time)
	
	-- import values from file
	importSolution(u, all_subsets, "dfof", file)
	
	-- plot
	if math.abs((time-startTime)/plotStep - math.floor((time-startTime)/plotStep+0.5)) < 1e-5 then
		out:print(fileName .. "vtk/fluorescence", u, math.floor((time-startTime)/plotStep+0.5), time)
	end

	time = time + dt
end

-- end timeseries, produce gathering file
out:write_time_pvd(fileName .. "vtk/fluorescence", u)



