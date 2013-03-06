----------------------------------------------------------------
--  Script for test purposes								  --
--                                                            --
--  Author: Markus Breit                                      --
----------------------------------------------------------------

-- load pre-implemented lua functions
ug_load_script("ug_util.lua")

-- dimension
dim = 3

-- choose dimension and algebra
InitUG(dim, AlgebraType("CPU", 1));

-- choice of grid
gridName = "simple_reticulum_3d.ugx"				-- for testing

-- refinements before distributing grid
numPreRefs = util.GetParamNumber("-numPreRefs", 0)

-- total refinements
numRefs = util.GetParamNumber("-numRefs",    0)

-- choose length of maximal time step during the whole simulation
timeStep = 0.001

-- choose length of time step at the beginning
-- if not timeStepStart = 2^(-n)*timeStep, take nearest lower number of that form
timeStepStart = timeStep

-- choose end time
nTimeSteps = 1
endTime = nTimeSteps*timeStep


---------------
-- constants --
---------------
-- diffusion coefficients
D_cac = 220.0
D_cae = 220.0

-- initial concentrations
ca_cyt_init = 4.0e-8
ca_er_init = 2.5e-4

---------------------------------------------------------------------
-- functions steering tempo-spatial parameterization of simulation --
---------------------------------------------------------------------
function CaCytStart(x, y, z, t)
    return ca_cyt_init
end

function CaERStart(x, y, z, t)
    return ca_er_init
end

function ourDiffTensorCAcyt(x, y, z, t)
    return	D_cac, 0, 0,
            0, D_cac, 0,
            0, 0, D_cac
end

function ourDiffTensorCAer(x, y, z, t)
    return	D_cae, 0, 0,
            0, D_cae, 0,
            0, 0, D_cae
end

function ourRhs(x, y, z, t)
    return 0;
end


-- burst of calcium influx for active synapses (~1200 ions)
function ourNeumannBndCA(x, y, z, t, si)
	if 	(si==6 and t<=0.01)
	then efflux = -2e-4
	else efflux = 0.0
	end
    return true, efflux
end


-------------------------------
-- setup approximation space --
-------------------------------
-- create, load, refine and distribute domain
print("create, refine and distribute domain")
neededSubsets = {}
distributionMethod = "metis"
weightingFct = InterSubsetPartitionWeighting()
weightingFct:set_default_weights(1,1)
weightingFct:set_inter_subset_weight(0, 1, 1000)
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets, distributionMethod, nil, nil, nil, weightingFct)

--[[
--print("Saving domain grid and hierarchy.")
--SaveDomain(dom, "refined_grid_p" .. GetProcessRank() .. ".ugx")
--SaveGridHierarchyTransformed(dom:grid(), "refined_grid_hierarchy_p" .. GetProcessRank() .. ".ugx", 20.0)
print("Saving parallel grid layout")
SaveParallelGridLayout(dom:grid(), "parallel_grid_layout_p"..GetProcessRank()..".ugx", 20.0)
--]]

-- create approximation space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)

cytVol = "cyt"

nucVol = "nuc"
nucMem = "mem_nuc"

erVol = "er"

plMem = "mem_cyt"
synapses = "syn1"
plMem = plMem ..", "..synapses

erMem = "mem_er"
measZonesERM = "measZoneERM"..1
erMem = erMem .. ", " .. measZonesERM

outerDomain = cytVol .. ", " .. nucVol .. ", " .. nucMem .. ", " .. plMem .. ", " .. erMem
innerDomain = erVol .. ", " .. erMem

approxSpace:add_fct("ca_er", "Lagrange", 1, innerDomain)
approxSpace:add_fct("ca_cyt", "Lagrange", 1, outerDomain)

approxSpace:init_levels()
approxSpace:print_layout_statistic()
approxSpace:print_statistic()

--------------------------
-- setup user functions --
--------------------------
print ("Setting up Assembling")

-- start value function setup
CaCytStartValue = LuaUserNumber3d("CaCytStart")
CaERStartValue = LuaUserNumber3d("CaERStart")

-- diffusion Tensor setup
diffusionMatrixCAcyt = LuaUserMatrix3d("ourDiffTensorCAcyt")
diffusionMatrixCAer = LuaUserMatrix3d("ourDiffTensorCAer")

-- rhs setup
rhs = LuaUserNumber3d("ourRhs")

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

elemDiscER = ConvectionDiffusion("ca_er", erVol, "fv1") 
elemDiscER:set_diffusion(diffusionMatrixCAer)
elemDiscER:set_source(rhs)
elemDiscER:set_upwind(upwind)

elemDiscCYT = ConvectionDiffusion("ca_cyt", cytVol..", "..nucVol, "fv1")
elemDiscCYT:set_diffusion(diffusionMatrixCAcyt)
elemDiscCYT:set_source(rhs)
elemDiscCYT:set_upwind(upwind)

------------------------------
-- setup Neumann boundaries --
------------------------------
-- synaptic activity
neumannDiscCA = NeumannBoundary("ca_cyt")
neumannDiscCA:add("ourNeumannBndCA", plMem, "cyt")

------------------------------------------
-- setup complete domain discretization --
------------------------------------------
domainDisc = DomainDiscretization(approxSpace)

-- diffusion discretizations
domainDisc:add(elemDiscER)
domainDisc:add(elemDiscCYT)

-- (outer) boundary conditions
domainDisc:add(neumannDiscCA)

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

base = LinearSolver()
base:set_convergence_check(baseConvCheck)
base:set_preconditioner(gs)

gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(timeDisc)
gmg:set_base_level(0)
--gmg:set_parallel_base_solver(false)
gmg:set_base_solver(base)
gmg:set_smoother(gs)
gmg:set_cycle_type(1)
gmg:set_num_presmooth(3)
gmg:set_num_postsmooth(3)

-- biCGstab --
convCheck = ConvCheck()
convCheck:set_maximum_steps(2000)
convCheck:set_minimum_defect(1e-24)
convCheck:set_reduction(1e-06)
convCheck:set_verbose(true)
bicgstabSolver = BiCGStab()
bicgstabSolver:set_preconditioner(ilu)
bicgstabSolver:set_convergence_check(convCheck)

-----------------------
-- non linear solver --
-----------------------
-- convergence check
newtonConvCheck = CompositeConvCheck3dCPU1(approxSpace)
newtonConvCheck:set_functions("")
newtonConvCheck:set_maximum_steps(10)
newtonConvCheck:set_minimum_defect({}, 1e-18)
newtonConvCheck:set_reduction({}, 1e-08)
newtonConvCheck:set_verbose(true)
newtonConvCheck:timeMeasurement(true)

-- Newton solver
newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(bicgstabSolver)
newtonSolver:set_convergence_check(newtonConvCheck)

newtonSolver:init(op)

-------------
-- solving --
-------------

-- get grid function
u = GridFunction(approxSpace)

-- set initial value
Interpolate(CaCytStartValue, u, "ca_cyt", 0.0)
Interpolate(CaERStartValue, u, "ca_er", 0.0)

-- timestep in seconds
dt = timeStepStart
time = 0.0
step = 0

-- create new grid function for old value
uOld = u:clone()

-- store grid function in vector of  old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)



while endTime-time > 0.001*dt do
	print("++++++ POINT IN TIME  " .. math.floor((time+dt)/dt+0.5)*dt .. "s  BEGIN ++++++")
	
	-- setup time Disc for old solutions and timestep
	timeDisc:prepare_step(solTimeSeries, dt)
	
	-- prepare newton solver
	if newtonSolver:prepare(u) == false then print ("Newton solver failed at step "..step.."."); exit(); end 
	
	-- apply newton solver
	if newtonSolver:apply(u) == false
	then
		time = endTime
	else
		-- update new time
		time = solTimeSeries:time(0) + dt
		
		-- get oldest solution
		oldestSol = solTimeSeries:oldest()
		
		-- copy values into oldest solution (we reuse the memory here)
		VecScaleAssign(oldestSol, 1.0, u)
		
		-- push oldest solutions with new values to front, oldest sol pointer is popped from end
		solTimeSeries:push_discard_oldest(oldestSol, time)
		
		print("++++++ POINT IN TIME  " .. math.floor(time/dt+0.5)*dt .. "s  END ++++++++");
	end

end
