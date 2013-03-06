----------------------------------------------------------
--
--   Example - Script to test simulation on 3d reconstructed neuron
--
--   Author: Markus Breit
--
----------------------------------------------------------

ug_load_script("ug_util.lua")

-- Dimension
dim = 3

-- choose dimension and algebra
InitUG(dim, AlgebraType("CPU", 1));

-- Grid
--gridName = "rc19_amp.ugx"
--gridName = "RC19amp_ug4_finished.ugx"
gridName = "simple_reticulum_3d.ugx"

-- Total refinements
numMaxRefs    = util.GetParamNumber("-numRefs",    0)

-- choose number of time steps
nTimeSteps =  util.GetParamNumber("-nTimeSteps", 1)

-- choose length of time step
timeStep =  util.GetParamNumber("-tstep", 0.001)


--------------------------------
-- constants
--------------------------------
-- total cytosolic calmodulin concentration
-- (double the real value in order to simulate two binding sites in one)
totalClmd = 2*21.9e-6;

-- calmodulin diffusion coefficient
D_clmd = 0.25;

-- calmodulin binding rates
k_bind_clmd_c = 	2.3e06
k_unbind_clmd_c = 	2.4
k_bind_clmd_n = 	1.6e08
k_unbind_clmd_n = 	405.0

-- initial concentrations
ca_cyt_init = 7.5e-8
ca_er_init = 2.5e-4
ip3_init = 4.0e-8
clmd_c_init = totalClmd / (k_bind_clmd_c/k_unbind_clmd_c*ca_cyt_init + 1)
clmd_n_init = totalClmd / (k_bind_clmd_n/k_unbind_clmd_n*ca_cyt_init + 1)


--------------------------------
-- User Data Functions (begin)
--------------------------------
function CaCytStart(x, y, z, t)
    return ca_cyt_init
end

function CaERStart(x, y, z, t)
    return ca_er_init
end

function IP3Start(x, y, z, t)
    return ip3_init --1.0e-6
end

function ClmdCStart(x, y, z, t)
    return clmd_c_init
end

function ClmdNStart(x, y, z, t)
    return clmd_n_init
end

function ourDiffTensorCAcyt(x, y, z, t)
    return	220, 0, 0,
            0, 220, 0,
            0, 0, 220
end

function ourDiffTensorCAer(x, y, z, t)
    return	220, 0, 0,
            0, 220, 0,
            0, 0, 220
end

function ourDiffTensorIP3(x, y, z, t)
    return	280, 0, 0,
            0, 280, 0,
            0, 0, 280
end

function ourDiffTensorClmd(x, y, z, t)
    return	D_clmd, 0, 0,
            0, D_clmd, 0,
            0, 0, D_clmd
end

function ourRhs(x, y, z, t)
    return 0;
end


-- firing pattern
syns = {}
for i=6,15 do
	syns["start"..i] = 0.005*(i-6)
	syns["end"..i] = 0.005*(i-6)+0.01
end

function ourNeumannBndCA(x, y, z, t, si)
	-- burst for active synapses
	if 	(si>=6 and si<=15 and syns["start"..si]<=t and t<syns["end"..si])
	--then efflux = -5e-6 * 11.0/16.0*(1.0+5.0/((10.0*(t-syns["start"..si])+1)*(10.0*(t-syns["start"..si])+1)))
	then efflux = -5e-4
	else efflux = 0.0
	end
	--]]
	
    return true, efflux
end

ip3EntryDelay = 0.005
ip3EntryDuration = 1.0
function ourNeumannBndIP3(x, y, z, t, si)
	-- burst for active synapses
	if 	(si>=6 and si<=15 and syns["start"..si]+ip3EntryDelay<=t and t<syns["start"..si]+ip3EntryDelay+ip3EntryDuration)
	then efflux = -2.1e-5
	else efflux = 0.0
	end
	
    return true, efflux
end




-- iteration over number of refinements
dSol = GridFunction(approxSpace)
for numRefs = 0, numMaxRefs do


--------------------------------
-- User Data Functions (end)
--------------------------------

-- create, load, refine and distribute domain
print("create, refine and distribute domain")
neededSubsets = {}
dom = util.CreateAndDistributeDomain(gridName, numRefs, 0, neededSubsets)

-- write grid to file for test purpose
SaveDomain(dom, "refined_grid.ugx")

-- create approximation space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)
innerDomain = "er, mem_er"
--outerDomain = "cyt, mem_cyt, mem_er"
outerDomain = "cyt, nuc, mem_cyt, mem_er, mem_nuc"
synapses = ""
---[[
for i=1,1 do
	synapses = synapses .. ", syn" .. i
end
outerDomain = outerDomain .. synapses
--]]
approxSpace:add_fct("ca_er", "Lagrange", 1, innerDomain)
approxSpace:add_fct("ca_cyt", "Lagrange", 1, outerDomain)
approxSpace:add_fct("ip3", "Lagrange", 1, outerDomain)
approxSpace:add_fct("clmd_c", "Lagrange", 1, outerDomain)
approxSpace:add_fct("clmd_n", "Lagrange", 1, outerDomain)

--OrderLex(approxSpace, "lr")
--OrderCuthillMcKee(approxSpace, true);
approxSpace:init_levels()
approxSpace:print_layout_statistic()
approxSpace:print_statistic()

-------------------------------------------
--  Setup User Functions
-------------------------------------------
print ("Setting up Assembling")

-- Start value function setup
    CaCytStartValue = LuaUserNumber3d("CaCytStart")
    CaERStartValue = LuaUserNumber3d("CaERStart")
    IP3StartValue = LuaUserNumber3d("IP3Start")
    ClmdCStartValue = LuaUserNumber3d("ClmdCStart")
    ClmdNStartValue = LuaUserNumber3d("ClmdNStart")

-- Diffusion Tensor setup
	diffusionMatrixCAcyt = LuaUserMatrix3d("ourDiffTensorCAcyt")
	diffusionMatrixCAer = LuaUserMatrix3d("ourDiffTensorCAer")
	diffusionMatrixIP3 = LuaUserMatrix3d("ourDiffTensorIP3")
	diffusionMatrixClmd = LuaUserMatrix3d("ourDiffTensorClmd")

-- rhs setup
	rhs = LuaUserNumber3d("ourRhs")
	--rhs = ConstUserNumber(0.0)

-- neumann setup
	neumannCA = LuaCondUserNumber3d("ourNeumannBndCA")
	--neumannCA = ConstUserNumber(0.0)
	neumannIP3 = LuaCondUserNumber3d("ourNeumannBndIP3")

--[[
-- dirichlet setup
	dirichlet = LuaCondUserNumber3d("ourDirichletBnd")
	
-- dirichlet setup
	membraneDirichlet = LuaCondUserNumber3d("membraneDirichletBnd")
--]]

-----------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
-----------------------------------------------------------------

-- Note: No VelocityField and Reaction is set. The assembling assumes default
--       zero values for them

if dim == 2 then 
    upwind = NoUpwind2d()
elseif dim == 3 then 
    upwind = NoUpwind3d()
end

elemDiscER = ConvectionDiffusion("ca_er", "er", "fv1") 
elemDiscER:set_diffusion(diffusionMatrixCAer)
elemDiscER:set_source(rhs)
elemDiscER:set_upwind(upwind)

elemDiscCYT = ConvectionDiffusion("ca_cyt", "cyt, nuc", "fv1")
elemDiscCYT:set_diffusion(diffusionMatrixCAcyt)
elemDiscCYT:set_source(rhs)
elemDiscCYT:set_upwind(upwind)

elemDiscIP3 = ConvectionDiffusion("ip3", "cyt, nuc", "fv1")
elemDiscIP3:set_diffusion(diffusionMatrixIP3)
elemDiscIP3:set_source(rhs)
elemDiscIP3:set_upwind(upwind)

elemDiscClmdC = ConvectionDiffusion("clmd_c", "cyt, nuc", "fv1")
elemDiscClmdC:set_diffusion(diffusionMatrixClmd)
elemDiscClmdC:set_source(rhs)
elemDiscClmdC:set_upwind(upwind)

elemDiscClmdN = ConvectionDiffusion("clmd_n", "cyt, nuc", "fv1")
elemDiscClmdN:set_diffusion(diffusionMatrixClmd)
elemDiscClmdN:set_source(rhs)
elemDiscClmdN:set_upwind(upwind)

-----------------------------------------------------------------
--  Setup reaction terms of buffering
-----------------------------------------------------------------
elemDiscBuffering = FV1Buffer("cyt")	-- where buffering occurs
elemDiscBuffering:add_reaction(
	"clmd_c",						-- the buffering substance
	"ca_cyt",						-- the buffered substance
	totalClmd,						-- total amount of buffer
	k_bind_clmd_c,					-- binding rate constant
	k_unbind_clmd_c)				-- unbinding rate constant
elemDiscBuffering:add_reaction(
	"clmd_n",
	"ca_cyt",				
	totalClmd,
	k_bind_clmd_n,
	k_unbind_clmd_n)
-----------------------------------------------------------------
--  Setup inner boundary (channels on ER membrane)
-----------------------------------------------------------------

-- We pass the function needed to evaluate the flux function here.
-- The order, in which the discrete fcts are passed, is crucial!
innerDisc = FV1InnerBoundaryCalciumER("ca_cyt, ca_er, ip3", "mem_er")

-----------------------------------------------------------------
--  Setup Neumann Boundary
-----------------------------------------------------------------

neumannDiscCA = NeumannBoundary("ca_cyt")
neumannDiscCA:add(neumannCA, "mem_cyt" .. synapses, "cyt")
neumannDiscIP3 = NeumannBoundary("ip3")
neumannDiscIP3:add(neumannIP3, "mem_cyt" .. synapses, "cyt")


--[[
-----------------------------------------------------------------
--  Setup Dirichlet Boundary
-----------------------------------------------------------------

dirichletBND = DirichletBoundary()
dirichletBND:add(dirichlet, "c", "Boundary, MembraneBnd")

membraneDirichletBND = DirichletBoundary()
membraneDirichletBND:add(membraneDirichlet, "c_membrane", "MembraneBnd")

--]]

-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization(approxSpace)

-- diffusion disc.s
domainDisc:add(elemDiscER)
domainDisc:add(elemDiscCYT)
domainDisc:add(elemDiscIP3)
domainDisc:add(elemDiscClmdC)
domainDisc:add(elemDiscClmdN)

-- buffering disc
domainDisc:add(elemDiscBuffering)

-- (outer) boundary conditions
domainDisc:add(neumannDiscCA)
domainDisc:add(neumannDiscIP3)

-- ER flux
domainDisc:add(innerDisc)
--domainDisc:add(dirichletBND)
--domainDisc:add(membraneDirichletBND)

-------------------------------------------
--  Setup Time Discretization
-------------------------------------------

-- create time discretization
timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(1.0) -- 1.0 is implicit euler

-- create operator from discretization
op = AssembledOperator()
op:set_discretization(timeDisc)
op:init()

-------------------------------------------
--  Algebra
-------------------------------------------
print ("Setting up Algebra Solver")

-- create algebraic Preconditioner
jac = Jacobi()
jac:set_damp(0.8)
gs = GaussSeidel()
sgs = SymmetricGaussSeidel()
bgs = BackwardGaussSeidel()
ilu = ILU()
ilut = ILUT()

-- exact Solver
exactSolver = LU()


-- create GMG ---
-----------------

	-- Base Solver
	baseConvCheck = ConvCheck()
	baseConvCheck:set_maximum_steps(10000)
	baseConvCheck:set_minimum_defect(1e-28)
	baseConvCheck:set_reduction(1e-2)
	baseConvCheck:set_verbose(false)
	--base = LU()
	---[[
	base = LinearSolver()
	base:set_convergence_check(baseConvCheck)
	base:set_preconditioner(gs)
	--]]
	-- Gemoetric Multi Grid
	gmg = GeometricMultiGrid(approxSpace)
	gmg:set_discretization(timeDisc)
	--gmg:set_surface_level(numRefs)
	gmg:set_base_level(0)
	gmg:set_base_solver(base)
	gmg:set_smoother(gs)
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(3)
	gmg:set_num_postsmooth(3)


-- create Convergence Check
convCheck = ConvCheck()
convCheck:set_maximum_steps(50)
convCheck:set_minimum_defect(1e-24)
convCheck:set_reduction(1e-6)
convCheck:set_verbose(true)

--[[
-- create Linear Solver
linSolver = LinearSolver()
linSolver:set_preconditioner(gmg)
linSolver:set_convergence_check(convCheck)
--]]

-- create BiCGStab Solver
bicgstabSolver = BiCGStab()
bicgstabSolver:set_preconditioner(gmg)
bicgstabSolver:set_convergence_check(convCheck)

-------------------------------------------
--  Non linear Solver
-------------------------------------------

-- convergence check
---[[
newtonConvCheck = CompositeConvCheck3dCPU1(approxSpace)
newtonConvCheck:set_functions("ip3")
newtonConvCheck:set_maximum_steps(20)
newtonConvCheck:set_minimum_defect({1e-20}, 1e-20)
newtonConvCheck:set_reduction({1e-08}, 1e-08)
newtonConvCheck:set_verbose(true)
newtonConvCheck:timeMeasurement(true)
--]]
--[[
newtonConvCheck = ConvCheck()
newtonConvCheck:set_maximum_steps(20)
newtonConvCheck:set_minimum_defect(1e-21)
newtonConvCheck:set_reduction(1e-08)
newtonConvCheck:set_verbose(true)
--]]
newtonLineSearch = StandardLineSearch()
newtonLineSearch:set_maximum_steps(8)
newtonLineSearch:set_lambda_start(1)
newtonLineSearch:set_accept_best(true)

-- create Newton Solver
newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(bicgstabSolver)
newtonSolver:set_convergence_check(newtonConvCheck)
--newtonSolver:set_line_search(newtonLineSearch)

newtonSolver:init(op)

-------------------------------------------
--  Solving
-------------------------------------------

-- get grid function
u = GridFunction(approxSpace)

-- set initial value
Interpolate(CaCytStartValue, u, "ca_cyt", 0.0)
Interpolate(CaERStartValue, u, "ca_er", 0.0)
Interpolate(IP3StartValue, u, "ip3", 0.0)
Interpolate(ClmdCStartValue, u, "clmd_c", 0.0)
Interpolate(ClmdNStartValue, u, "clmd_n", 0.0)

-- timestep in seconds
dt = timeStep
time = 0.0
step = 0

-- filename
filename = "/Users/markus/Developing/ug4/trunk/bin/retic/result"

-- write start solution
print("Writing start values")
out = VTKOutput()
out:print(filename, u, step, time)
--takeMeasurement(u, approxSpace, time, "nuc", "ca_cyt", "measurements")
--exportSolution(u, approxSpace, time, "mem_cyt", "ca_cyt", "solution/solution");

-- some info output
print( "   numPreRefs is   " .. numPreRefs ..     ",  numRefs is         " .. numRefs)
print( "   nTimeSteps is " .. nTimeSteps)

-- create new grid function for old value
uOld = u:clone()

-- store grid function in vector of  old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)

for step = 1, nTimeSteps do
	print("++++++ TIMESTEP " .. step .. " BEGIN ++++++")
	
	-- choose time step (currently constant)
	do_dt = dt
	
	-- setup time Disc for old solutions and timestep
	timeDisc:prepare_step(solTimeSeries, do_dt)
	
	-- prepare newton solver
	if newtonSolver:prepare(u) == false then print ("Newton solver failed at step "..step.."."); exit(); end 
	
	-- apply newton solver
	if newtonSolver:apply(u) == false then print ("Newton solver failed at step "..step.."."); exit(); end 
	
	-- update new time
	time = solTimeSeries:time(0) + do_dt
	
	-- plot solution
	print("Time: "..time)
	if math.abs(time/0.001 - math.floor(time/0.001+0.5)) < 1e-9
	then out:print(filename, u, math.floor(time/0.001+0.5), time)
	end
	
	-- take measurement in nucleus
	--takeMeasurement(u, approxSpace, time, "nuc", "ca_cyt", "measurements")
	
	-- export solution of ca on mem_er
	--exportSolution(u, approxSpace, time, "mem_cyt", "ca_cyt", "solution/solution");
	
	
	-- get oldest solution
	oldestSol = solTimeSeries:oldest()
	
	-- copy values into oldest solution (we reuse the memory here)
	VecScaleAssign(oldestSol, 1.0, u)
	
	-- push oldest solutions with new values to front, oldest sol pointer is poped from end
	solTimeSeries:push_discard_oldest(oldestSol, time)
	
	print("++++++ TIMESTEP " .. step .. "  END ++++++");
end

-- end timeseries, produce gathering file
out:write_time_pvd(filename, u)
