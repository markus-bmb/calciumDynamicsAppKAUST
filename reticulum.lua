----------------------------------------------------------------
--  Example script for simulation on 3d reconstructed neuron  --
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
--gridName = "rc19_amp.ugx"
gridName = "RC19amp_ug4_finished.ugx" -- defect!?
--gridName = "simple_reticulum_3d.ugx"

-- refinements before distributing grid
numPreRefs = util.GetParamNumber("-numPreRefs", 0)

-- total refinements
numRefs    = util.GetParamNumber("-numRefs",    0)

-- choose number of time steps
nTimeSteps =  util.GetParamNumber("-nTimeSteps", 5)

-- choose length of time step
timeStep =  util.GetParamNumber("-tstep", 0.001)

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
ca_cyt_init = 7.5e-8
ca_er_init = 2.5e-4
ip3_init = 4.0e-8
clb_init = totalClb / (k_bind_clb/k_unbind_clb*ca_cyt_init + 1)


-- calmodulin --
totalClmd = 2*21.9e-6;

-- calmodulin diffusion coefficient
D_clm = 0.25;

-- calmodulin binding rates
k_bind_clmd_c = 	2.3e06
k_unbind_clmd_c = 	2.4
k_bind_clmd_n = 	1.6e08
k_unbind_clmd_n = 	405.0

-- initial concentrations
clmC_init = totalClmd / (k_bind_clmd_c/k_unbind_clmd_c*ca_cyt_init + 1)
clmN_init = totalClmd / (k_bind_clmd_n/k_unbind_clmd_n*ca_cyt_init + 1)


-- reaction reate IP3
reactionRateIP3 = 0.11

-- equilibrium concentration IP3
equilibriumIP3 = 4.0e-08

-- reation term IP3
reactionTermIP3 = -reactionRateIP3 * equilibriumIP3

---------------------------------------------------------------------
-- functions steering tempo-spatial parameterization of simulation --
---------------------------------------------------------------------
function CaCytStart(x, y, z, t)
    return ca_cyt_init
end

function CaERStart(x, y, z, t)
    return ca_er_init
end

function IP3Start(x, y, z, t)
    return ip3_init
end

function clbStart(x, y, z, t)
    return clb_init
end

function clmCStart(x, y, z, t)
    return clmC_init
end

function clmNStart(x, y, z, t)
    return clmN_init
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

function ourDiffTensorIP3(x, y, z, t)
    return	D_ip3, 0, 0,
            0, D_ip3, 0,
            0, 0, D_ip3
end

function ourDiffTensorClb(x, y, z, t)
    return	D_clb, 0, 0,
            0, D_clb, 0,
            0, 0, D_clb
end

function ourDiffTensorClm(x, y, z, t)
    return D_clm, 0, 0,
           0, D_clm, 0,
           0, 0, D_clm
end

function ourRhs(x, y, z, t)
    return 0;
end


-- firing pattern of the synapses
syns = {}
for i=6,13 do
	syns["start"..i] = 0.005*(i-6)
	syns["end"..i] = 0.005*(i-6)+0.01
end

-- burst of calcium influx for active synapses (~1200 ions)
function ourNeumannBndCA(x, y, z, t, si)
	if 	(si>=6 and si<=13 and syns["start"..si]<=t and t<syns["end"..si])
	--then efflux = -5e-6 * 11.0/16.0*(1.0+5.0/((10.0*(t-syns["start"..si])+1)*(10.0*(t-syns["start"..si])+1)))
	then efflux = -2e-4
	else efflux = 0.0
	end	
    return true, efflux
end

-- burst of ip3 at active synapses (triangular, immediate)
ip3EntryDelay = 0.000
ip3EntryDuration = 2.0
function ourNeumannBndIP3(x, y, z, t, si)
	if 	(si>=6 and si<=13 and syns["start"..si]+ip3EntryDelay<=t
	     and t<syns["start"..si]+ip3EntryDelay+ip3EntryDuration)
	then efflux = - 2.1e-6/1.188 * (1.0 - (t-syns["start"..si])/ip3EntryDuration)
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
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets)


-- create approximation space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)
innerDomain = "er, mem_er"
outerDomain = "cyt, nuc, mem_cyt, mem_er, mem_nuc"
synapses = ""
for i=1,8 do
	synapses = synapses .. ", syn" .. i
end
outerDomain = outerDomain .. synapses

approxSpace:add_fct("ca_er", "Lagrange", 1, innerDomain)
approxSpace:add_fct("ca_cyt", "Lagrange", 1, outerDomain)
approxSpace:add_fct("ip3", "Lagrange", 1, outerDomain)
approxSpace:add_fct("clb", "Lagrange", 1, outerDomain)
--approxSpace:add_fct("clm_c", "Lagrange", 1, outerDomain)
--approxSpace:add_fct("clm_n", "Lagrange", 1, outerDomain)

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
IP3StartValue = LuaUserNumber3d("IP3Start")
ClbStartValue = LuaUserNumber3d("clbStart")
ClmCStartValue = LuaUserNumber3d("clmCStart")
ClmNStartValue = LuaUserNumber3d("clmNStart")

-- diffusion Tensor setup
diffusionMatrixCAcyt = LuaUserMatrix3d("ourDiffTensorCAcyt")
diffusionMatrixCAer = LuaUserMatrix3d("ourDiffTensorCAer")
diffusionMatrixIP3 = LuaUserMatrix3d("ourDiffTensorIP3")
diffusionMatrixClb = LuaUserMatrix3d("ourDiffTensorClb")
diffusionMatrixClm = LuaUserMatrix3d("ourDiffTensorClm")

-- rhs setup
rhs = LuaUserNumber3d("ourRhs")

-- Neumann setup (Neumann-0 represented by declaring nothing)
neumannCA = LuaCondUserNumber3d("ourNeumannBndCA")
neumannIP3 = LuaCondUserNumber3d("ourNeumannBndIP3")

--[[
-- dirichlet setup
	dirichlet = LuaCondUserNumber3d("ourDirichletBnd")
	
-- dirichlet setup
	membraneDirichlet = LuaCondUserNumber3d("membraneDirichletBnd")
--]]

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

elemDiscER = ConvectionDiffusion("ca_er", "er") 
elemDiscER:set_disc_scheme("fv1")
elemDiscER:set_diffusion(diffusionMatrixCAer)
elemDiscER:set_source(rhs)
elemDiscER:set_upwind(upwind)

elemDiscCYT = ConvectionDiffusion("ca_cyt", "cyt, nuc")
elemDiscCYT:set_disc_scheme("fv1")
elemDiscCYT:set_diffusion(diffusionMatrixCAcyt)
elemDiscCYT:set_source(rhs)
elemDiscCYT:set_upwind(upwind)

elemDiscIP3 = ConvectionDiffusion("ip3", "cyt, nuc")
elemDiscIP3:set_disc_scheme("fv1")
elemDiscIP3:set_diffusion(diffusionMatrixIP3)
elemDiscIP3:set_reaction_rate(reactionRateIP3)
elemDiscIP3:set_reaction(reactionTermIP3)
elemDiscIP3:set_source(rhs)
elemDiscIP3:set_upwind(upwind)

elemDiscClb = ConvectionDiffusion("clb", "cyt, nuc")
elemDiscClb:set_disc_scheme("fv1")
elemDiscClb:set_diffusion(diffusionMatrixClb)
elemDiscClb:set_source(rhs)
elemDiscClb:set_upwind(upwind)

--[[
elemDiscClmC = ConvectionDiffusion("clm_c", "cyt, nuc")
elemDiscClb:set_disc_scheme("fv1")
elemDiscClb:set_diffusion(diffusionMatrixClm)
elemDiscClb:set_source(rhs)
elemDiscClb:set_upwind(upwind)

elemDiscClmN = ConvectionDiffusion("clm_n", "cyt, nuc")
elemDiscClb:set_disc_scheme("fv1")
elemDiscClb:set_diffusion(diffusionMatrixClm)
elemDiscClb:set_source(rhs)
elemDiscClb:set_upwind(upwind)
--]]
---------------------------------------
-- setup reaction terms of buffering --
---------------------------------------
elemDiscBuffering = FV1Buffer("cyt")	-- where buffering occurs
elemDiscBuffering:add_reaction(
	"clb",						    -- the buffering substance
	"ca_cyt",						-- the buffered substance
	totalClb,						-- total amount of buffer
	k_bind_clb,					    -- binding rate constant
	k_unbind_clb)				    -- unbinding rate constant

--[[ Calmodulin
elemDiscBuffering_clm = FV1Buffer("cyt")
elemDiscBuffering_clm:add_reaction(
	"clm_c",
	"ca_cyt",
	totalClmd,
	k_bind_clmd_c,
	k_unbind_clmd_c)
elemDiscBuffering_clm:add_reaction(
	"clm_n",
	"ca_cyt",				
	totalClmd,
	k_bind_clmd_n,
	k_unbind_clmd_n)
--]]

----------------------------------------------------
-- setup inner boundary (channels on ER membrane) --
----------------------------------------------------

-- We pass the function needed to evaluate the flux function here.
-- The order, in which the discrete fcts are passed, is crucial!
innerDisc = FV1InnerBoundaryCalciumER("ca_cyt, ca_er, ip3", "mem_er")

------------------------------
-- setup Neumann boundaries --
------------------------------
neumannDiscCA = NeumannBoundary("cyt")
neumannDiscCA:add(neumannCA, "ca_cyt", "mem_cyt" .. synapses)
neumannDiscIP3 = NeumannBoundary("cyt")
neumannDiscIP3:add(neumannIP3, "ip3", "mem_cyt" .. synapses)

--[[
-----------------------------------------------------------------
--  Setup Dirichlet Boundary
-----------------------------------------------------------------

dirichletBND = DirichletBoundary()
dirichletBND:add(dirichlet, "c", "Boundary, MembraneBnd")

membraneDirichletBND = DirichletBoundary()
membraneDirichletBND:add(membraneDirichlet, "c_membrane", "MembraneBnd")

--]]

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

-- ER flux
domainDisc:add(innerDisc)

--domainDisc:add(dirichletBND)
--domainDisc:add(membraneDirichletBND)

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
baseConvCheck:set_maximum_steps(10000)
baseConvCheck:set_minimum_defect(1e-28)
baseConvCheck:set_reduction(1e-2)
baseConvCheck:set_verbose(false)
base = LinearSolver()
base:set_convergence_check(baseConvCheck)
base:set_preconditioner(gs)

gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(timeDisc)
gmg:set_base_level(0)
gmg:set_base_solver(base)
gmg:set_smoother(gs)
gmg:set_cycle_type(1)
gmg:set_num_presmooth(3)
gmg:set_num_postsmooth(3)

-- biCGstab --
convCheck = ConvCheck()
convCheck:set_maximum_steps(50)
convCheck:set_minimum_defect(1e-24)
convCheck:set_reduction(1e-04)
convCheck:set_verbose(true)
bicgstabSolver = BiCGStab()
bicgstabSolver:set_preconditioner(gmg)
bicgstabSolver:set_convergence_check(convCheck)

-----------------------
-- non linear solver --
-----------------------
-- convergence check
newtonConvCheck = CompositeConvCheck3dCPU1(approxSpace)
newtonConvCheck:set_functions("")
newtonConvCheck:set_maximum_steps(20)
newtonConvCheck:set_minimum_defect({}, 1e-18)
newtonConvCheck:set_reduction({}, 1e-08)
newtonConvCheck:set_verbose(true)
newtonConvCheck:timeMeasurement(true)
--[[
newtonConvCheck = ConvCheck()
newtonConvCheck:set_maximum_steps(20)
newtonConvCheck:set_minimum_defect(1e-21)
newtonConvCheck:set_reduction(1e-08)
newtonConvCheck:set_verbose(true)
--]]

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

-- set initial values
Interpolate(CaCytStartValue, u, "ca_cyt", 0.0)
Interpolate(CaERStartValue, u, "ca_er", 0.0)
Interpolate(IP3StartValue, u, "ip3", 0.0)
Interpolate(ClbStartValue, u, "clb", 0.0)
--Interpolate(ClmCStartValue, u, "clm_c", 0.0)
--Interpolate(ClmNStartValue, u, "clm_n", 0.0)

-- timestep in seconds
dt = timeStep
time = 0.0
step = 0

-- filename for output
filename = "retic/result"

-- write start solution
out = VTKOutput()
out:print(filename, u, step, time)
takeMeasurement(u, approxSpace, time, "nuc", "ca_cyt", "measurements")
--exportSolution(u, approxSpace, time, "mem_cyt", "ca_cyt", "solution/solution");

-- some info output
print( "   numPreRefs is   " .. numPreRefs ..     ",  numRefs is         " .. numRefs)
print( "   nTimeSteps is " .. nTimeSteps)

-- create new grid function for old value
uOld = u:clone()

-- store grid function in vector of old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)

for step = 1, nTimeSteps do
	print("++++++ TIMESTEP " .. step .. " BEGIN ++++++")
	
	-- choose time step (currently constant)
	do_dt = dt
	
	-- setup time disc for old solutions and timestep
	timeDisc:prepare_step(solTimeSeries, do_dt)
	
	-- prepare Newton solver
	if newtonSolver:prepare(u) == false
	then print ("Newton solver failed at step "..step.."."); exit();
	end 
	
	-- apply Newton solver
	if newtonSolver:apply(u) == false
	then print ("Newton solver failed at step "..step.."."); exit();
	end 
	
	-- update new time
	time = solTimeSeries:time(0) + do_dt
	
	-- plot solution every 0.005 seconds of simulated time
	if math.abs(time/0.005 - math.floor(time/0.005+0.5)) < 1e-9
	then out:print(filename, u, math.floor(time/0.005+0.5), time)
	end
	
	-- take measurement in nucleus
	takeMeasurement(u, approxSpace, time, "nuc", "ca_cyt", "measurements")
	
	-- export solution of ca on mem_er
	--exportSolution(u, approxSpace, time, "mem_cyt", "ca_cyt", "solution/solution");
	
	-- get oldest solution
	oldestSol = solTimeSeries:oldest()
	
	-- copy values into oldest solution (we reuse the memory here)
	VecScaleAssign(oldestSol, 1.0, u)
	
	-- push oldest solutions with new values to front, oldest solution pointer is popped from end
	solTimeSeries:push_discard_oldest(oldestSol, time)
	
	print("++++++ TIMESTEP " .. step .. "  END ++++++");
end

-- end timeseries, produce gathering file and save
out:write_time_pvd(filename, u)
