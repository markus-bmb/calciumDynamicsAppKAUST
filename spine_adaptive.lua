----------------------------------------------------------------
--  Example script for simulation on 3d spine model	      --
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
gridName = "spines/normSpine.ugx"
--gridName = "spines/movedApp.ugx"
--gridName = "spines/bigSpine.ugx"
--gridName = "spines/bigSpineMovedApp.ugx"
--gridName = "spines/bigSpineBigApp.ugx"
--gridName = "spines/bigSpineBigMovedApp.ugx"
--gridName = "spines/bigSpine_smallER_longDend_scaled.ugx"

-- refinements before distributing grid
numPreRefs = util.GetParamNumber("-numPreRefs", 0)

-- total refinements
numRefs    = util.GetParamNumber("-numRefs",    0)

-- choose number of time steps
nTimeSteps =  util.GetParamNumber("-nTimeSteps", 5)

-- choose length of time step
timeStep =  util.GetParamNumber("-tstep", 0.01)

-- chose plotting interval
plotStep = util.GetParamNumber("-pstep", 0.01)

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
ca_cyt_init = 4.0e-8
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
for i=9,9 do
	syns["start"..i] = 0.005*(i-9)
	syns["end"..i] = 0.005*(i-9)+0.01
end


-- burst of calcium influx for active synapses (~1200 ions)
function ourNeumannBndCA(x, y, z, t, si)
	if 	(si==9 and syns["start"..si]<t and t<=syns["end"..si])
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
	if 	(si==9 and syns["start"..si]+ip3EntryDelay<t
	     and t<=syns["start"..si]+ip3EntryDelay+ip3EntryDuration)
	then efflux = - 2.1e-5/1.188 * (1.0 - (t-syns["start"..si])/ip3EntryDuration)
	else efflux = 0.0
	end
    return true, efflux
end


-- PMCA pumps (h2b variant)
PMCA =
{
	Kd		= 6.0e-08,	-- mol*dm^-3 (Elwess et al.)
--  Kd		= 3.4e-07,	-- mol*dm^-3 (Graupner) 
	maxFlux	= 1.7e-23,	-- mol*s^-1
	density	= 3.0e03		-- �m^-2
}

function PMCAfluxDensity(caCyt)
	local gatingFactor = 1.0 / (PMCA['Kd']*PMCA['Kd']/caCyt/caCyt + 1.0)
	local dimCorrFactor = 1.0e15	-- for dimensional correction: to mol*�m^2*dm^-3*s^-1
	return PMCA['density'] * PMCA['maxFlux'] * gatingFactor * dimCorrFactor;
end

function PMCAfluxDensityDeriv(caCyt)
	local d_gatingFactor = 2*PMCA['Kd']*PMCA['Kd']*caCyt / math.pow(PMCA['Kd']*PMCA['Kd'] + caCyt*caCyt, 2)
	local dimCorrFactor = 1.0e15	-- for dimensional correction: to mol*�m^2*dm^-3*s^-1
	return PMCA['density'] * PMCA['maxFlux'] * d_gatingFactor * dimCorrFactor;
end

-- NCX pumps
-- params according to Graupner et al.
NCX =
{
	Kd		= 1.8e-06,	-- mol*dm^-3
	maxFlux	= 2.5e-21,	-- mol*s^-1
	density	= 1.0e02		-- �m^-2
}
	
function NCXfluxDensity(caCyt)
	local gatingFactor = 1.0 / (NCX['Kd']/caCyt + 1.0)
	local dimCorrFactor = 1.0e15	-- for dimensional correction: to mol*�m^2*dm^-3*s^-1
	return NCX['density'] * NCX['maxFlux'] * gatingFactor * dimCorrFactor;
end

function NCXfluxDensityDeriv(caCyt)
	local d_gatingFactor = NCX['Kd'] / ((NCX['Kd'] + caCyt) * (NCX['Kd'] + caCyt))
	local dimCorrFactor = 1.0e15	-- for dimensional correction: to mol*�m^2*dm^-3*s^-1
	return NCX['density'] * NCX['maxFlux'] * d_gatingFactor * dimCorrFactor;
end

-- plasma membrane leakage
function leakPM(x, y, z, t, si)
	return true, 2.1e-20
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
innerDomain = "er, app, mem_er, mem_app"
outerDomain = "cyt, head, neck, dend, mem_cyt, mem_er, mem_app"
synapses = ", syn"
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
neumannPMCA = LuaUserFunctionNumber("PMCAfluxDensity", 1, false)
	neumannPMCA:set_deriv(0, "PMCAfluxDensityDeriv")
neumannNCX = LuaUserFunctionNumber("NCXfluxDensity", 1, false)
	neumannNCX:set_deriv(0, "NCXfluxDensityDeriv")
neumannLeak = LuaCondUserNumber3d("leakPM")

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

elemDiscER = ConvectionDiffusion("ca_er", "er, app") 
elemDiscER:set_disc_scheme("fv1")
elemDiscER:set_diffusion(diffusionMatrixCAer)
elemDiscER:set_source(rhs)
elemDiscER:set_upwind(upwind)

elemDiscCYT = ConvectionDiffusion("ca_cyt", "cyt, head, neck, dend")
elemDiscCYT:set_disc_scheme("fv1")
elemDiscCYT:set_diffusion(diffusionMatrixCAcyt)
elemDiscCYT:set_source(rhs)
elemDiscCYT:set_upwind(upwind)

elemDiscIP3 = ConvectionDiffusion("ip3", "cyt, head, neck, dend")
elemDiscIP3:set_disc_scheme("fv1")
elemDiscIP3:set_diffusion(diffusionMatrixIP3)
elemDiscIP3:set_reaction_rate(reactionRateIP3)
elemDiscIP3:set_reaction(reactionTermIP3)
elemDiscIP3:set_source(rhs)
elemDiscIP3:set_upwind(upwind)

elemDiscClb = ConvectionDiffusion("clb", "cyt, head, neck, dend")
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
elemDiscBuffering = FV1Buffer("cyt, head, neck, dend")	-- where buffering occurs
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
innerDiscIP3R = FV1InnerBoundaryIP3R("ca_cyt, ca_er, ip3", "mem_er, mem_app")
innerDiscRyR = FV1InnerBoundaryRyR("ca_cyt, ca_er", "mem_er")
innerDiscSERCA = FV1InnerBoundarySERCA("ca_cyt, ca_er", "mem_er, mem_app")
innerDiscLeak = FV1InnerBoundaryERLeak("ca_cyt, ca_er", "mem_er, mem_app")

------------------------------
-- setup Neumann boundaries --
------------------------------
neumannDiscCA = NeumannBoundary("cyt, head, neck, dend")
neumannDiscCA:add(neumannCA, "ca_cyt", "mem_cyt" .. synapses)

--[[
neumannPMCA:set_input(0, elemDiscCYT:value())
--neumannDiscCA:add(neumannPMCA, "ca_cyt", "mem_cyt" .. synapses)
neumannNCX:set_input(0, elemDiscCYT:value())
--neumannDiscCA:add(neumannNCX, "ca_cyt", "mem_cyt" .. synapses)
neumannDiscCA:add(neumannLeak, "ca_cyt", "mem_cyt" .. synapses)
--]]
neumannDiscPMCA = FV1BoundaryPMCA("ca_cyt", "mem_cyt" .. synapses)
neumannDiscNCX = FV1BoundaryNCX("ca_cyt", "mem_cyt" .. synapses)
neumannDiscLeak = FV1BoundaryPMLeak("", "mem_cyt" .. synapses)

neumannDiscIP3 = NeumannBoundary("cyt, head, neck, dend")
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
domainDisc:add(neumannDiscPMCA)
domainDisc:add(neumannDiscNCX)
domainDisc:add(neumannDiscLeak)
domainDisc:add(neumannDiscIP3)

-- ER flux
domainDisc:add(innerDiscIP3R)
domainDisc:add(innerDiscRyR)
domainDisc:add(innerDiscSERCA)
domainDisc:add(innerDiscLeak)

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
baseConvCheck:set_maximum_steps(1000)
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
convCheck:set_reduction(1e-06)
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
newtonConvCheck:set_maximum_steps(10)
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

-- set initial value
Interpolate(CaCytStartValue, u, "ca_cyt", 0.0)
Interpolate(CaERStartValue, u, "ca_er", 0.0)
Interpolate(IP3StartValue, u, "ip3", 0.0)
Interpolate(ClbStartValue, u, "clb", 0.0)


-- timestep in seconds
dt = timeStep
time = 0.0
step = 0

-- filename
fileName = "normSpine/"
--fileName = "movedApp/"
--fileName = "bigSpine/"
--fileName = "bigSpineMovedApp/"
--fileName = "bigSpineBigApp/"
--fileName = "bigSpineBigMovedApp/"
--fileName = "bigSpine2/"

-- write start solution
print("Writing start values")
out = VTKOutput()
out:print(fileName .. "result", u, step, time)
takeMeasurement(u, approxSpace, time, "head", "ca_cyt", fileName .. "measurements_head")
takeMeasurement(u, approxSpace, time, "neck", "ca_cyt", fileName .. "measurements_neck")
takeMeasurement(u, approxSpace, time, "dend", "ca_cyt", fileName .. "measurements_dend")
--exportSolution(u, approxSpace, time, "mem_cyt", "ca_cyt", "solution/solution");

-- some info output
print( "   numPreRefs is   " .. numPreRefs ..     ",  numRefs is         " .. numRefs)
print( "   nTimeSteps is " .. nTimeSteps)

-- create new grid function for old value
uOld = u:clone()

-- store grid function in vector of  old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)


min_dt = timeStep / 2048.0
cb_interval = 10
lv = 0
cb_counter = {}
cb_counter[0] = 0
while time < timeStep*nTimeSteps do
	print("++++++ POINT IN TIME  " .. time+dt .. "s  BEGIN ++++++")
	
	-- setup time Disc for old solutions and timestep
	timeDisc:prepare_step(solTimeSeries, dt)
	
	-- prepare newton solver
	if newtonSolver:prepare(u) == false then print ("Newton solver failed at step "..step.."."); exit(); end 
	
	-- apply newton solver
	if newtonSolver:apply(u) == false
	then
		-- in case of failure:
		print ("Newton solver failed at point in time " .. time .. " with time step " .. dt)
		dt = dt/2
		lv = lv + 1
		VecScaleAssign(u, 1.0, solTimeSeries:latest())
		-- halve time step and try again unless time step below minimum
		if dt < min_dt
		then 
			print ("Time step below minimum. Aborting. Failed at point in time " .. time .. ".")
			time = timeStep*nTimeSteps
		else
			print ("Trying with half the time step...")
			cb_counter[lv] = 0
		end
	else
		-- update new time
		time = solTimeSeries:time(0) + dt
		
		-- update check-back counter and if applicable, reset dt
		cb_counter[lv] = cb_counter[lv] + 1
		while cb_counter[lv] % (2*cb_interval) == 0 and lv > 0 do
			dt = 2*dt;
			lv = lv - 1
			cb_counter[lv] = cb_counter[lv] + cb_counter[lv+1] / 2
			cb_counter[lv+1] = 0
		end
		
		-- plot solution every plotStep seconds
		if math.abs(time/plotStep - math.floor(time/plotStep+0.5)) < 1e-5
		then out:print(fileName .. "result", u, math.floor(time/plotStep+0.5), time)
		end
		
		-- take measurement in nucleus every timeStep seconds 
		--if math.abs(time/timeStep - math.floor(time/timeStep+0.5)) < 1e-5
		--then
		takeMeasurement(u, approxSpace, time, "head", "ca_cyt", fileName.."measurements_head")
		takeMeasurement(u, approxSpace, time, "neck", "ca_cyt", fileName.."measurements_neck")
		takeMeasurement(u, approxSpace, time, "dend", "ca_cyt", fileName.."measurements_dend")
		--end
				
		-- export solution of ca on mem_er
		--exportSolution(u, approxSpace, time, "mem_cyt", "ca_cyt", "solution/solution");
		
		-- get oldest solution
		oldestSol = solTimeSeries:oldest()
		
		-- copy values into oldest solution (we reuse the memory here)
		VecScaleAssign(oldestSol, 1.0, u)
		
		-- push oldest solutions with new values to front, oldest sol pointer is popped from end
		solTimeSeries:push_discard_oldest(oldestSol, time)
		
		print("++++++ POINT IN TIME  " .. time .. "s  END ++++++++");
	end

end

-- end timeseries, produce gathering file
out:write_time_pvd(fileName .. "result", u)