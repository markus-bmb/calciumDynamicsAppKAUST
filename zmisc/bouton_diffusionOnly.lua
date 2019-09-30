--------------------------------------------------------------------------------
-- Bouton example simulation for Ranjita                                      --
-- This script simulates calcium dynamics within a half-sperical bouton that  --
-- is activated by a VDCC influx modeled according to Ranjita Dutta-Roy.      --
-- Apart from the influx, only PMCA efflux and diffusion is considered and    --
-- there is a Dirichlet constraint at the boundary, which is probably a bad   --
-- idea.                                                                      --
--                                                                            --
-- Author: Markus Breit                                                       --
-- Date:   2018-02-09                                                         --
--------------------------------------------------------------------------------

-- for profiler output
--SetOutputProfileStats(true)

-- load pre-implemented lua functions
ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

-- choose dimension and algebra
InitUG(3, AlgebraType("CPU", 1));
 
-- choose outfile directory
fileName = util.GetParam("-outName", "boutonTest")
fileName = fileName.."/"

-- whether to use discrete RyR channel
defaultGrid = "../apps/calciumDynamics_app/grids/bouton_noER.ugx"

-- choice of grid
gridName = util.GetParam("-grid", defaultGrid)

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

-- specify -vtk to generate vtk output (and plotting interval)
generateVTKoutput = util.HasParamOption("-vtk")
plotStep = util.GetParamNumber("-pstep", 0.01)

-- choose solver setup
solverID = util.GetParam("-solver", "GS")
solverID = string.upper(solverID)
validSolverIDs = {}
validSolverIDs["GMG"] = 0;
validSolverIDs["GS"] = 0;
validSolverIDs["ILU"] = 0;
if (validSolverIDs[solverID] == nil) then
    error("Unknown solver identifier " .. solverID)
end

-- solver output verbosity
verbose = util.HasParamOption("-verbose")

-- stimulation frequency
stimFreq = util.GetParamNumber("-endTime", 20.0, "stimulation frequency (Hz)")



-------------------------------
-- setup approximation space --
-------------------------------
-- create, load, refine and distribute domain
print("create, refine and distribute domain")
reqSubsets = {"cyt", "pm", "az", "diri"}
dom = util.CreateDomain(gridName, numRefs, reqSubsets)

-- in parallel environments: use a load balancer to distribute the grid
balancer.partitioner = "parmetis"
balancer.staticProcHierarchy = true
balancer.firstDistLvl = -1
balancer.redistSteps = 0

balancer.ParseParameters()
balancer.PrintParameters()
loadBalancer = balancer.CreateLoadBalancer(dom)
if loadBalancer ~= nil then
	balancer.Rebalance(dom, loadBalancer)
	print("Edge cut on base level: "..balancer.defaultPartitioner:edge_cut_on_lvl(0))
	loadBalancer:estimate_distribution_quality()
	loadBalancer:print_quality_records()
end

print(dom:domain_info():to_string())
--[[
--print("Saving domain grid and hierarchy.")
--SaveDomain(dom, "refined_grid_p" .. ProcRank() .. ".ugx")
--SaveGridHierarchyTransformed(dom:grid(), "refined_grid_hierarchy_p" .. ProcRank() .. ".ugx", 2.0)
print("Saving parallel grid layout")
SaveParallelGridLayout(dom:grid(), "parallel_grid_layout_p"..ProcRank()..".ugx", 2.0)
--]]


-- create approximation space
approxSpace = ApproximationSpace(dom)

cytVol = "cyt"
plMem = "pm, az"
plMem_vec = {"pm", "az"}
bnd = "diri"

outerDomain = cytVol .. ", " .. plMem .. ", " .. bnd

approxSpace:add_fct("ca_cyt", "Lagrange", 1)

approxSpace:init_top_surface()
approxSpace:init_levels()
approxSpace:print_layout_statistic()
approxSpace:print_statistic()

OrderCuthillMcKee(approxSpace, true)

---------------
-- constants --
---------------
-- diffusion coefficients
D_cac = 220.0

-- initial concentrations
ca_cyt_init = 5.0e-08

vdccDensity = 1128.0
pmcaDensity = 500.0

-- iso potential is assumed; calculated from HH point neuron model
gk = 360.0  -- C/(V s m^2)
gna = 1200.0  -- C/(V s m^2)
gl = 3.0  -- C/(V s m^2)

ek = -0.077  -- V
ena = 0.05  -- V

cm = 0.01  -- F/m^2

ie = 0.7  -- C / (m^2 s)

function alpha_n(v)
	local x = -(v + 0.055)
	if (math.abs(x) > 1e-10) then 
		return 1e4*x / (math.exp(100.0*x) - 1.0)
	end
	return 1e4 * (0.01 - 0.5*x)
end
function beta_n(v)
	return 125.0 * math.exp(-(v + 0.065) / 0.08)
end

function alpha_m(v)
	local x = -(v + 0.04)
	if (math.abs(x) > 1e-10) then
  		return 1e5*x / (math.exp(100.0*x) - 1.0)
	end
	return 1e5 * (0.01 - 0.5*x);
end
function beta_m(v)
	return 4e3 * math.exp(-(v + 0.065) / 0.018)
end
function alpha_h(v)
	return 70.0 * math.exp(-(v + 0.065) / 0.02)
end
function beta_h(v)
	return 1e3 / (math.exp(-(v + 0.035) / 0.01) + 1.0)
end

alpha10 = 4.04e3  -- 1/s
beta10 = 2.88e3
alpha20 = 6.70e3
beta20 = 6.30e3
alpha30 = 4.39e3
beta30 = 8.16e3
alpha40 = 17.33e3
beta40 = 1.84e3

V1 = 0.04914  -- V
V2 = 0.04208
V3 = 0.05531
V4 = 0.02655

g_O = 3e-13  -- S

-- initial (equilibrium) conditions
vm = -0.07  -- V  
n = alpha_n(vm) / (alpha_n(vm) + beta_n(vm))
m = alpha_m(vm) / (alpha_m(vm) + beta_m(vm))
h = alpha_h(vm) / (alpha_h(vm) + beta_h(vm))

local alpha1 = alpha10*math.exp(vm/V1)
local beta1 = beta10*math.exp(-vm/V1)
local alpha2 = alpha20*math.exp(vm/V2)
local beta2 = beta20*math.exp(-vm/V2)
local alpha3 = alpha30*math.exp(vm/V3)
local beta3 = beta30*math.exp(-vm/V3)
local alpha4 = alpha40*math.exp(vm/V4)
local beta4 = beta40*math.exp(-vm/V4)
K1 = beta1/alpha1
K2 = beta2/alpha2
K3 = beta3/alpha3
K4 = beta4/alpha4

denom = 1.0 /(1.0 + K4 + K3*K4 + K2*K3*K4 + K1*K2*K3*K4)
C0 = K1*K2*K3*K4 * denom
C1 = K2*K3*K4 * denom
C2 = K3*K4 * denom
C3 = K4 * denom
O = denom

el = vm + gk/gl*n*n*n*n*(vm - ek) + gna/gl*m*m*m*h*(vm - ena)
print("el = " .. el)


nChannel = 40
azArea = compute_volume_of_subset(approxSpace, 2)
function vdccInflux(x, y, z, t, si)
	if t < 0.005 then -- only one AP
		return O * g_O * (0.138 - vm) / (2*96485.0) * nChannel / azArea * 1e15
	end
	return 0
end


--------------------------
-- setup discretization --
--------------------------
-- diffusion --
elemDiscCYT = ConvectionDiffusion("ca_cyt", cytVol, "fv1")
elemDiscCYT:set_diffusion(D_cac)


-- plasma membrane transport systems
--[[
vdcc = VDCC_BG_UserData({"ca_cyt", ""}, {"az"}, approxSpace)
vdcc:set_constant(1, 1.5)
vdcc:set_scale_inputs({1e3,1.0})
vdcc:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)
vdcc:set_channel_type_N()
vdcc:set_potential_function("stimAPs")
vdcc:init(0.0)

discVDCC = MembraneTransportFV1(plMem, vdcc)
discVDCC:set_density_function(vdccDensity)
--]]
caInflux = UserFluxBoundaryFV1("ca_cyt", "az")
caInflux:set_flux_function("vdccInflux")


pmca = PMCA({"ca_cyt", ""})
pmca:set_constant(1, 1.0)
pmca:set_scale_inputs({1e3,1.0})
pmca:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

discPMCA = MembraneTransportFV1(plMem, pmca)
discPMCA:set_density_function(pmcaDensity)

-- dirichlet bnd
diri = DirichletBoundary()
diri:add(ca_cyt_init, "ca_cyt", "diri")


------------------------------------------
-- setup complete domain discretization --
------------------------------------------
domainDisc = DomainDiscretization(approxSpace)

domainDisc:add(elemDiscCYT)
--domainDisc:add(discVDCC)
domainDisc:add(caInflux)
--domainDisc:add(diri)

domainDisc:add(discPMCA)

-- setup time discretization --
timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(1.0) -- 1.0 is implicit Euler

-- create operator from discretization
op = AssembledOperator()
op:set_discretization(timeDisc)
op:init()

------------------
-- solver setup --
------------------
-- debug writer
dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_base_dir(fileName)
dbgWriter:set_vtk_output(false)

-- biCGstab --
convCheck = ConvCheck()
convCheck:set_minimum_defect(1e-50)
convCheck:set_reduction(1e-6)
convCheck:set_verbose(verbose)

if (solverID == "ILU") then
    bcgs_steps = 10000
    ilu = ILU()
    ilu:set_sort(true)
    bcgs_precond = ilu
elseif (solverID == "GS") then
    bcgs_steps = 10000
    bcgs_precond = GaussSeidel()
else -- (solverID == "GMG")
	gmg = GeometricMultiGrid(approxSpace)
	gmg:set_discretization(domainDisc)
	gmg:set_base_level(0)
	gmg:set_gathered_base_solver_if_ambiguous(true)
	
	-- treat SuperLU problems with Dirichlet constraints by using constrained version
	gmg:set_base_solver(SuperLU())
	
	--ilu_gmg = ILU()
    --ilu_gmg:set_sort(true)
	gmg:set_smoother(GaussSeidel())
	gmg:set_smooth_on_surface_rim(true)
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(3)
	gmg:set_num_postsmooth(3)
	--gmg:set_rap(true)
	--gmg:set_debug(GridFunctionDebugWriter(approxSpace))
	
    bcgs_steps = 100
	bcgs_precond = gmg
end

convCheck:set_maximum_steps(bcgs_steps)

bicgstabSolver = BiCGStab()
bicgstabSolver:set_preconditioner(bcgs_precond)
bicgstabSolver:set_convergence_check(convCheck)
--bicgstabSolver:set_debug(dbgWriter)

--- non-linear solver ---
-- convergence check
newtonConvCheck = CompositeConvCheck(approxSpace, 10, 1e-16, 1e-08)
newtonConvCheck:set_verbose(true)
newtonConvCheck:set_time_measurement(true)
--newtonConvCheck:set_adaptive(true)

-- Newton solver
newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(bicgstabSolver)
newtonSolver:set_convergence_check(newtonConvCheck)
--newtonSolver:set_debug(dbgWriter)

newtonSolver:init(op)


-------------
-- solving --
-------------

-- get grid function
u = GridFunction(approxSpace)

-- set initial value
InterpolateInner(ca_cyt_init, u, "ca_cyt", 0.0)

-- timestep in seconds
dt = timeStepStart
time = 0.0
step = 0


-- create new grid function for old value
uOld = u:clone()

-- store grid function in vector of  old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)

if generateVTKoutput then
	out = VTKOutput()
	out:print(fileName .. "vtk/result", u, step, time)
end

compute_volume(approxSpace, "cyt, az")

min_dt = timeStep / math.pow(2,15)
cb_interval = 10
lv = startLv
levelUpDelay = 0.0;
cb_counter = {}
for i = 0, startLv do cb_counter[i] = 0 end
while endTime-time > 0.001*dt do
	print("++++++ POINT IN TIME  " .. math.floor((time+dt)/dt+0.5)*dt .. "s  BEGIN ++++++")
	
	-- update potential
	smallDT = dt
	if smallDT > 1e-5 then
		smallDT = 1e-5
	end
	
	smallT = 0.0
	while smallT < dt do
		if smallT + smallDT > dt then
			smallDT = dt - smallT
			smallT = dt
		else
			smallT = smallT + smallDT
		end
		
		n = n + smallDT*(alpha_n(vm) * (1.0 - n) - beta_n(vm) * n)
		m = m + smallDT*(alpha_m(vm) * (1.0 - m) - beta_m(vm) * m)
		h = h + smallDT*(alpha_h(vm) * (1.0 - h) - beta_h(vm) * h)
		vm = vm + smallDT/cm * (ie - (gk * n*n*n*n * (vm - ek) + gna * m*m*m * h * (vm - ena) + gl * (vm - el)))
		
		-- update Ca entry
		local alpha1 = alpha10*math.exp(vm/V1)
		local beta1 = beta10*math.exp(-vm/V1)
		local alpha2 = alpha20*math.exp(vm/V2)
		local beta2 = beta20*math.exp(-vm/V2)
		local alpha3 = alpha30*math.exp(vm/V3)
		local beta3 = beta30*math.exp(-vm/V3)
		local alpha4 = alpha40*math.exp(vm/V4)
		local beta4 = beta40*math.exp(-vm/V4)
		
		local dC0dt = -alpha1*C0 + beta1*C1
		local dC1dt = alpha1*C0 - beta1*C1 - alpha2*C1 + beta2*C2
		local dC2dt = alpha2*C1 - beta2*C2 - alpha3*C2 + beta3*C3
		local dC3dt = alpha3*C2 - beta3*C3 - alpha4*C3 + beta4*O
		local dOdt = alpha4*C3 - beta4*O
		
		C0 = C0 + smallDT*dC0dt
		C1 = C1 + smallDT*dC1dt
		C2 = C2 + smallDT*dC2dt
		C3 = C3 + smallDT*dC3dt
		O = O + smallDT*dOdt
	end
	
	--[[
	print("vm = " .. vm)
	print("n = " .. n)
	print("m = " .. m)
	print("h = " .. h)
	print("C0 = " .. C0)
	print("C1 = " .. C1)
	print("C2 = " .. C2)
	print("C3 = " .. C3)
	print("O = " .. O)
	--]]
	
	-- setup time Disc for old solutions and timestep
	timeDisc:prepare_step(solTimeSeries, dt)
	
	-- apply newton solver
	if newtonSolver:apply(u) == false
	then
		-- in case of failure:
		print ("Newton solver failed at point in time " .. time .. " with time step " .. dt)
		
		-- end timeseries, produce gathering file
		if (generateVTKoutput) then out:write_time_pvd(fileName .. "vtk/result", u) end
		
		exit()
		
		dt = dt/2
		lv = lv + 1
		VecScaleAssign(u, 1.0, solTimeSeries:latest())
		
		-- halve time step and try again unless time step below minimum
		if dt < min_dt
		then 
			print ("Time step below minimum. Aborting. Failed at point in time " .. time .. ".")
			time = endTime
		else
			print ("Trying with half the time step...")
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
		
		-- plot solution every plotStep seconds
		if (generateVTKoutput) then
			if math.abs(time/plotStep - math.floor(time/plotStep+0.5)) < 1e-5 then
				out:print(fileName .. "vtk/result", u, math.floor(time/plotStep+0.5), time)
			end
		end
		
		-- get oldest solution
		oldestSol = solTimeSeries:oldest()
		
		-- copy values into oldest solution (we reuse the memory here)
		VecScaleAssign(oldestSol, 1.0, u)
		
		-- push oldest solutions with new values to front, oldest sol pointer is popped from end
		solTimeSeries:push_discard_oldest(oldestSol, time)
		
		print("++++++ POINT IN TIME  " .. math.floor(time/dt+0.5)*dt .. "s  END ++++++++");
	end

end

-- end timeseries, produce gathering file
if (generateVTKoutput) then out:write_time_pvd(fileName .. "vtk/result", u) end
