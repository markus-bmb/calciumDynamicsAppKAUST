--------------------------------------------------------------
-- Bouton simulation for Ranjita                            --
--                                                          --
-- Author: Markus Breit                                     --
--------------------------------------------------------------

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
discreteRyR = util.HasParamOption("-discreteRyR")
defaultGrid = "../apps/calciumDynamics_app/grids/bouton.ugx"
if discreteRyR then
	defaultGrid = "../apps/calciumDynamics_app/grids/bouton_discreteRyR.ugx"
end

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
reqSubsets = {"cyt", "er", "pm", "az", "erm", "diri"}
dom = util.CreateDomain(gridName, numRefs, reqSubsets)


-- in parallel environments: use a load balancer to distribute the grid
balancer.partitioner = "parmetis"

-- protect ER membrane from being cut by partitioning
ccw = SubsetCommunicationWeights(dom)
ccw:set_weight_on_subset(100000.0, 4) -- ER membrane
balancer.communicationWeights = ccw

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
erVol = "er"
erMem = "erm"
erMemVec = {"erm"}
ryrSubset = "erm"
if discreteRyR then
	erMem = "erm, ryr"
	erMemVec = {"erm", "ryr"}
	ryrSubset = "ryr"
end
plMem = "pm, az"
plMem_vec = {"pm", "az"}
bnd = "diri"

outerDomain = cytVol .. ", " .. plMem .. ", " .. erMem .. ", " .. bnd
innerDomain = erVol .. ", " .. erMem

approxSpace:add_fct("ca_er", "Lagrange", 1, innerDomain)
approxSpace:add_fct("ca_cyt", "Lagrange", 1, outerDomain)
--approxSpace:add_fct("ip3", "Lagrange", 1, outerDomain)
approxSpace:add_fct("clb", "Lagrange", 1, outerDomain)

approxSpace:init_top_surface()
approxSpace:init_levels()
approxSpace:print_layout_statistic()
approxSpace:print_statistic()

--OrderCuthillMcKee(approxSpace, true) -- (not applicable)


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
k_bind_clb = 27.0e06
k_unbind_clb = 19

-- initial concentrations
ca_cyt_init = 5.0e-08
ca_er_init = 2.5e-4
--ip3_init = 4.0e-8
clb_init = totalClb / (k_bind_clb/k_unbind_clb*ca_cyt_init + 1)

-- reaction reate IP3
--reactionRateIP3 = 0.11

-- equilibrium concentration IP3
--equilibriumIP3 = 4.0e-08

-- reation term IP3
--reactionTermIP3 = -reactionRateIP3 * equilibriumIP3

--IP3Rdensity = 17.3
RYRdensity = 1.0 / compute_volume_of_subset(approxSpace, 4) -- 4 is ER membrane
if discreteRyR then
	RYRdensity = 1.0 / compute_volume_of_subset(approxSpace, 5) -- 5 is ryr subset
end
leakERconstant = 1.0e-17

-- this is a little bit more complicated, since it must be ensured that
-- the net flux for equilibrium concentrations is zero
-- MUST be adapted whenever any parameterization of ER flux mechanisms is changed!
local v_s = 6.5e-27						-- V_S param of SERCA pump
local k_s = 1.8e-7						-- K_S param of SERCA pump
--local j_ip3r = 3.7606194166520605e-23   -- single channel IP3R flux (mol/s) - to be determined via gdb
local j_ryr = 1.1201015633466695e-21    -- single channel RyR flux (mol/s) - to be determined via gdb
				  							-- ryr1: 1.1204582669024472e-21	
local j_leak = ca_er_init-ca_cyt_init	-- leak proportionality factor
	
SERCAdensity = leakERconstant * j_leak
--	SERCAdensity = SERCAdensity + IP3Rdensity * j_ip3r
SERCAdensity = SERCAdensity + RYRdensity * j_ryr
if discreteRyR then
	-- calculate as if RyR were continuously distributed
	ryr_area = compute_volume_of_subset(approxSpace, 5)
	erm_area = compute_volume_of_subset(approxSpace, 4)
	pseudoRYRdensity = RYRdensity * ryr_area / (ryr_area + erm_area)
	SERCAdensity = SERCAdensity + pseudoRYRdensity * j_ryr
end
SERCAdensity = SERCAdensity / (v_s/(k_s/ca_cyt_init+1.0)/ca_er_init)
print("SERCA density: " .. SERCAdensity)

pmcaDensity = 500.0
ncxDensity  = 15.0
vdccDensity = 1128.0

function pmLeakage(x, y, z, t, si)
	factor = 1e12 / (1.0-1e3*ca_cyt_init)
	if si == 3 then -- active zone
		return factor * (pmcaDensity * 6.9672131147540994e-24	-- single pump PMCA flux (mol/s)
				         + ncxDensity *  6.7567567567567566e-23	-- single pump NCX flux (mol/s)
				         + vdccDensity * (-1.5752042094823713e-25))    -- single channel VDCC flux (mol/s)
	end
	
	return factor * (pmcaDensity * 6.9672131147540994e-24	-- single pump PMCA flux (mol/s)
			         + ncxDensity *  6.7567567567567566e-23)	-- single pump NCX flux (mol/s)
end

-- AP firing pattern (20Hz / 1Hz, 20APs)
stimFreq = 20 -- 1
onset = 0.001
numAP = 20
function stimAPs(x, y, z, t, si)	
	-- spike train
	if (t >= onset and t < onset + 1.0/stimFreq * numAP) then
        phase = (t-onset) % (1.0/stimFreq)
        
        if phase < 0.001 or phase > 0.008 then
        	return -0.075
		
    	elseif phase < 0.0025 then
			return -0.075 + 8.0 * (phase - 0.001)
		
		elseif phase < 0.003 then
        	return -0.053 + 160.0 * (phase - 0.0025)
		
		elseif phase < 0.0045 then
        	return 0.027 - 62.0 * (phase - 0.003)
        	
        elseif phase < 0.007 then
        	return -0.066 + 0.4 * (phase - 0.0045)
		end
	end
		
	return -0.075
end


--------------------------
-- setup discretization --
--------------------------
-- diffusion --
elemDiscCYT = ConvectionDiffusion("ca_cyt", cytVol, "fv1")
elemDiscCYT:set_diffusion(D_cac)

elemDiscER = ConvectionDiffusion("ca_er", erVol, "fv1") 
elemDiscER:set_diffusion(D_cae)

--[[
elemDiscIP3 = ConvectionDiffusion("ip3", cytVol, "fv1")
elemDiscIP3:set_diffusion(D_ip3)
elemDiscIP3:set_reaction_rate(reactionRateIP3)
elemDiscIP3:set_reaction(reactionTermIP3)
--]]

elemDiscClb = ConvectionDiffusion("clb", cytVol, "fv1")
elemDiscClb:set_diffusion(D_clb)


-- buffering --
elemDiscBuffering = BufferFV1(cytVol)	-- where buffering occurs
elemDiscBuffering:add_reaction(
	"clb",						    -- the buffering substance
	"ca_cyt",						-- the buffered substance
	totalClb,						-- total amount of buffer
	k_bind_clb,					    -- binding rate constant
	k_unbind_clb)				    -- unbinding rate constant


-- er membrane transport systems

-- We pass the function needed to evaluate the flux function here.
-- The order, in which the discrete fcts are passed, is crucial!
--[[
ip3r = IP3R({"ca_cyt", "ca_er", "ip3"})
ip3r:set_scale_inputs({1e3,1e3,1e3})
ip3r:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)
--]]

--ryr = RyR({"ca_cyt", "ca_er"})
ryr = RyRinstat({"ca_cyt", "ca_er"}, {ryrSubset}, approxSpace)
ryr:set_scale_inputs({1e3,1e3})
ryr:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

serca = SERCA({"ca_cyt", "ca_er"})
serca:set_scale_inputs({1e3,1e3})
serca:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

leakER = Leak({"ca_er", "ca_cyt"})
leakER:set_scale_inputs({1e3,1e3})
leakER:set_scale_fluxes({1e3}) -- from mol/(m^2 s) to (mol um)/(dm^3 s)


--discIP3R = MembraneTransportFV1(erMem, ip3r)
--discIP3R:set_density_function(IP3Rdensity)

discRyR = MembraneTransportFV1(ryrSubset, ryr)
discRyR:set_density_function(RYRdensity)

discSERCA = MembraneTransportFV1(erMem, serca)
discSERCA:set_density_function(SERCAdensity)

discERLeak = MembraneTransportFV1(erMem, leakER)
discERLeak:set_density_function(1e12*leakERconstant/(1e3)) -- from mol/(um^2 s M) to m/s


-- plasma membrane transport systems
pmca = PMCA({"ca_cyt", ""})
pmca:set_constant(1, 1.0)
pmca:set_scale_inputs({1e3,1.0})
pmca:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

ncx = NCX({"ca_cyt", ""})
ncx:set_constant(1, 1.0)
ncx:set_scale_inputs({1e3,1.0})
ncx:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

leakPM = Leak({"", "ca_cyt"})
leakPM:set_constant(0, 1.0)
leakPM:set_scale_inputs({1.0,1e3})
leakPM:set_scale_fluxes({1e3}) -- from mol/(m^2 s) to (mol um)/(dm^3 s)

vdcc = VDCC_BG_UserData({"ca_cyt", ""}, {"az"}, approxSpace)
vdcc:set_constant(1, 1.5)
vdcc:set_scale_inputs({1e3,1.0})
vdcc:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)
vdcc:set_channel_type_N()
vdcc:set_potential_function("stimAPs")
vdcc:init(0.0)


discPMCA = MembraneTransportFV1(plMem, pmca)
discPMCA:set_density_function(pmcaDensity)

discNCX = MembraneTransportFV1(plMem, ncx)
discNCX:set_density_function(ncxDensity)

discPMLeak = MembraneTransportFV1(plMem, leakPM)
discPMLeak:set_density_function("pmLeakage")

discVDCC = MembraneTransportFV1(plMem, vdcc)
discVDCC:set_density_function(vdccDensity)


-- dirichlet bnd
diri = DirichletBoundary()
diri:add(ca_cyt_init, "ca_cyt", "diri")
diri:add(clb_init, "clb", "diri")


------------------------------------------
-- setup complete domain discretization --
------------------------------------------
domainDisc = DomainDiscretization(approxSpace)

domainDisc:add(elemDiscER)
domainDisc:add(elemDiscCYT)
--domainDisc:add(elemDiscIP3)
domainDisc:add(elemDiscClb)

domainDisc:add(elemDiscBuffering)

--domainDisc:add(discIP3R)
domainDisc:add(discRyR)
domainDisc:add(discSERCA)
domainDisc:add(discERLeak)

domainDisc:add(discPMCA)
domainDisc:add(discNCX)
domainDisc:add(discPMLeak)
domainDisc:add(discVDCC)

domainDisc:add(diri)

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
	
	ilu_gmg = ILU()
    ilu_gmg:set_sort(true)
	gmg:set_smoother(ilu_gmg)
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
--newtonConvCheck:set_component_check("ca_cyt, ca_er, clb", 1e-18, 1e-08)
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
InterpolateInner(ca_er_init, u, "ca_er", 0.0)
--InterpolateInner(ip3_init, u, "ip3", 0.0)
InterpolateInner(clb_init, u, "clb", 0.0)

-- timestep in seconds
dt = timeStepStart
time = 0.0
step = 0

--take_measurement(u, time, measZones, "ca_cyt, ip3, clb", fileName .. "meas/data")


-- create new grid function for old value
uOld = u:clone()

-- store grid function in vector of  old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)



-- calculate start solution
opStat = AssembledOperator()
opStat:set_discretization(domainDisc)
newtonSolver:init(opStat)
newtonConvCheck:set_maximum_steps(1)
for i = 1,20 do
	timeDisc:prepare_step(solTimeSeries, 0.0)
	newtonSolver:apply(u)
	oldestSol = solTimeSeries:oldest()
	VecScaleAssign(oldestSol, 1.0, u)
	solTimeSeries:push_discard_oldest(oldestSol, time)
end
newtonConvCheck:set_maximum_steps(10)
newtonSolver:init(op)


if generateVTKoutput then
	out = VTKOutput()
	out:print(fileName .. "vtk/result", u, step, time)
end

compute_volume(approxSpace, "cyt, er, az, erm")

min_dt = timeStep / math.pow(2,15)
cb_interval = 10
lv = startLv
levelUpDelay = 0.0;
cb_counter = {}
for i = 0, startLv do cb_counter[i] = 0 end
while endTime-time > 0.001*dt do
	print("++++++ POINT IN TIME  " .. math.floor((time+dt)/dt+0.5)*dt .. "s  BEGIN ++++++")
	
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
		
		--take_measurement(u, time, measZones, "ca_cyt, ip3, clb", fileName .. "meas/data")
		
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
