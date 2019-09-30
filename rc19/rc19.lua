--------------------------------------------------------------------------------
-- This script performs a 3d calcium simulation on a reconstructed plasma     --
-- membrane of a dendrite and soma morphology that is augmented by an ER.     --
-- It is supposed to recreate experiments published in:                       --
-- Branco et al.: "Dendritic discrimination of temporal input sequences in    --
-- cortical neurons", Science (2010).                                         --
--                                                                            --
-- Author: Markus Breit                                                       --
-- Date:   2012-12-31                                                         --
--------------------------------------------------------------------------------

-- load pre-implemented lua functions
ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")
AssertPluginsLoaded({"neuro_collection", "ConvectionDiffusion", "Parmetis"})

-- choose dimension and algebra
InitUG(3, AlgebraType("CPU", 1))

-- choice of grid
gridName = util.GetParam("-grid", "calciumDynamics_app/grids/rc19amp.ugx")
--gridName = "simple_reticulum_3d.ugx"  -- for testing

-- refinements before distributing grid
numPreRefs = util.GetParamNumber("-numPreRefs", 0)

-- total refinements
numRefs = util.GetParamNumber("-numRefs", 0)

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

-- choose solver setup
solverID = util.GetParam("-solver", "GMG")
solverID = string.upper(solverID)
validSolverIDs = {}
validSolverIDs["GMG"] = 0
validSolverIDs["GS"] = 0
validSolverIDs["ILU"] = 0
if (validSolverIDs[solverID] == nil) then
    error("Unknown solver identifier " .. solverID)
end

-- choose order of synaptic activity
order = util.GetParamNumber("-synOrder", 0)

-- choose time between synaptic stimulations (in ms)
jump = util.GetParamNumber("-jumpTime", 5)

-- specify "-verbose" to output linear solver convergence
verbose	= util.HasParamOption("-verbose")

-- choose outfile directory
outPath = util.GetParam("-outName", "rc19test")
outPath = outPath.."/"

-- specify -vtk to generate vtk output
generateVTKoutput = util.HasParamOption("-vtk")

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


-- reaction reate IP3
reactionRateIP3 = 0.11

-- equilibrium concentration IP3
equilibriumIP3 = 4.0e-08

-- reation term IP3
reactionTermIP3 = -reactionRateIP3 * equilibriumIP3

---------------------------------------------------------------------
-- functions steering tempo-spatial parameterization of simulation --
---------------------------------------------------------------------
-- density correction factor (simulates larger surface area of ER caused by stacking etc.)
dcf = 1.0

-- project coordinates on dendritic length from soma (approx.)
function dendLengthPos(x,y,z)
	return (0.91*(x-1.35) +0.40*(y+6.4) -0.04*(z+2.9)) / 111.0
end

function IP3Rdensity(x,y,z,t,si)
	local dens = math.abs(dendLengthPos(x,y,z))
	
	if math.abs(dens) > 1e3 or dens~=dens then
		print("problem: d = " .. dens .. ". (x="..x..", y="..y..", z="..z..")")
	end
	
	-- fourth order polynomial, distance to soma
	dens = 1.4 -2.8*dens +6.6*math.pow(dens,2) -7.0*math.pow(dens,3) +2.8*math.pow(dens,4)
	dens = dens * 17.3
	-- cluster for branching points
	if (si==22 or si==23 or si==24) then dens = dens * 10 end 
	return dens
end

function RYRdensity(x,y,z,t,si)
	local dens = math.abs(dendLengthPos(x,y,z))
	-- fourth order polynomial, distance to soma
	dens = 1.5 -3.5*dens +9.1*math.pow(dens,2) -10.5*math.pow(dens,3) +4.3*math.pow(dens,4)
	dens = dens * 0.86; 
	return dens
end

leakERconstant = 3.8e-17

-- this is a little bit more complicated, since it must be ensured that
-- the net flux for equilibrium concentrations is zero
-- MUST be adapted whenever any parameterization of ER flux mechanisms is changed!
function SERCAdensity(x,y,z,t,si)
	local v_s = 6.5e-27						-- V_S param of SERCA pump
	local k_s = 1.8e-7						-- K_S param of SERCA pump
	local j_ip3r = 3.7606194166520605e-23	-- single channel IP3R flux (mol/s) - to be determined via gdb
	local j_ryr = 1.1204582669024472e-21	-- single channel RyR flux (mol/s) - to be determined via gdb
	local j_leak = ca_er_init-ca_cyt_init	-- leak proportionality factor
	
	local dens =  IP3Rdensity(x,y,z,t,si) * j_ip3r
				+ RYRdensity(x,y,z,t,si) * j_ryr
				+ leakERconstant * j_leak
	dens = dens / (v_s/(k_s/ca_cyt_init+1.0)/ca_er_init)
	return dens
end

function LEAKERconstant(x,y,z,t,si)
	return dcf*3.4e-08
end

pmcaDensity = 500.0
ncxDensity  = 15.0
vgccDensity = 1.0

leakPMconstant =  pmcaDensity * 6.9672131147540994e-24	-- single pump PMCA flux (mol/s)
				+ ncxDensity *  6.7567567567567566e-23	-- single pump NCX flux (mol/s)
				+ vgccDensity * (-1.5752042094823713e-25)    -- single channel VGCC flux (mol/s)
				-- *1.5 // * 0.5 for L-type // T-type
if (leakPMconstant < 0) then error("PM leak flux is outward for these density settings!") end


-- firing pattern of the synapses
syns = {}
synStart = 6
synStop = 13
caEntryDuration = 0.01
if (order==0) then
	for i=synStart,synStop do
		syns[i] = jump/1000.0*(i-synStart)
	end
else
	for i=synStart,synStop do
		syns[i] = jump/1000.0*(synStop-i)
	end
end

-- burst of calcium influx for active synapses (~1200 ions)
function ourNeumannBndCA(x, y, z, t, si)
	local efflux = 0.0
	if si >= synStart and si <= synStop and syns[si] < t and t <= syns[si] + caEntryDuration then
		-- efflux = -5e-6 * 11.0/16.0*(1.0+5.0/((10.0*(t-syns["start"..si])+1)*(10.0*(t-syns["start"..si])+1)))
		efflux = -2e-4
	end
    return true, efflux
end


-- burst of ip3 at active synapses (triangular, immediate)
ip3EntryDelay = 0.000
ip3EntryDuration = 2.0
function ourNeumannBndIP3(x, y, z, t, si)
	local efflux = 0.0
	if si >= synStart and si <= synStop and syns[si] + ip3EntryDelay < t and t <= syns[si] + ip3EntryDelay + ip3EntryDuration then
		efflux = - 2.1e-5/1.188 * (1.0 - (t-syns[si])/ip3EntryDuration)
	end
    return true, efflux
end


-------------------------------
-- setup approximation space --
-------------------------------
-- create, load, refine and distribute domain
dom = util.CreateDomain(gridName, 0)

-- in parallel environments: use a load balancer to distribute the grid
balancer.partitioner = "parmetis"
balancer.staticProcHierarchy = true
balancer.firstDistLvl = 0
balancer.firstDistProcs = 64
balancer.redistSteps = 2
balancer.redistProcs = 32

loadBalancer = balancer.CreateLoadBalancer(dom)
if loadBalancer ~= nil then
	mu = ManifoldUnificator(dom)
	mu:add_protectable_subsets("mem_er")
	cdgm = ClusteredDualGraphManager()
	cdgm:add_unificator(SiblingUnificator())
	cdgm:add_unificator(mu)
	balancer.defaultPartitioner:set_dual_graph_manager(cdgm)
	balancer.qualityRecordName = "coarse"
	balancer.Rebalance(dom, loadBalancer)
end

print(dom:domain_info():to_string())

if load_balancer ~= nil then
	loadBalancer:print_quality_records()
end

--[[
SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), "refined_grid_hierarchy_p" .. ProcRank() .. ".ugx", 20.0)
SaveParallelGridLayout(dom:grid(), "parallel_grid_layout_p"..ProcRank()..".ugx", 20.0)
--]]

-- create approximation space
cytVol = "cyt"
measZones = ""
--for i=1,15 do
--	measZones = measZones .. ", measZone" .. i
--end
cytVol = cytVol .. measZones

nucVol = "nuc"
nucMem = "mem_nuc"

erVol = "er"

plMem = "mem_cyt"
plMem_vec = {"mem_cyt"}
synapses = ""
for i=1,8 do
	synapses = synapses .. ", syn" .. i
	plMem_vec[#plMem_vec+1] = "syn"..i
end
plMem = plMem .. synapses

erMem = "mem_er"
measZonesERM = "measZoneERM"..1
for i=2,15 do
	measZonesERM = measZonesERM .. ", measZoneERM" .. i
end
erMem = erMem .. ", " .. measZonesERM

outerDomain = cytVol .. ", " .. nucVol .. ", " .. nucMem .. ", " .. plMem .. ", " .. erMem
innerDomain = erVol .. ", " .. erMem

approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("ca_er", "Lagrange", 1, innerDomain)
approxSpace:add_fct("ca_cyt", "Lagrange", 1, outerDomain)
approxSpace:add_fct("ip3", "Lagrange", 1, outerDomain)
approxSpace:add_fct("clb", "Lagrange", 1, outerDomain)

approxSpace:init_levels()
approxSpace:print_statistic()


----------------------------------------------------------
-- setup FV convection-diffusion element discretization --
----------------------------------------------------------
elemDiscER = ConvectionDiffusionFV1("ca_er", erVol) 
elemDiscER:set_diffusion(D_cae)

elemDiscCYT = ConvectionDiffusionFV1("ca_cyt", cytVol..", "..nucVol)
elemDiscCYT:set_diffusion(D_cac)

elemDiscIP3 = ConvectionDiffusionFV1("ip3", cytVol..", "..nucVol)
elemDiscIP3:set_diffusion(D_ip3)
elemDiscIP3:set_reaction_rate(reactionRateIP3)
elemDiscIP3:set_reaction(reactionTermIP3)

elemDiscClb = ConvectionDiffusionFV1("clb", cytVol..", "..nucVol)
elemDiscClb:set_diffusion(D_clb)


---------------------------------------
-- setup reaction terms of buffering --
---------------------------------------
elemDiscBuffering = BufferFV1(cytVol)	-- where buffering occurs
elemDiscBuffering:add_reaction(
	"clb",						    -- the buffering substance
	"ca_cyt",						-- the buffered substance
	totalClb,						-- total amount of buffer
	k_bind_clb,					    -- binding rate constant
	k_unbind_clb)				    -- unbinding rate constant


----------------------------------------------------
-- setup inner boundary (channels on ER membrane) --
----------------------------------------------------
-- We pass the function needed to evaluate the flux function here.
-- The order, in which the discrete fcts are passed, is crucial!
ip3r = IP3R({"ca_cyt", "ca_er", "ip3"})
ip3r:set_scale_inputs({1e3,1e3,1e3})
ip3r:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

ryr = RyR({"ca_cyt", "ca_er"})
ryr:set_scale_inputs({1e3,1e3})
ryr:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

serca = SERCA({"ca_cyt", "ca_er"})
serca:set_scale_inputs({1e3,1e3})
serca:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

leakER = Leak({"ca_er", "ca_cyt"})
leakER:set_scale_inputs({1e3,1e3})
leakER:set_scale_fluxes({1e3}) -- from mol/(m^2 s) to (mol um)/(dm^3 s)


innerDiscIP3R = MembraneTransportFV1(erMem, ip3r)
innerDiscIP3R:set_density_function("IP3Rdensity")

innerDiscRyR = MembraneTransportFV1(erMem, ryr)
innerDiscRyR:set_density_function("RYRdensity")

innerDiscSERCA = MembraneTransportFV1(erMem, serca)
innerDiscSERCA:set_density_function("SERCAdensity")

innerDiscLeak = MembraneTransportFV1(erMem, leakER)
innerDiscLeak:set_density_function(1e12*leakERconstant/(1e3)) -- from mol/(um^2 s M) to m/s

----------------------------
-- setup outer boundaries --
----------------------------
-- synaptic activity
neumannDiscCA = NeumannBoundaryFV1("ca_cyt")
neumannDiscCA:add("ourNeumannBndCA", plMem,"cyt")
neumannDiscIP3 = NeumannBoundaryFV1("ip3")
neumannDiscIP3:add("ourNeumannBndIP3", plMem, "cyt")

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

vdcc = VDCC_BG_VM2UG({"ca_cyt", ""}, plMem_vec, approxSpace,
					 "neuronRes/timestep".."_order".. 0 .."_jump"..string.format("%1.1f", 5.0).."_",
					 "%.3f", ".dat", false)
vdcc:set_constant(1, 1.5)
vdcc:set_scale_inputs({1e3,1.0})
vdcc:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)
vdcc:set_channel_type_L() --default, but to be sure
vdcc:set_file_times(0.001, 0.0)
vdcc:init(0.0)

neumannDiscPMCA = MembraneTransportFV1(plMem, pmca)
neumannDiscPMCA:set_density_function(pmcaDensity)

neumannDiscNCX = MembraneTransportFV1(plMem, ncx)
neumannDiscNCX:set_density_function(ncxDensity)

neumannDiscLeak = MembraneTransportFV1(plMem, leakPM)
neumannDiscLeak:set_density_function(1e12*leakPMconstant / (1.0-1e3*ca_cyt_init))

neumannDiscVGCC = MembraneTransportFV1(plMem, vdcc)
neumannDiscVGCC:set_density_function(vgccDensity)


------------------------------------------
-- setup complete domain discretization --
------------------------------------------
domainDisc = DomainDiscretization(approxSpace)

-- diffusion discretizations
domainDisc:add(elemDiscER)
domainDisc:add(elemDiscCYT)
domainDisc:add(elemDiscIP3)
domainDisc:add(elemDiscClb)

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
-- debug writer
dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_vtk_output(false)

-- gmg
gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(timeDisc)
gmg:set_base_level(0)
gmg:set_base_solver(SuperLU())
gmg:set_gathered_base_solver_if_ambiguous(true)
gmg:set_smoother(GaussSeidel())
gmg:set_smooth_on_surface_rim(true)
gmg:set_cycle_type(1)
gmg:set_num_presmooth(3)
gmg:set_num_postsmooth(3)
gmg:set_rap(true)
--gmg:set_debug(dbgWriter)

-- biCGstab --
convCheck = ConvCheck()
convCheck:set_minimum_defect(1e-24)
convCheck:set_reduction(1e-06)
convCheck:set_verbose(verbose)
convCheck:set_maximum_steps(10000)

linearSolver = BiCGStab()
if solverID == "GS" then
	linearSolver:set_preconditioner(GaussSeidel())
elseif solverID == "ILU" then
	linearSolver:set_preconditioner(ILU())
else
	linearSolver:set_preconditioner(gmg)
	convCheck:set_maximum_steps(100)
end
linearSolver:set_convergence_check(convCheck)


-- nonlinear convergence check
newtonConvCheck = CompositeConvCheck(approxSpace, 10, 1e-21, 1e-08)
newtonConvCheck:set_verbose(true)
newtonConvCheck:set_time_measurement(true)
newtonConvCheck:set_adaptive(true)

-- Newton solver
newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(linearSolver)
newtonSolver:set_convergence_check(newtonConvCheck)

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
Interpolate(clb_init, u, "clb", 0.0)


-- timestep in seconds
dt = timeStepStart
time = 0.0
step = 0


if generateVTKoutput then
	out = VTKOutput()
	out:print(outPath .. "vtk/result", u, step, time)
end


take_measurement(u, time, measZonesERM, "ca_cyt, ca_er, ip3, clb", outPath .. "meas/data")


-- create new grid function for old value
uOld = u:clone()

-- store grid function in vector of  old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)


min_dt = timeStep / math.pow(2,15)
cb_interval = 10
lv = startLv
levelUpDelay = (synStop-synStart)*jump/1000.0 + caEntryDuration
cb_counter = {}
for i=0,startLv do cb_counter[i]=0 end
while endTime-time > 0.001*dt do
	print("++++++ POINT IN TIME  " .. math.floor((time+dt)/dt+0.5)*dt .. "s  BEGIN ++++++")
	
	-- setup time Disc for old solutions and timestep
	timeDisc:prepare_step(solTimeSeries, dt)
	
	-- apply newton solver
	if not newtonSolver:apply(u) then
		-- in case of failure:
		print ("Newton solver failed at point in time " .. time .. " with time step " .. dt)
		
		dt = dt/2
		lv = lv + 1
		VecScaleAssign(u, 1.0, solTimeSeries:latest())
		
		-- halve time step and try again unless time step below minimum
		if dt < min_dt
		then 
			print ("Time step below minimum. Aborting. Failed at point in time " .. time .. ".")
			time = endTime
		else
			print ("Trying with half the time step " .. dt)
			cb_counter[lv] = 0
		end
	else
		-- update new time
		time = solTimeSeries:time(0) + dt
		
		-- update check-back counter and if applicable, reset dt
		cb_counter[lv] = cb_counter[lv] + 1
		while cb_counter[lv] % (2*cb_interval) == 0 and lv > 0 and (time >= levelUpDelay or lv > startLv) do
			print ("Doubling time due to continuing convergence; now: " .. 2*dt)
			dt = 2*dt
			lv = lv - 1
			cb_counter[lv] = cb_counter[lv] + cb_counter[lv+1] / 2
			cb_counter[lv+1] = 0
		end
		
		-- plot solution every plotStep seconds
		if generateVTKoutput then
			if math.abs(time/plotStep - math.floor(time/plotStep+0.5)) < 1e-5 then
				out:print(outPath .. "vtk/result", u, math.floor(time/plotStep+0.5), time)
				out:write_time_pvd(outPath .. "vtk/result", u)
			end
		end

		-- take measurement		
		take_measurement(u, time, measZonesERM, "ca_cyt, ca_er, ip3, clb", outPath .. "meas/data")
		
		-- get oldest solution
		oldestSol = solTimeSeries:oldest()
		
		-- copy values into oldest solution (we reuse the memory here)
		VecScaleAssign(oldestSol, 1.0, u)
		
		-- push oldest solutions with new values to front, oldest sol pointer is popped from end
		solTimeSeries:push_discard_oldest(oldestSol, time)
		
		print("++++++ POINT IN TIME  " .. math.floor(time/dt+0.5)*dt .. "s  END ++++++++")
	end

end

