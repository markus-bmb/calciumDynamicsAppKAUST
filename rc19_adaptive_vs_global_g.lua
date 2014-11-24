----------------------------------------------------------------
--  Example script for simulation on 3d reconstructed neuron  --
--                                                            --
--  Author: Markus Breit                                      --
----------------------------------------------------------------

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
gridName = "rc19/rc19amp.ugx"
--gridName = "rc19/rc19_amp_singleDend.ugx"
--gridName = "rc19/rc19_amp_measZones.ugx"
--gridName = "rc19/rc19_amp_new.ugx"
--gridName = "rc19/rc19_amp.ugx"					-- deprecated
--gridName = "rc19/RC19amp_ug4_finished.ugx"		-- dead
--gridName = "simple_reticulum_3d.ugx"				-- for testing

-- refinements before distributing grid
numPreRefs = util.GetParamNumber("-numPreRefs", 0)

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
endTime = util.GetParamNumber("-endTime", -1)
if (endTime == -1)
then
	-- choose number of time steps
	nTimeSteps = util.GetParamNumber("-nTimeSteps", 1)
	endTime = nTimeSteps*timeStep
end

-- chose plotting interval
plotStep = util.GetParamNumber("-pstep", 0.01)

-- choose solver setup
solverID = util.GetParam("-solver", "GS")
solverID = string.upper(solverID)
validSolverIDs = {}
validSolverIDs["GMG-GS"] = 0
validSolverIDs["GMG-ILU"] = 0
validSolverIDs["GMG-LU"] = 0
validSolverIDs["GMG-SLU"] = 0
--validSolverIDs["AMG-LU"] = 0
validSolverIDs["GS"] = 0
validSolverIDs["ILU"] = 0
validSolverIDs["JAC"] = 0
if (validSolverIDs[solverID] == nil) then
    error("Unknown solver identifier " .. solverID)
end

-- choose order of synaptic activity
order = util.GetParamNumber("-synOrder", 0)

-- choose time between synaptic stimulations (in ms)
jump = util.GetParamNumber("-jumpTime", 5)

-- choose outfile directory
fileName = util.GetParam("-outName", "rc19test")
fileName = fileName.."/"

-- specify -vtk to generate vtk output
generateVTKoutput = util.HasParamOption("-vtk")


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
	dens = dens * dcf * 17.3
	-- cluster for branching points
	if (si==22 or si==23 or si==24) then dens = dens * 10 end 
	return dens
end

function RYRdensity(x,y,z,t,si)
	local dens = math.abs(dendLengthPos(x,y,z))
	-- fourth order polynomial, distance to soma
	dens = 1.5 -3.5*dens +9.1*math.pow(dens,2) -10.5*math.pow(dens,3) +4.3*math.pow(dens,4)
	dens = dens * dcf * 0.86; 
	return dens
end

-- this is a little bit more complicated, since it must be ensured that
-- the net flux for equilibrium concentrations is zero
-- MUST be adapted whenever any parameterization of ER flux mechanisms is changed!
function SERCAdensity(x,y,z,t,si)
	local v_s = 6.5e-27						-- V_S param of SERCA pump
	local k_s = 1.8e-7						-- K_S param of SERCA pump
	local j_ip3r = 2.7817352713488838e-23 -- 3.7606194166520605e-23	-- single channel IP3R flux (mol/s) - to be determined via gdb
	local j_ryr = 4.6047720062808216e-22 -- 1.1204582669024472e-21	-- single channel RyR flux (mol/s) - to be determined via gdb
	local j_leak = ca_er_init-ca_cyt_init	-- leak proportionality factor
	
	local dens =  IP3Rdensity(x,y,z,t,si) * j_ip3r
				+ RYRdensity(x,y,z,t,si) * j_ryr
				+ LEAKERconstant(x,y,z,t,si) * j_leak
	dens = dens / (v_s/(k_s/ca_cyt_init+1.0)/ca_er_init)
	return dens
end

function LEAKERconstant(x,y,z,t,si)
	return dcf*3.4e-17
end

function PMCAdensity(x,y,z,t,si)
	return 500.0
end

function NCXdensity(x,y,z,t,si)
	return 15.0
end

function VGCCdensity(x,y,z,t,si)
	return 5.0
end

function LEAKPMconstant(x,y,z,t,si)
	local j_pmca = - 5.230769230769231e-24 -- - 6.9672131147540994e-24	-- single pump PMCA flux (mol/s) - to be determined via gdb
	local j_ncx = - 5.4347826086956515e-23 -- - 6.7567567567567566e-23	-- single pump NCX flux (mol/s) - to be determined via gdb
	local j_vgcc = 1.5752042094823713e-25	-- single channel VGCC flux (mol/s) - to be determined via gdb
		-- *1.5 // * 0.5 for L-type // T-type
		--j_vgcc = j_vgcc*1.5;
	local flux =  PMCAdensity(x,y,z,t,si) * j_pmca
				+ NCXdensity(x,y,z,t,si) * j_ncx
				+ VGCCdensity(x,y,z,t,si) * j_vgcc

	if (-flux < 0) then error("PM leak flux is outward for these density settings!") end
	
	return -flux -- 6.85e-22
end



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
	if 	(si>=synStart and si<=synStop and syns[si]<t and t<=syns[si]+caEntryDuration)
	--then efflux = -5e-6 * 11.0/16.0*(1.0+5.0/((10.0*(t-syns["start"..si])+1)*(10.0*(t-syns["start"..si])+1)))
	then efflux = -2e-4
	else efflux = 0.0
	end
    return -efflux
end


-- burst of ip3 at active synapses (triangular, immediate)
ip3EntryDelay = 0.000
ip3EntryDuration = 2.0
function ourNeumannBndIP3(x, y, z, t, si)
	if 	(si>=synStart and si<=synStop and syns[si]+ip3EntryDelay<t and t<=syns[si]+ip3EntryDelay+ip3EntryDuration)
	then efflux = - 2.1e-5/1.188 * (1.0 - (t-syns[si])/ip3EntryDuration)
	else efflux = 0.0
	end
    return -efflux
end

-------------------------------
-- setup approximation space --
-------------------------------
-- create, load, refine and distribute domain
print("create, refine and distribute domain")
neededSubsets = {}
distributionMethod = "metisReweigh"
weightingFct = InterSubsetPartitionWeighting()
weightingFct:set_default_weights(1,1)
weightingFct:set_inter_subset_weight(0, 1, 1000)
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets, distributionMethod, nil, nil, nil, weightingFct)

--[[
dom = util.CreateDomain(gridName, 0, neededSubsets)

balancer.partitioner = "dynBisection"
balancer.staticProcHierarchy = true
balancer.firstDistLvl = 0
balancer.firstDistProcs = 16
balancer.redistSteps = 1
balancer.redistProcs = 32
balancer.RefineAndRebalanceDomain(dom, numRefs)

-- in parallel environments: use a load balancer to distribute the grid
balancer.ParseParameters()
balancer.PrintParameters()
loadBalancer = balancer.CreateLoadBalancer(dom)

print(dom:domain_info():to_string())

if load_balancer ~= nil then
	loadBalancer:print_quality_records()
end
--]]

---[[
--print("Saving domain grid and hierarchy.")
--SaveDomain(dom, "refined_grid_p" .. ProcRank() .. ".ugx")
SaveGridHierarchyTransformed(dom:grid(), "refined_grid_hierarchy_p" .. ProcRank() .. ".ugx", 20.0)
print("Saving parallel grid layout")
SaveParallelGridLayout(dom:grid(), "parallel_grid_layout_p"..ProcRank()..".ugx", 20.0)
--]]

-- create approximation space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)

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
synapses = ""
for i=1,8 do
	synapses = synapses .. ", syn" .. i
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

approxSpace:add_fct("ca_er", "Lagrange", 1, innerDomain)
approxSpace:add_fct("ca_cyt", "Lagrange", 1, outerDomain)
approxSpace:add_fct("ip3", "Lagrange", 1, outerDomain)
approxSpace:add_fct("clb", "Lagrange", 1, outerDomain)

approxSpace:init_levels()
approxSpace:print_layout_statistic()
approxSpace:print_statistic()

----------------------------
-- setup error estimators --
----------------------------
eeCaCyt = SideAndElemErrEstData(2, 2, cytVol..", "..nucVol)
eeCaER 	= SideAndElemErrEstData(2, 2, erVol)
eeIP3 	= SideAndElemErrEstData(2, 2, cytVol..", "..nucVol)
eeClb 	= SideAndElemErrEstData(2, 2, cytVol..", "..nucVol)

----------------------------------------------------------
-- setup FV convection-diffusion element discretization --
----------------------------------------------------------
print ("Setting up Assembling")

-- Note: No VelocityField and Reaction is set. The assembling assumes default
--       zero values for them

if dim == 2 then 
    upwind = NoUpwind2d()
elseif dim == 3 then 
    upwind = NoUpwind3d()
end

elemDiscER = ConvectionDiffusion("ca_er", erVol, "fv1") 
elemDiscER:set_diffusion(D_cae)
elemDiscER:set_upwind(upwind)

elemDiscCYT = ConvectionDiffusion("ca_cyt", cytVol..", "..nucVol, "fv1")
elemDiscCYT:set_diffusion(D_cac)
elemDiscCYT:set_upwind(upwind)

elemDiscIP3 = ConvectionDiffusion("ip3", cytVol..", "..nucVol, "fv1")
elemDiscIP3:set_diffusion(D_ip3)
elemDiscIP3:set_reaction_rate(reactionRateIP3)
elemDiscIP3:set_reaction(reactionTermIP3)
elemDiscIP3:set_upwind(upwind)

elemDiscClb = ConvectionDiffusion("clb", cytVol..", "..nucVol, "fv1")
elemDiscClb:set_diffusion(D_clb)
elemDiscClb:set_upwind(upwind)

-- error estimators
elemDiscER:set_error_estimator(eeCaER)
elemDiscCYT:set_error_estimator(eeCaCyt)
elemDiscIP3:set_error_estimator(eeIP3)
elemDiscClb:set_error_estimator(eeClb)

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

-- error estimator
eeBuffering = MultipleSideAndElemErrEstData()
eeBuffering:add(eeClb)
eeBuffering:add(eeCaCyt)
eeBuffering:set_consider_me(false)

elemDiscBuffering:set_error_estimator(eeBuffering)

----------------------------------------------------
-- setup inner boundary (channels on ER membrane) --
----------------------------------------------------

-- We pass the function needed to evaluate the flux function here.
-- The order, in which the discrete fcts are passed, is crucial!
innerDiscIP3R = TwoSidedIP3RFV1("ca_cyt, ca_er, ip3", erMem)
innerDiscIP3R:set_density_function("IP3Rdensity")

innerDiscRyR = TwoSidedRyRFV1("ca_cyt, ca_er", erMem)
innerDiscRyR:set_density_function("RYRdensity")

innerDiscSERCA = TwoSidedSERCAFV1("ca_cyt, ca_er", erMem)
innerDiscSERCA:set_density_function("SERCAdensity")

innerDiscLeak = TwoSidedERCalciumLeakFV1("ca_cyt, ca_er", erMem)
innerDiscLeak:set_density_function("LEAKERconstant")

-- error estimators
eeERM = MultipleSideAndElemErrEstData()
eeERM:add(eeCaCyt)
eeERM:add(eeCaER)
eeERM:add(eeIP3)
eeERM:set_consider_me(false)

innerDiscIP3R:set_error_estimator(eeERM)
innerDiscRyR:set_error_estimator(eeERM)
innerDiscSERCA:set_error_estimator(eeERM)
innerDiscLeak:set_error_estimator(eeERM)

------------------------------
-- setup Neumann boundaries --
------------------------------
-- synaptic activity
neumannDiscCA = UserFluxBoundaryFV1("ca_cyt", plMem)
neumannDiscCA:set_flux_function("ourNeumannBndCA")
neumannDiscIP3 = UserFluxBoundaryFV1("ip3", plMem)
neumannDiscIP3:set_flux_function("ourNeumannBndIP3")

-- error estimators
eeNeumannCA = MultipleSideAndElemErrEstData()
eeNeumannCA:add(eeCaCyt)
eeNeumannCA:set_consider_me(false)
neumannDiscCA:set_error_estimator(eeNeumannCA)
eeNeumannIP3 = MultipleSideAndElemErrEstData()
eeNeumannIP3:add(eeIP3)
eeNeumannIP3:set_consider_me(false)
neumannDiscIP3:set_error_estimator(eeNeumannIP3)


-- plasme membrane transport systems
neumannDiscPMCA = OneSidedPMCAFV1("ca_cyt", plMem)
neumannDiscPMCA:set_density_function("PMCAdensity")

neumannDiscNCX = OneSidedNCXFV1("ca_cyt", plMem)
neumannDiscNCX:set_density_function("NCXdensity")

neumannDiscLeak = OneSidedPMCalciumLeakFV1("ca_cyt", plMem)
neumannDiscLeak:set_density_function("LEAKPMconstant")
--[[
neumannDiscVGCC = OneSidedBorgGrahamFV1WithVM2UG("ca_cyt", plMem, approxSpace,
		"neuronRes/timestep".."_order".. 0 .."_jump"..string.format("%1.1f", 5.0).."_", "%.3f", ".dat", false)
neumannDiscVGCC:set_channel_type_L() --default, but to be sure
neumannDiscVGCC:set_density_function("VGCCdensity")
neumannDiscVGCC:init(0.0)

voltageFilesInterval = 0.001;
--]]
-- error estimators
eePM = MultipleSideAndElemErrEstData()
eePM:add(eeCaCyt)
eePM:set_consider_me(false)


neumannDiscPMCA:set_error_estimator(eePM)
neumannDiscNCX:set_error_estimator(eePM)
neumannDiscLeak:set_error_estimator(eePM)
--neumannDiscVGCC:set_error_estimator(eePM)

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
domainDisc:add(neumannDiscPMCA)
domainDisc:add(neumannDiscNCX)
domainDisc:add(neumannDiscLeak)
--domainDisc:add(neumannDiscVGCC)

-- ER flux
domainDisc:add(innerDiscIP3R)
domainDisc:add(innerDiscRyR)
domainDisc:add(innerDiscSERCA)
domainDisc:add(innerDiscLeak)

-- constraints for adatptivity
domainDisc:add(OneSideP1Constraints())

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
superLU = SuperLU()

-- geometric multi-grid --
-- base solver
baseConvCheck = ConvCheck()
baseConvCheck:set_maximum_steps(1000)
baseConvCheck:set_minimum_defect(1e-28)
baseConvCheck:set_reduction(1e-1)
baseConvCheck:set_verbose(false)

if (solverID == "GMG-LU") then
    base = exactSolver
elseif (solverID == "GMG-SLU") then
    base = superLU
else
    base = LinearSolver()
    base:set_convergence_check(baseConvCheck)
    if (solverID == "GMG-ILU" or solverID == "SLU") then
        base:set_preconditioner(ilu)
    else
        base:set_preconditioner(gs)
    end
end

gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(timeDisc)
gmg:set_base_level(0)
if (solverID == "GMG-LU" or solverID == "SLU") then
    gmg:set_gathered_base_solver_if_ambiguous(true)
end
gmg:set_base_solver(base)
if (solverID == "GMG-ILU") then
    gmg:set_smoother(ilu)
else
    gmg:set_smoother(gs)
end 
gmg:set_cycle_type(1)
gmg:set_num_presmooth(10)
gmg:set_num_postsmooth(3)
gmg:set_rap(true)

--[[
-- AMG --
amg = RSAMGPreconditioner()
amg:set_num_presmooth(2)
amg:set_num_postsmooth(2)
amg:set_cycle_type(1)
amg:set_presmoother(gs)
amg:set_postsmoother(gs)
amg:set_base_solver(base)
--amg:set_debug(u)
--]]

-- biCGstab --
convCheck = ConvCheck()
convCheck:set_minimum_defect(1e-24)
convCheck:set_reduction(1e-06)
convCheck:set_verbose(false)
bicgstabSolver = BiCGStab()
if (solverID == "ILU") then
    convCheck:set_maximum_steps(2000)
    bicgstabSolver:set_preconditioner(ilu)
elseif (solverID == "GS") then
    convCheck:set_maximum_steps(2000)
    bicgstabSolver:set_preconditioner(gs)
elseif (solverID == "JAC") then
    convCheck:set_maximum_steps(2000)
    bicgstabSolver:set_preconditioner(jac)
else
    convCheck:set_maximum_steps(500)
    bicgstabSolver:set_preconditioner(gmg)
end
bicgstabSolver:set_convergence_check(convCheck)

-----------------------
-- non linear solver --
-----------------------
-- convergence check
newtonConvCheck = CompositeConvCheck3dCPU1(approxSpace, 20, 1e-20, 1e-08)
newtonConvCheck:set_component_check("ip3", 1e-20, 1e-08)
newtonConvCheck:set_verbose(true)
newtonConvCheck:set_time_measurement(true)
newtonConvCheck:set_adaptive(true)

newtonLineSearch = StandardLineSearch()
newtonLineSearch:set_maximum_steps(5)
newtonLineSearch:set_accept_best(true)
newtonLineSearch:set_verbose(false)

-- Newton solver
newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(bicgstabSolver)
newtonSolver:set_convergence_check(newtonConvCheck)
--newtonSolver:set_line_search(newtonLineSearch)

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


if (generateVTKoutput) then
	out = VTKOutput()
	out:print(fileName .. "vtk/result", u, step, time)
end

-- refiner setup
refiner = GlobalDomainRefiner(dom)
hn_refiner = HangingNodeDomainRefiner(dom)

TOL = 1e-15
refineFrac = 0.01
coarseFrac = 0.9
maxLevel = 6
maxElem = 7e6

-- set indicators for refinement in space and time to 0
space_refine_ind = 0.0
time_refine_ind = 0.0


outRefinement = VTKOutput()
approxSpace_vtk = ApproximationSpace(dom)
approxSpace_vtk:add_fct("eta_squared", "piecewise-constant");
u_vtk = GridFunction(approxSpace_vtk)
out_error = VTKOutput()
out_error:clear_selection()
out_error:select_all(false)
out_error:select_element("eta_squared", "error")

--takeMeasurement(u, time, measZonesERM, "ca_cyt, ca_er, ip3, clb", fileName .. "meas/data")

-- create new grid function for old value
uOld = u:clone()

-- store grid function in vector of  old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)


newTime = true
min_dt = timeStep / math.pow(2,15)
cb_interval = 10
lv = startLv
levelUpDelay = (synStop-synStart)*jump/1000.0 + caEntryDuration;
cb_counter = {}
for i=0,startLv do cb_counter[i]=0 end
while endTime-time > 0.001*dt do
	if newTime then
		print("++++++ POINT IN TIME  " .. math.floor((time+dt)/dt+0.5)*dt .. "s  BEGIN ++++++")
		newTime = false
		numCoarsenOld = -1.0;
		n=0
	end
	
	-- setup time disc for old solutions and timestep
	timeDisc:prepare_step(solTimeSeries, dt)
	
	-- prepare Newton solver
	if newtonSolver:prepare(u) == false then print ("Newton solver failed at step "..step.."."); exit(); end 
	
	-- prepare BG channel state
	--[[ never update, so that always -72mV
	if (time+dt<0.2) then
		vm_time = math.floor((time+dt)/voltageFilesInterval)*voltageFilesInterval	-- truncate to last time that data exists for
		neumannDiscVGCC:update_potential(vm_time)
	end
	neumannDiscVGCC:update_gating(time+dt)
	--]]
	
	-- apply Newton solver
	newton_fail = false
	error_fail = false
	
	if newtonSolver:apply(u) == false then
		-- in case of Newton convergence failure:
		newton_fail = true
		print ("Newton solver failed at point in time " .. time .. " with time step " .. dt)
	else
		-- check error
		-- timeDisc:mark_error(new solution, refiner, tolerance for total squared l2 error,
		--					   refinement fraction of max error, coarsening fraction (-1) of min error,
		--					   max #refinements)
		timeDisc:calc_error(u, u_vtk)
		
		changedGrid = false
		
		-- refining
		timeDisc:mark_for_refinement(hn_refiner, TOL, refineFrac, maxLevel)
		if hn_refiner:num_marked_elements() > 0 then
			hn_refiner:clear_marks()
			error_fail = true
			print ("Error estimator is above required error.")
		end
	end
	
	if (error_fail or newton_fail)
	then
		-- GRID REFINEMENT ------------------------------------------------
		
		-- if too many elements: 
		numElemBeforeRefinement = dom:domain_info():num_elements()
		if (numElemBeforeRefinement > maxElem) then
			print ("Adaptive refinement failed - too many elements. Aborting.")
			print ("Failed at point in time " .. time .. ".")
			time = endTime
		else
			refiner:refine()
			
			numElemAfterRefinement = dom:domain_info():num_elements()
			
			-- error is invalid, since grid has changed
			timeDisc:invalidate_error()
			
			--SaveDomain(dom, "refined_grid_".. math.floor((time+dt)/dt+0.5)*dt .. "_" .. n ..".ugx")
			n = n+1
			
			-- solve again on adapted grid
			VecScaleAssign(u, 1.0, solTimeSeries:latest())
			
			print ("Retrying with refined grid...")
		end
			-------------------------------------------------------------------
	else
		-- NORMAL PROCEEDING ----------------------------------------------
			
		numCoarsenOld = -1.0;
		-- update new time
		time = solTimeSeries:time(0) + dt
		newTime = true
		
		-- plot solution (& error) every plotStep seconds
		if (generateVTKoutput) then
			if math.abs(time/plotStep - math.floor(time/plotStep+0.5)) < 1e-5 then
				out:print(fileName .. "vtk/result", u, math.floor(time/plotStep+0.5), time)
				out_error:print(fileName .. "vtk/error_estimator_"..n, u_vtk, math.floor(time/plotStep+0.5), time)
			end
		end
		
		-- take measurement in nucleus every timeStep seconds 
		--takeMeasurement(u, time, measZonesERM, "ca_cyt, ca_er, ip3, clb", fileName .. "meas/data")
		
		-- export solution of ca on mem_er
		--exportSolution(u, approxSpace, time, "mem_cyt", "ca_cyt", fileName .. "sol/sol");
		
		print("++++++ POINT IN TIME  " .. math.floor(time/dt+0.5)*dt .. "s  END ++++++++");
	
		-- get oldest solution
		oldestSol = solTimeSeries:oldest()
			
		-- copy values into oldest solution (we reuse the memory here)
		VecScaleAssign(oldestSol, 1.0, u)
		
		-- push oldest solutions with new values to front, oldest sol pointer is popped from end
		solTimeSeries:push_discard_oldest(oldestSol, time)
		-------------------------------------------------------------------
	end
end

-- output of load balancing quality statistics
if loadBalancer ~= nil then
	loadBalancer:print_quality_records()
end

out_error:write_time_pvd(fileName .. "vtk/error_estimator", u_vtk)

-- end timeseries, produce gathering file
if (generateVTKoutput) then out:write_time_pvd(fileName .. "vtk/result", u) end

--[[
-- check if profiler is available
if GetProfilerAvailable() == true then
    print("")
    -- get node
    pn = GetProfileNode("main")
--    pn2 = GetProfileNode("GMG_lmgc")
    -- check if node is valid
    if pn:is_valid() then
	    print(pn:call_tree(0.0))
	    print(pn:groups())
--        print(pn2:total_time_sorted())
    else
        print("main is not known to the profiler.")
    end
else
    print("Profiler not available.")
end
--]]
