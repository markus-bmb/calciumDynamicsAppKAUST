------------------------------------------------------------------
-- Examination of calcium wave prerequisites                    --
--                                                              --
-- This script is intended to be used for simulations on 2d     --
-- representations of perfectly rotationally symmetric model    --
-- dendrites. It is supposed to be called by the script         --
-- 'wave_exam_batch.lua' during a bisection to find thresholds  --
-- for ER radius and RyR channel density above which a calcium  --
-- wave is elicited.                                            --
-- The LIMEX method is used for time stepping.                  --
--                                                              --
-- Author: mbreit                                               --
-- Date:   2017-08-29                                           --
------------------------------------------------------------------

-- for profiler output
--SetOutputProfileStats(true)

-- load pre-implemented lua functions
ug_load_script("ug_util.lua")
ug_load_script("plugins/Limex/limex_util.lua")

AssertPluginsLoaded({"neuro_collection", "Limex"})

EnableLUA2C(true)  -- speed up evaluation of lua functions by c program
--SetDebugLevel(debugID.LUACompiler, 0) 


-------------------------------------
-- parse command line parameters  ---
-------------------------------------
-- choice of grid name
gridName = util.GetParam("-grid", "modelDendrite1d.ugx")

-- grid parameters
dendLength = util.GetParamNumber("-dendLength", 50.0)
dendRadius = util.GetParamNumber("-dendRadius", dendRadius or 0.5)
erRadius = util.GetParamNumber("-erRadius", erRadius or 0.158)
nSeg = util.GetParamNumber("-nSeg", 96)

-- refinements (global and at ERM)
numRefs = util.GetParamNumber("-numRefs", 0)

-- which ER mechanisms are to be activated?
setting = util.GetParam("-setting", "ryr")
setting = string.lower(setting)
validSettings = {}
validSettings["all"] = 0;
validSettings["none"] = 0;
validSettings["ip3r"] = 0;
validSettings["ryr"] = 0;
if (validSettings[setting] == nil) then
    error("Unknown setting " .. setting)
end

-- densities
ryrDens = util.GetParamNumber("-ryrDens", ryrDens or 0.86)

-- buffer
totalBuffer = util.GetParamNumber("-totBuf", 4*40.0e-6)

-- whether to scale synaptic influx with dendritic radius
scaledInflux = util.HasParamOption("-scaledInflux")

-- choice of algebra
useBlockAlgebra = util.HasParamOption("-block")

-- choice of solver setup
solverID = util.GetParam("-solver", "GMG")
solverID = string.upper(solverID)
validSolverIDs = {}
validSolverIDs["GMG"] = 0
validSolverIDs["GS"] = 0
validSolverIDs["ILU"] = 0
if (validSolverIDs[solverID] == nil) then
    error("Unknown solver identifier " .. solverID)
end

-- error tolerance for Limex iteration
toleratedError = util.GetParamNumber("-tol", 0.01)
nstages = util.GetParamNumber("-nst", 3)

-- specify "-verbose" to output linear solver convergence
verbose = util.HasParamOption("-verbose")

-- parameters for instationary simulation
dt = util.GetParamNumber("-dt", 1e-2)
endTime = util.GetParamNumber("-endTime", 1.0)

-- choose outfile directory
outDir = util.GetParam("-outName", "caWaveExploration1d")
outDir = outDir .. "/"

-- specify -vtk to generate vtk output
generateVTKoutput = util.HasParamOption("-vtk")
pstep = util.GetParamNumber("-pstep", dt, "plotting interval")



-- init with dimension and algebra
if useBlockAlgebra then
	InitUG(1, AlgebraType("CPU", 7))
else
	InitUG(1, AlgebraType("CPU", 1))
end

-----------------------
-- geometry creation --
-----------------------
if ProcRank() == 0 then
	gen = DendriteGenerator()
	gen:set_dendrite_length(dendLength)
	gen:set_dendrite_radius(dendRadius)
	gen:set_er_radius(erRadius)
	gen:set_num_segments(nSeg)
	
	gridName = outDir .. "grid/" .. gridName
	gen:create_dendrite_1d(gridName)
end


-------------------------
--  problem constants  --
-------------------------
-- setting-dependent variables
withIP3R = true
withRyR = true
withSERCAandLeak = true

if setting == "none" then 
	withIP3R = false
	withRyR = false
	withSERCAandLeak = false
end

if setting == "ip3r" then
	withRyR = false
end

if setting == "ryr" then
	withIP3R = false
end


-- total cytosolic calbindin concentration
-- (four times the real value in order to simulate four binding sites in one)
totalClb = totalBuffer

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
ip3_init = 4.0e-8
clb_init = totalClb / (k_bind_clb/k_unbind_clb*ca_cyt_init + 1)


-- IP3 constants
reactionRateIP3 = 0.11
equilibriumIP3 = 4.0e-08
reactionTermIP3 = -reactionRateIP3 * equilibriumIP3

-- ER densities
IP3Rdensity = 17.3
RYRdensity = ryrDens --0.86
local v_s = 6.5e-27  -- V_S param of SERCA pump
local k_s = 1.8e-7   -- K_S param of SERCA pump
local j_ip3r = 3.7606194166520605e-23   -- single channel IP3R flux (mol/s) - to be determined via gdb
local j_ryr = 1.1201015633466695e-21    -- single channel RyR flux (mol/s) - to be determined via gdb
				  						-- ryr1: 1.1204582669024472e-21	
---[[
-- equilibration using SERCA
leakERconstant = 3.8e-17
local j_leak = ca_er_init-ca_cyt_init	-- leak proportionality factor


SERCAfluxDensity = leakERconstant * j_leak
if withIP3R then 
	SERCAfluxDensity = SERCAfluxDensity + IP3Rdensity * j_ip3r
end
if withRyR then
	SERCAfluxDensity = SERCAfluxDensity + RYRdensity * j_ryr
end
SERCAdensity = SERCAfluxDensity / (v_s/(k_s/ca_cyt_init+1.0)/ca_er_init)
if (SERCAdensity < 0) then error("SERCA flux density is outward for these density settings!") end
--]]
--[[
-- equilibration using leakage
SERCAdensity = 1973.0
SERCAflux = v_s / (k_s / ca_cyt_init + 1.0) / ca_er_init

netEquilFlux = SERCAdensity*SERCAflux
if withIP3R then 
	netEquilFlux = netEquilFlux - IP3Rdensity * j_ip3r
end
if withRyR then
	netEquilFlux = netEquilFlux - RYRdensity * j_ryr
end

leakERconstant = netEquilFlux / (ca_er_init - ca_cyt_init)
if (leakERconstant < 0) then
	error("ER leakage flux density is outward for these density settings!")
end
--]]

-- PM densities
pmcaDensity = 500.0
ncxDensity  = 15.0
vdccDensity = 0.0  -- 1.0
leakPMconstant =  pmcaDensity * 6.9672131147540994e-24	-- single pump PMCA flux (mol/s)
				+ ncxDensity *  6.7567567567567566e-23	-- single pump NCX flux (mol/s)
				--+ vdccDensity * (-1.5752042094823713e-25)    -- single channel VGCC flux (mol/s)
				-- *1.5 // * 0.5 for L-type // T-type
if (leakPMconstant < 0) then error("PM leak flux is outward for these density settings!") end


volScaleER = math.pi * erRadius*erRadius
volScaleCyt = math.pi * dendRadius*dendRadius - volScaleER

-- activation pattern
caEntryDuration = 0.001
function synCurrentDensityCa(x, t, si)
	-- single spike (~1200 ions)
	local influx = 0.0
	if t <= caEntryDuration then
		influx = 2.5e-3 * (1.0 - t/caEntryDuration)
	end
	
    return - volScaleCyt * influx
end

ip3EntryDelay = 0.000
ip3EntryDuration = 0.2
function synCurrentDensityIP3(x, t, si)
	local influx = 0.0
	if t > ip3EntryDelay and t <= ip3EntryDelay+ip3EntryDuration then
		influx = 7.5e-5 * (1.0 - t/ip3EntryDuration)
	end
	
    return - volScaleCyt * influx
end


-------------------------------
-- setup approximation space --
-------------------------------
-- load domain
reqSubsets = {"dend", "act", "meas"}
dom = util.CreateDomain(gridName, numRefs, reqSubsets)

-- create approximation space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("ca_cyt", "Lagrange", 1)
approxSpace:add_fct("ca_er", "Lagrange", 1)
approxSpace:add_fct("clb", "Lagrange", 1)
if withIP3R then
	approxSpace:add_fct("ip3", "Lagrange", 1)
end
if withRyR then
	approxSpace:add_fct("o2", "Lagrange", 1)
	approxSpace:add_fct("c1", "Lagrange", 1)
	approxSpace:add_fct("c2", "Lagrange", 1)
end
approxSpace:init_levels()
approxSpace:init_surfaces()
approxSpace:init_top_surface()
approxSpace:print_layout_statistic()
approxSpace:print_statistic()

OrderCuthillMcKee(approxSpace, true)


print(dom:domain_info():to_string())
--SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), outDir .. "grid/refined_grid_hierarchy_p" .. ProcRank() .. ".ugx", 1.0)
--SaveParallelGridLayout(dom:grid(), outDir .. "grid/parallel_grid_layout_p"..ProcRank()..".ugx", 1.0)


--------------------------
-- setup discretization --
--------------------------
-- diffusion --
diffCaCyt = ConvectionDiffusion("ca_cyt", "dend", "fv1")
diffCaCyt:set_mass_scale(volScaleCyt)
diffCaCyt:set_diffusion(D_cac*volScaleCyt)

diffCaER = ConvectionDiffusion("ca_er", "dend", "fv1")
diffCaER:set_mass_scale(volScaleER)
diffCaER:set_diffusion(D_cae*volScaleER)

diffClb = ConvectionDiffusion("clb", "dend", "fv1")
diffClb:set_mass_scale(volScaleCyt)
diffClb:set_diffusion(D_clb*volScaleCyt)

if withIP3R then
	diffIP3 = ConvectionDiffusion("ip3", "dend", "fv1")
	diffIP3:set_mass_scale(volScaleCyt)
	diffIP3:set_diffusion(D_ip3*volScaleCyt)
	diffIP3:set_reaction_rate(reactionRateIP3*volScaleCyt)
	diffIP3:set_reaction(reactionTermIP3*volScaleCyt)
end

-- buffering --
discBuffer = BufferFV1("dend") -- where buffering occurs
discBuffer:add_reaction(
	"clb",                     -- the buffering substance
	"ca_cyt",                  -- the buffered substance
	totalClb,                  -- total amount of buffer
	k_bind_clb*volScaleCyt,    -- binding rate constant
	k_unbind_clb*volScaleCyt   -- unbinding rate constant
)


-- er membrane transport systems
if withIP3R then
	ip3r = IP3R({"ca_cyt", "ca_er", "ip3"})
	ip3r:set_scale_inputs({1e3,1e3,1e3})
	ip3r:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)
	
	discIP3R = MembraneTransport1d("dend", ip3r)
	discIP3R:set_density_function(IP3Rdensity)
	discIP3R:set_radius(erRadius)
end

if withRyR then
	ryr = RyRImplicit({"ca_cyt", "ca_er", "o2", "c1", "c2"}, {"dend"})
	ryr:set_scale_inputs({1e3, 1e3, 1.0, 1.0, 1.0})
	ryr:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)
	ryrStateDisc = RyRImplicit_1drotsym({"ca_cyt", "ca_er", "o2", "c1", "c2"}, {"dend"})
	ryrStateDisc:set_calcium_scale(1e3)
	
	discRyR = MembraneTransport1d("dend", ryr)
	discRyR:set_density_function(RYRdensity)
	discRyR:set_radius(erRadius)
end

if withSERCAandLeak then
	serca = SERCA({"ca_cyt", "ca_er"})
	serca:set_scale_inputs({1e3,1e3})
	serca:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

	discSERCA = MembraneTransport1d("dend", serca)
	discSERCA:set_density_function(SERCAdensity)
	discSERCA:set_radius(erRadius)
	
	leakER = Leak({"ca_er", "ca_cyt"})
	leakER:set_scale_inputs({1e3,1e3})
	leakER:set_scale_fluxes({1e3}) -- from mol/(m^2 s) to (mol um)/(dm^3 s)
	
	discERLeak = MembraneTransport1d("dend", leakER)
	discERLeak:set_density_function(1e12*leakERconstant/(1e3)) -- from mol/(um^2 s M) to m/s
	discERLeak:set_radius(erRadius)
end


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


discPMCA = MembraneTransport1d("dend", pmca)
discPMCA:set_density_function(pmcaDensity)
discPMCA:set_radius(dendRadius)

discNCX = MembraneTransport1d("dend", ncx)
discNCX:set_density_function(ncxDensity)
discNCX:set_radius(dendRadius)

discPMLeak = MembraneTransport1d("dend", leakPM)
discPMLeak:set_density_function(1e12*leakPMconstant / (1.0-1e3*ca_cyt_init))
discPMLeak:set_radius(dendRadius)


-- synaptic activity
synapseInfluxCa = NeumannBoundary("ca_cyt", "fv1")
synapseInfluxCa:add("synCurrentDensityCa", "act", "dend")
if withIP3R then
	synapseInfluxIP3 = NeumannBoundary("ip3", "fv1")
	synapseInfluxIP3:add("synCurrentDensityIP3", "act", "dend")
end

-- domain discretization --
domDisc = DomainDiscretization(approxSpace)

domDisc:add(diffCaCyt)
domDisc:add(diffCaER)
domDisc:add(diffClb)

domDisc:add(discBuffer)

if withIP3R then
	domDisc:add(diffIP3)
	domDisc:add(discIP3R)
	domDisc:add(synapseInfluxIP3)
end
if withRyR then
	domDisc:add(discRyR)
	domDisc:add(ryrStateDisc) -- also add ryr as elem disc (for state variables)
end
if withSERCAandLeak then
	domDisc:add(discSERCA)
	domDisc:add(discERLeak)
end

domDisc:add(discPMCA)
domDisc:add(discNCX)
domDisc:add(discPMLeak)

domDisc:add(synapseInfluxCa)



-- setup time discretization --
timeDisc = ThetaTimeStep(domDisc)
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
dbgWriter:set_base_dir(outDir)
dbgWriter:set_vtk_output(false)

-- biCGstab --
convCheck = ConvCheck()
convCheck:set_minimum_defect(1e-50)
convCheck:set_reduction(1e-8)
convCheck:set_verbose(verbose)

if (solverID == "ILU") then
    bcgs_steps = 1000
    ilu = ILU()
    ilu:set_sort(true)
    bcgs_precond = ilu
elseif (solverID == "GS") then
    bcgs_steps = 1000
    bcgs_precond = GaussSeidel()
else -- (solverID == "GMG")
	gmg = GeometricMultiGrid(approxSpace)
	gmg:set_discretization(timeDisc)
	gmg:set_base_level(0)
	gmg:set_gathered_base_solver_if_ambiguous(true)
	
	-- treat SuperLU problems with Dirichlet constraints by using constrained version
	gmg:set_base_solver(SuperLU())
	
	smoother = GaussSeidel()
	gmg:set_smoother(smoother)
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

-- Newton solver
newtonSolver = LimexNewtonSolver()
newtonSolver:set_linear_solver(bicgstabSolver)
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
InterpolateInner(clb_init, u, "clb", 0.0)
if withIP3R then
	InterpolateInner(ip3_init, u, "ip3", 0.0)
end
if withRyR then
	ryrStateDisc:calculate_steady_state(u)
end

-- timestep in seconds
dtmin = 1e-9
dtmax = 1e-2
time = 0.0
step = 0

-- initial vtk output
if generateVTKoutput then
	out = VTKOutput()
	out:print(outDir .. "vtk/solution", u, step, time)
end


------------------
--  LIMEX setup --
------------------
stageNSteps = {}    -- number of time steps for each stage
for i = 1, nstages do stageNSteps[i] = i end

limex = LimexTimeIntegrator(nstages)
for i = 1, nstages do
	limex:add_stage(stageNSteps[i], newtonSolver, domDisc)
end

limex:set_tolerance(toleratedError)
limex:set_time_step(dt)
limex:set_dt_min(dtmin)
limex:set_dt_max(dtmax)
limex:set_increase_factor(2.0)
limex:set_reduction_factor(0.1)
limex:set_stepsize_greedy_order_factor(1)
limex:set_stepsize_safety_factor(0.25)

-- GridFunction error estimator (relative norm)
limexEstimator = ScaledGridFunctionEstimator()

errorEvalCaCyt = H1ComponentSpace("ca_cyt", "dend", 3)  -- fct names, subset names, order
errorEvalCaER = H1ComponentSpace("ca_er", "dend", 3)
errorEvalClb = H1ComponentSpace("clb", "dend", 3)
limexEstimator:add(errorEvalCaCyt)
limexEstimator:add(errorEvalCaER)
limexEstimator:add(errorEvalClb)

if withIP3R then
	errorEvalIP3 = H1ComponentSpace("ip3", "dend", 3)
	limexEstimator:add(errorEvalIP3)
end
if withRyR then
	errorEvalO2 = L2ComponentSpace("o2", 3, 1.0, "dend")  -- fct names, order, weight, subset names
	errorEvalC1 = L2ComponentSpace("c1", 3, 1.0, "dend")
	errorEvalC2 = L2ComponentSpace("c2", 3, 1.0, "dend")
	limexEstimator:add(errorEvalO2)
	limexEstimator:add(errorEvalC1)
	limexEstimator:add(errorEvalC2)
end

limex:add_error_estimator(limexEstimator)


-- for vtk output
if generateVTKoutput then 
	local vtkObserver = VTKOutputObserver(outDir .."vtk/solution", out, pstep)
	limex:attach_observer(vtkObserver)
end

-- for measurements
waveHitEnd = false
waveGotStuck = false
maxRyRFluxDens = 0.0
stuckWaveXPos = 0.0
interruptTime = 0.0


-- prepare output
waveFrontPosFile = outDir.."meas/waveFrontX.dat"
waveFrontPosFH = assert(io.open(waveFrontPosFile, "a"))

measInterval = 1e-4
lastMeasPt = -1
wpe = WaveProfileExporter(approxSpace, "ca_cyt", "dend", outDir .. "meas/waveProfile")


function measWaveActivity(step, time, dt)
	curSol = measObserver:get_current_solution()

	-- measure free Ca in dendrite
	take_measurement(curSol, time, "dend", "ca_cyt", outDir.."meas/meanCaCyt.dat")
	
	-- measure concentration at right end
	local measConc = take_measurement(curSol, time, "meas", "ca_cyt", outDir.."meas/caAtRightEnd_erRad"..erRadius.."_ryrDens"..ryrDens)
	if measConc > 4*ca_cyt_init then
		waveHitEnd = true
		interruptTime = time
		limex:interrupt()
	end
	
	-- measure wave front x position
	if withRyR then
		local measPt = math.floor(time/measInterval)
		waveFrontXPos = wave_front_x(curSol, "c1, c2", "dend", 0.1)
		if measPt > lastMeasPt then
			-- write current wave front location to file
			if ProcRank() == 0 then
				if waveFrontXPos < -1e10 then
					waveFrontPosFH:write(time, "\t", 0.0, "\n")
				else
					waveFrontPosFH:write(time, "\t", waveFrontXPos + 25.0, "\n")
				end
				waveFrontPosFH:flush()
			end
			
			-- write complete wave profile to file
			wpe:exportWaveProfileX(curSol, time)
		
			lastMeasPt = measPt	
		end
		
		
		-- measure maximal RyR flux density
		local ryrFluxDens = max_ryr_flux_density(curSol, "ca_cyt, ca_er, c1, c2", "dend", ryr)
		if ryrFluxDens > maxRyRFluxDens then
			maxRyRFluxDens = ryrFluxDens
		elseif ryrFluxDens < 0.25 * maxRyRFluxDens then -- wave got stuck
			waveGotStuck = true
			stuckWaveXPos = waveFrontXPos
			interruptTime = time
			limex:interrupt()
		end
	end
		
	print("Current (real) time: " .. time .. ",   last dt: " .. dt)
		
	return 0.0
end

measObserver = LuaCallbackObserver()
measObserver:set_callback("measWaveActivity")
limex:attach_observer(measObserver)

	
-- solve problem
limex:apply(u, endTime, u, time)

-- output outcome
print("")
if waveHitEnd then
	avgVelocity = dendLength / interruptTime / 1000.0
	print("#######################################")
	print("## A calcium wave has been detected. ##")
	print(string.format("## Average velocity was %5.3f um/ms. ##", avgVelocity))
	print("#######################################")
elseif waveGotStuck then
	interruptTimeMs = interruptTime * 1000.0
	print("##################################")
	print("## A calcium wave terminated at ##")
	print(string.format("## t = %5.2f ms, x = %5.2f um.  ##", interruptTimeMs, stuckWaveXPos+0.5*dendLength))
	print("##################################")
else
	print("###########################################")
	print("## A calcium wave has been elicited, but ##")
	print("## neither terminated nor hit the end.   ##")
	print("###########################################")
end


-- close output files
waveFrontPosFH:close()


if generateVTKoutput then 
	out:write_time_pvd(outDir .. "vtk/solution", u)
end


if doProfiling then
	WriteProfileData(outDir .."pd.pdxml")
end
