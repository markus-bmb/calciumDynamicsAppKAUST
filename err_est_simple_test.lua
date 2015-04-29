----------------------------------------------------------------
--  Example script for error estimation	in 3D				  --
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
gridName = "test3d.ugx"

-- total refinements
numRefs = util.GetParamNumber("-numRefs",    0)
numPreRefs = util.GetParamNumber("-numPreRefs",    0)
adaptive_steps = util.GetParamNumber("-adaptSteps",    1)

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
nTimeSteps = util.GetParamNumber("-nTimeSteps", 1)
def_endTime = nTimeSteps*timeStep
endTime = util.GetParamNumber("-endTime", def_endTime)

-- choose plotting interval
plotStep = util.GetParamNumber("-pstep", timeStep)

-- choose solver setup
solverID = util.GetParam("-solver", "ILUC")
solverID = string.upper(solverID)
validSolverIDs = {}
validSolverIDs["GMG-GS"] = 0;
validSolverIDs["GMG-ILU"] = 0;
validSolverIDs["GMG-LU"] = 0;
validSolverIDs["GMG-SLU"] = 0;
validSolverIDs["JAC"] = 0;
validSolverIDs["GS"] = 0;
validSolverIDs["SLU"] = 0;
validSolverIDs["ILU"] = 0;
validSolverIDs["ILUC"] = 0;
validSolverIDs["ILUT"] = 0;
if (validSolverIDs[solverID] == nil) then
    error("Unknown solver identifier " .. solverID)
end
 
-- choose outfile directory
fileName = util.GetParam("-outName", "error_estimator")
fileName = fileName.."/"

-- specify -vtk to generate vtk output
generateVTKoutput = util.HasParamOption("-vtk")


-------------------------------
-- setup approximation space --
-------------------------------
-- create, load, refine and distribute domain
print("create, refine and distribute domain")
neededSubsets = {}
dom = util.CreateDomain(gridName, numPreRefs, neededSubsets)

balancer.partitioner = "dynBisection"
--balancer.firstDistLvl = 0
balancer.redistProcs = 64
balancer.parallelElementThreshold = 8
balancer.maxLvlsWithoutRedist = 1
balancer.ParseParameters()
balancer.PrintParameters()

write(">> distributing and refining domain ...\n")
-- in parallel environments: use a load balancer to distribute the grid
loadBalancer = balancer.CreateLoadBalancer(dom)
balancer.RefineAndRebalanceDomain(dom, numRefs, loadBalancer)
write(">> distributing done\n")

print(dom:domain_info():to_string())

if load_balancer ~= nil then
	loadBalancer:print_quality_records()
end

---[[
--print("Saving domain grid and hierarchy.")
SaveDomain(dom, fileName.."/refined_grid_p" .. ProcRank() .. ".ugx")
SaveGridHierarchyTransformed(dom:grid(), fileName.."/refined_grid_hierarchy_p" .. ProcRank() .. ".ugx", 2.0)
--print("Saving parallel grid layout")
SaveParallelGridLayout(dom:grid(), fileName.."/parallel_grid_layout_p"..ProcRank()..".ugx", 2.0)
--]]

-- create approximation space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)

approxSpace:add_fct("c", "Lagrange", 1)

approxSpace:init_levels()
approxSpace:print_layout_statistic()
approxSpace:print_statistic()

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

laplaceDisc = ConvectionDiffusion("c", "Inner", "fv1")
laplaceDisc:set_diffusion(1.0)
laplaceDisc:set_upwind(upwind)

ee = SideAndElemErrEstData(2, 2, "Inner")
laplaceDisc:set_error_estimator(ee)

--[[
-- Neumann bnd
neumannDisc0 = UserFluxBoundaryFV1("c", "Boundary0")
neumannDisc0:set_flux_function(1)
neumannDisc1 = UserFluxBoundaryFV1("c", "Boundary1")
neumannDisc1:set_flux_function(-0.2)

eeNeumann0 = MultipleSideAndElemErrEstData()
eeNeumann0:add(ee)
eeNeumann0:set_consider_me(false)
neumannDisc0:set_error_estimator(eeNeumann0)
eeNeumann1 = MultipleSideAndElemErrEstData()
eeNeumann1:add(ee)
eeNeumann1:set_consider_me(false)
neumannDisc1:set_error_estimator(eeNeumann1)
--]]

-- Dirichlet node
diri = DirichletBoundary()
diri:add(0, "c", "Boundary0")
diri:add(1, "c", "Boundary1")

------------------------------------------
-- setup complete domain discretization --
------------------------------------------
domainDisc = DomainDiscretization(approxSpace)

domainDisc:add(laplaceDisc)
domainDisc:add(diri)

-- constraints for adatptivity
adapt_constraint = OneSideP1Constraints()
domainDisc:add(adapt_constraint)

-- create stationary operator from domain discretization
op_stat = AssembledOperator(domainDisc)

-------------------------------
-- setup time discretization --
-------------------------------
timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(1.0) -- 1.0 is implicit Euler

-- create operator from discretization
op = AssembledOperator()
op:set_discretization(timeDisc)
op:init()

------------------------------------------------------
-- solver setup --------------------------------------
------------------------------------------------------

-- create algebraic preconditioner
jac = Jacobi()
jac:set_damp(0.8)
gs = GaussSeidel()
sgs = SymmetricGaussSeidel()
bgs = BackwardGaussSeidel()
ilu = ILU()
--iluc = ILUC(approxSpace)
--iluc:add_constraint(adapt_constraint)
ilut = ILUT(1e-4)

-- exact solver
exactSolver = LU()
superLU = SuperLU()


-- geometric multi-grid --
-- base solver
-- baseConvCheck = ConvCheck()
-- baseConvCheck:set_maximum_steps(1000)
-- baseConvCheck:set_minimum_defect(1e-30)
-- baseConvCheck:set_reduction(1e-8)
-- baseConvCheck:set_verbose(false)

-- -- debug writer
-- dbgWriter = GridFunctionDebugWriter(approxSpace)
-- dbgWriter:set_vtk_output(false)

-- if (solverID == "GMG-LU") then
--     base = exactSolver
-- elseif (solverID == "GMG-SLU") then
--     base = superLU
-- else
--     base = LinearSolver()
--     base:set_convergence_check(baseConvCheck)
--     if (solverID == "GMG-ILU") then
--         base:set_preconditioner(ilu)
--     else
--         base:set_preconditioner(gs)
--     end
-- end

-- gmg
-- gmg = GeometricMultiGrid(approxSpace)
-- gmg:set_discretization(timeDisc)
-- gmg:set_base_level(0)
-- if (solverID == "GMG-LU" or solverID == "GMG-SLU") then
--     gmg:set_gathered_base_solver_if_ambiguous(true)
-- end
-- gmg:set_base_solver(base)
-- gmg:set_smoother(iluc)
-- gmg:set_cycle_type(1)
-- gmg:set_num_presmooth(3)
-- gmg:set_num_postsmooth(3)
--gmg:set_rap(true) -- causes error in base solver!!
--gmg:set_debug(dbgWriter)

-- biCGstab --
convCheck = ConvCheck()
convCheck:set_minimum_defect(1e-50)
convCheck:set_reduction(1e-12)
convCheck:set_verbose(true)
bicgstabSolver = BiCGStab()
if (solverID == "ILU") then
    bcgs_steps = 10000
    bcgs_precond = ilu
elseif (solverID == "ILUC") then
    bcgs_steps = 10000
    bcgs_precond = iluc
elseif (solverID == "ILUT") then
    bcgs_steps = 10000
    bcgs_precond = ilut
elseif (solverID == "JAC") then
    bcgs_steps = 10000
    bcgs_precond = jac
elseif (solverID == "GS") then
    bcgs_steps = 10000
    bcgs_precond = gs
elseif (solverID == "SLU") then
    bcgs_steps = 1000
    bcgs_precond = superLU
else
    bcgs_steps = 1000
    bcgs_precond = gmg
end
convCheck:set_maximum_steps(bcgs_steps)
bicgstabSolver:set_preconditioner(bcgs_precond)
bicgstabSolver:set_convergence_check(convCheck)
--print(bicgstabSolver:config_string())
--bicgstabSolver:set_debug(GridFunctionDebugWriter(approxSpace))
--bicgstabSolver:set_restart(2000)

--- non-linear solver ---
-- convergence check
newtonConvCheck = CompositeConvCheck(approxSpace, 20, 5e-17, 1e-12)
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
u:set(0.0);


-- stationary problem

out = VTKOutput()
out:select_all(true)


--- ADAPTIVE REFINEMENT SETUP ---
-- refiner setup
refiner = HangingNodeDomainRefiner(dom)
TOL = 5e-2			-- maximally tolerated overall error
refineFrac = 0.01	-- without effect in the current implementation (I believe)
coarseFrac = 0.9	-- same here
maxLevel = 6		-- maximal number of adaptive refinement levels
maxElem = 1e10		-- maximal number of elements in the geometry

-- for visualization of adaptive refinement:
-- the outError can be used to write vtu files containing the error indicator
-- eta_squared (named "error" in the vtu file then) on each element.
approxSpace_vtk = ApproximationSpace(dom)
approxSpace_vtk:add_fct("eta_squared", "piecewise-constant");
u_vtk = GridFunction(approxSpace_vtk)
out_error = VTKOutput()
out_error:clear_selection()
out_error:select_all(false)
out_error:select_element("eta_squared", "error")
Interpolate(0.0, u_vtk, "eta_squared", 0.0)


n = 0
error_fail = true

while error_fail and n <= adaptive_steps do
	
	newtonSolver:init(op_stat)
	newtonSolver:prepare(u)
	
	error_fail = false
	
	if newtonSolver:apply(u) == false then
	       	print("Newton Solver failed."); exit();
	end
	
	-- calculate the error indicators (for solution u; output errors in u_vtk)
	domainDisc:calc_error(u, u_vtk)
	
	-- mark elements with high indicator for refinement
	domainDisc:mark_for_refinement(refiner, TOL, refineFrac, maxLevel)
	if refiner:num_marked_elements() > 0 then
		error_fail = true
		print ("Error estimator is above required error.")
	end

	-- refine the marked elements
	if error_fail then 
		-- print error and solution to vtu files
		out_error:print(fileName .. "vtk/error_failed_"..n, u_vtk, 0, 0.0)
		out:print(fileName .."vtk/solution_failed_"..n, u)

		SaveVectorForConnectionViewer(u_vtk, fileName.."err"..n..".vec")
		SaveVectorForConnectionViewer(u, fileName.."sol"..n..".vec")
		
		-- save failed grid hierarchy
		---[[
		SaveGridHierarchyTransformed(dom:grid(), fileName .. "refined_grid_hierarchy_failed"..n.."_p" .. ProcRank() .. ".ugx", 30)
		SaveParallelGridLayout(dom:grid(), fileName .. "parallel_grid_layout_failed"..n.."_p"..ProcRank()..".ugx", 2)
		--]]
		
		-- export solution
		--exportSolution(u, approxSpace, 0, "Inner", "c", fileName.."sol/exp_sol_failed_"..n)
		
		-- refine
		refiner:refine()
		
		-- rebalancing (should distribute the geometry evenly among the processors)
		if loadBalancer ~= nil then
			print("rebalancing...")
			balancer.Rebalance(dom, loadBalancer)
			loadBalancer:create_quality_record("time_".. n)
		end
		
		print(dom:domain_info():to_string())
		
		-- error is invalid, since grid has changed
		domainDisc:invalidate_error()
		
		n = n+1
		if n <= adaptive_steps then print ("Retrying with refined grid...") end	
	else
		out_error:print(fileName .. "vtk/error_success_"..n, u_vtk, 0, 0.0)
		out:print(fileName .."vtk/solution_success_"..n, u, 0, 0.0)
	end
end


if loadBalancer ~= nil then
	loadBalancer:print_quality_records()
end
