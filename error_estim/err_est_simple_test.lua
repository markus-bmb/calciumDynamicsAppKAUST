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

-- choose dimension and algebra
InitUG(3, AlgebraType("CPU", 1))


-- choice of grid
gridName = "calciumDynamics_app/grids/test3d.ugx"

-- total refinements
numRefs = util.GetParamNumber("-numRefs", 0)
numPreRefs = util.GetParamNumber("-numPreRefs", 0)
adaptive_steps = util.GetParamNumber("-adaptSteps", 8)
	
-- choose outfile directory
outPath = util.GetParam("-outName", "error_estimator")
outPath = outPath.."/"

-- specify -vtk to generate vtk output
generateVTKoutput = util.HasParamOption("-vtk")


-------------------------------
-- setup approximation space --
-------------------------------
-- create, load, refine and distribute domain
neededSubsets = {"Inner", "Boundary1"}
dom = util.CreateDomain(gridName, numPreRefs, neededSubsets)

balancer.partitioner = "dynBisection"
balancer.firstDistLvl = 0
balancer.redistProcs = 64
balancer.parallelElementThreshold = 8
balancer.maxLvlsWithoutRedist = 1
balancer.ParseParameters()
balancer.PrintParameters()

-- in parallel environments: use a load balancer to distribute the grid
loadBalancer = balancer.CreateLoadBalancer(dom)
balancer.RefineAndRebalanceDomain(dom, numRefs, loadBalancer)

print(dom:domain_info():to_string())

if load_balancer ~= nil then
	loadBalancer:print_quality_records()
end

--[[
SaveGridHierarchyTransformed(dom:grid(), outPath.."/refined_grid_hierarchy_p" .. ProcRank() .. ".ugx", 2.0)
SaveParallelGridLayout(dom:grid(), outPath.."/parallel_grid_layout_p"..ProcRank()..".ugx", 2.0)
--]]

-- create approximation space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:init_levels()
approxSpace:print_layout_statistic()
approxSpace:print_statistic()


----------------------------
-- setup error estimators --
----------------------------
ee = SideAndElemErrEstData(1, 1, "Inner")

eeMult = MultipleSideAndElemErrEstData(approxSpace)
eeMult:add(ee, "c")
eeMult:set_consider_me(false)


----------------------------------------------------------
-- setup FV convection-diffusion element discretization --
----------------------------------------------------------
laplaceDisc = ConvectionDiffusionFV1("c", "Inner")
laplaceDisc:set_diffusion(1.0)
laplaceDisc:set_source(-1.0)

laplaceDisc:set_error_estimator(ee)

--[[
-- Neumann bnd
neumannDisc0 = UserFluxBoundaryFV1("c", "BoundaryN")
neumannDisc0:set_flux_function(1)

neumannDisc0:set_error_estimator(eeMult)
neumannDisc1:set_error_estimator(eeMult)
--]]

---[[
-- Dirichlet node
diri = DirichletBoundary()
diri:add(1.0, "c", "Boundary1")
diri:set_error_estimator(eeMult)
--]]


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
op_stat = AssembledLinearOperator(domainDisc)
op_stat:init()


------------------------------------------------------
-- solver setup --------------------------------------
------------------------------------------------------
-- geometric multi-grid --
gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_base_level(0)
gmg:set_gathered_base_solver_if_ambiguous(true)
gmg:set_base_solver(SuperLU())
gmg:set_smoother(GaussSeidel())
gmg:set_cycle_type(1)
gmg:set_num_presmooth(2)
gmg:set_num_postsmooth(2)
gmg:set_rap(true)
--gmg:set_debug(dbgWriter)

-- biCGstab --
convCheck = ConvCheck()
convCheck:set_minimum_defect(1e-50)
convCheck:set_reduction(1e-12)
convCheck:set_verbose(true)
convCheck:set_maximum_steps(100)

linearSolver = BiCGStab()
linearSolver:set_preconditioner(gmg)
linearSolver:set_convergence_check(convCheck)
--linearSolver:set_debug(GridFunctionDebugWriter(approxSpace))


-------------
-- solving --
-------------
-- get grid function
u = GridFunction(approxSpace)
u:set(0.0)
b = u:clone()

out = VTKOutput()
out:select_all(true)

--- ADAPTIVE REFINEMENT SETUP ---
-- refiner setup
refiner = HangingNodeDomainRefiner(dom)
TOL = 5e-4			-- maximally tolerated overall error
maxLevel = 8		-- maximal number of adaptive refinement levels
refStrat = StdRefinementMarking(TOL, maxLevel)

-- for visualization of adaptive refinement:
-- the outError can be used to write vtu files containing the error indicator
-- eta_squared (named "error" in the vtu file then) on each element.
if generateVTKoutput then
	approxSpace_vtk = ApproximationSpace(dom)
	approxSpace_vtk:add_fct("eta_squared", "piecewise-constant");
	u_vtk = GridFunction(approxSpace_vtk)
	out_error = VTKOutput()
	out_error:clear_selection()
	out_error:select_all(false)
	out_error:select_element("eta_squared", "error")
	Interpolate(0.0, u_vtk, "eta_squared", 0.0)
end

n = 0
error_fail = true

while error_fail and n <= adaptive_steps do
	error_fail = false

	AssembleLinearOperatorRhsAndSolution(op_stat, u, b)		
	if not ApplyLinearSolver(op_stat, u, b, linearSolver) then
		print("Linear solver failed."); exit();
	end
	
	-- calculate the error indicators (for solution u; output errors in u_vtk)
	if generateVTKoutput then
		domainDisc:calc_error(u, u_vtk)
	else
		domainDisc:calc_error(u)
	end
	
	-- mark elements with high indicator for refinement
	domainDisc:mark_with_strategy(refiner, refStrat)
	if refiner:num_marked_elements() > 0 then
		error_fail = true
		print ("Error estimator is above required error.")
	end

	-- refine the marked elements
	if error_fail then 
		-- print error and solution to vtu files
		if generateVTKoutput then
			out_error:print(outPath .. "vtk/error_failed_"..n, u_vtk, 0, 0.0)
			out:print(outPath .."vtk/solution_failed_"..n, u)
		end
		
		--[[		
		SaveVectorForConnectionViewer(u_vtk, outPath.."err"..n..".vec")
		SaveVectorForConnectionViewer(u, outPath.."sol"..n..".vec")
		
		-- save failed grid hierarchy
		SaveGridHierarchyTransformed(dom:grid(), outPath .. "refined_grid_hierarchy_failed"..n.."_p" .. ProcRank() .. ".ugx", 30)
		SaveParallelGridLayout(dom:grid(), outPath .. "parallel_grid_layout_failed"..n.."_p"..ProcRank()..".ugx", 2)
		--]]
		
		-- refine
		if n ~= adaptive_steps then
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
			
			print ("Retrying with refined grid...")
		end
		
		n = n+1
	else
		if generateVTKoutput then
			out_error:print(outPath .. "vtk/error_success_"..n, u_vtk, 0, 0.0)
			out:print(outPath .."vtk/solution_success_"..n, u, 0, 0.0)
		end
	end
end


if loadBalancer ~= nil then
	loadBalancer:print_quality_records()
end
