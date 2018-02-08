------------------------------------------------------------------
-- Evaluating the L2 distance of two grid functions             --
--                                                              --
-- Author: mbreit                                               --
-- Date:   2018-02-08                                           --
------------------------------------------------------------------

-- load pre-implemented lua functions
ug_load_script("ug_util.lua")

-- choice of grid name
gridName = util.GetParam("-grid", "simpleGrid.ugx")

-- init with dimension and algebra
InitUG(2, AlgebraType("CPU", 1))

-- create domain
dom = util.CreateDomain(gridName, 0, {})

-- create approximation space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("u", "Lagrange", 1)
approxSpace:init_levels()
approxSpace:init_surfaces()
approxSpace:init_top_surface()

-- generate two grid functions on given grid and approximation space
u1 = GridFunction(approxSpace)
u2 = GridFunction(approxSpace)

-- two arbitrary functions
function myU1(x,y)
	return x*y
end
function myU2(x,y)
	return (x-1)*y
end

-- set grid function values to represent the previously defined functions on the grid
InterpolateInner("myU1", u1, "u", 0.0)
InterpolateInner("myU2", u2, "u", 0.0)

-- calculate L2 distance ||u1-u2||
-- param1: first grid function
-- param2: component of first grid function (we only have "u")
-- param3: second grid function
-- param4: component of second grid function (we only have "u")
-- param5 is quadrature order (we need 2 since we have element-wise linear functions)
l2err = L2Error(u1, "u", u2, "u", 2)

print("L2 distance is " .. l2err .. ".")

