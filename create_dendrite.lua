----------------------------------------------------------
--  Example script for simulation on 3d spine model		--
--														--
--  Author:	Markus Breit								--
--  Date:	21-08-2014									--
----------------------------------------------------------

-- load pre-implemented lua functions
ug_load_script("ug_util.lua")

-- dimension
dim = 3

-- initialize ug with dimension and algebra
InitUG(dim, AlgebraType("CPU", 1));
--EnableLUA2C(true)	-- speed up evaluation of lua functions by c program
--SetDebugLevel(debugID.LUACompiler, 0) 




BuildDendrite(1.0, 1.0/3.0, 20.0, 5, 0.25, 0.6, 0.45, 0.5, "testCircle.ugx")
