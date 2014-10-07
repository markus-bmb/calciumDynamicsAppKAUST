----------------------------------------------------------
--  Example script for simulation on 3d spine model		--
--														--
--  Author:	Marcus Kessler								--
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

--setSpineNeckRadius(0.2)

--void BuildDendrite(number cyt_radius, number er_radius, number dend_length, int pos_app, number app_neck_radius, number app_neck_length, number app_head_radius, number app_head_height, number spine_neck_radius, number spine_neck_length, number spine_head_radius, number spine_head_height, string fileName){
--BuildDendrite(1.0, 1.0/3.0, 20.0, 5, 0.25, 0.6, 0.1, 0.1, "testCircle.ugx")
BuildDendrite(2.0, 1.0, 20.0, 2, 0.4, 0.4, 0.1, 0.1, "testCircle.ugx")
