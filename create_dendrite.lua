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

--void BuildDendrite(number cyt_radius = 2.0, number er_radius = 0.5, number dend_length = 10.0, number pos_app = 5.0, number app_neck_radius = 0.4, number app_neck_length = 1.0, number app_head_radius = 0.3, number app_head_length = 0.3,
--		number spine_neck_radius = 1.0, number spine_neck_length = 0.5, number spine_head_radius = 0.5, number spine_head_length = 1.5, string fileName = "dendrite"){

--for var=1,10,1 do
--    BuildDendrite(0.25, 0.08, 6.0, 3, 0.02, 0.78, 0.041, 0.18, 0.075, 0.667, 0.185, 0.52, "vlachos" .. var .. ".ugx")
--end

--    BuildDendrite(0.3, 0.08, 6.0, 3, 0.02, 0.04, 0.041, 0.18, 0.125, 0.1, 0.185, 0.52, "vlachos_neu1.ugx")
--    BuildDendrite(0.29, 0.078, 6.0, 3, 0.02, 0.1, 0.041, 0.18, 0.125, 0.14, 0.185, 0.52, "vlachos_neu2.ugx")
--    BuildDendrite(0.28, 0.076, 6.0, 3, 0.02, 0.16, 0.041, 0.18, 0.125, 0.18, 0.185, 0.52, "vlachos_neu3.ugx")
--    BuildDendrite(0.27, 0.074, 6.0, 3, 0.02, 0.22, 0.041, 0.18, 0.125, 0.22, 0.185, 0.52, "vlachos_neu4.ugx")
--    BuildDendrite(0.26, 0.072, 6.0, 3, 0.02, 0.28, 0.041, 0.18, 0.125, 0.26, 0.185, 0.52, "vlachos_neu5.ugx")
--    BuildDendrite(0.25, 0.07, 6.0, 3, 0.02, 0.33, 0.041, 0.18, 0.125, 0.3, 0.185, 0.52, "vlachos_neu6.ugx")
--    BuildDendrite(0.24, 0.068, 6.0, 3, 0.02, 0.38, 0.041, 0.18, 0.125, 0.34, 0.185, 0.52, "vlachos_neu7.ugx")
--    BuildDendrite(0.23, 0.066, 6.0, 3, 0.02, 0.44, 0.041, 0.18, 0.125, 0.38, 0.185, 0.52, "vlachos_neu8.ugx")
--    BuildDendrite(0.22, 0.064, 6.0, 3, 0.02, 0.5, 0.041, 0.18, 0.125, 0.44, 0.185, 0.52, "vlachos_neu9.ugx")
--    BuildDendrite(0.21, 0.062, 6.0, 3, 0.04, 0.56, 0.041, 0.18, 0.125, 0.7, 0.225, 0.60, "bigheadwideapp1.ugx")
--    BuildDendrite(0.21, 0.062, 6.0, 3, 0.04, 0.56, 0.045, 0.18, 0.125, 0.7, 0.225, 0.60, "bigheadwideapp2.ugx")
--    BuildDendrite(0.21, 0.062, 6.0, 3, 0.04, 0.56, 0.049, 0.18, 0.125, 0.7, 0.225, 0.60, "bigheadwideapp3.ugx")
--    BuildDendrite(0.21, 0.062, 6.0, 3, 0.04, 0.56, 0.053, 0.18, 0.125, 0.7, 0.225, 0.60, "bigheadwideapp4.ugx")
--    BuildDendrite(0.45, 0.11, 10.0, 5, 0.01, 0.01, 0.0, 0.0, 0.074, 0.667, 0.186, 0.519, "paperERproto1.ugx")
--    BuildDendrite(0.45, 0.11, 10.0, 5.0, 0.028, 0.94, 0.03, 0.12, 0.074, 0.667, 0.186, 0.519, "paperSAproto1.ugx")
    BuildDendrite_Vector({0.45, 0.11, 10.0, 5.0, 0.028, 0.94, 0.03, 0.12, 0.074, 0.667, 0.186, 0.519}, {true, true, true}, "dendrite.ugx")
    
--    configurator = BDC("dendrite.ugx")
--    configurator.set_cyt_rad(2.0)
--	configurator.set_er_rad(0.5)
--	configurator.set_dend_len(10.0)
--	configurator.set_pos(5.0)
--	configurator.set_app_neck_rad(0.4)
--	configurator.set_app_neck_len(1.0)
--	configurator.set_app_head_rad(0.3)
--	configurator.set_app_head_len(0.3)
--	configurator.set_spine_neck_rad(1.0)
--	configurator.set_spine_neck_len(0.5)
--	configurator.set_spine_head_rad(0.5)
--	configurator.set_spine_head_len(1.5)
--    build_dendrite_with_conf(configurator)
    
    

