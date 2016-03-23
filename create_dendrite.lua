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
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.01, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var1.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.04, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var2.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.07, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var3.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.1, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var4.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.13, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var5.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.16, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var6.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.19, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var7.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.22, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var8.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.25, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var9.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.28, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var10.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.31, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var11.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.34, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var12.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.37, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var13.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.4, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var14.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.43, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var15.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.46, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var16.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.49, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var17.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.52, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var18.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.55, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var19.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.58, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var20.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.61, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var21.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.64, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var22.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.67, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var23.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.7, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var24.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.73, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var25.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.76, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var26.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.79, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var27.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.82, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var28.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.85, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var29.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.88, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var30.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.91, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var31.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.94, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var32.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 0.97, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var33.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 1.0, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var34.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 1.03, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var35.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 1.06, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var36.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 1.09, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var37.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 1.12, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var38.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 1.15, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var39.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 1.18, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var40.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 1.21, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var41.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 1.24, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var42.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 1.27, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var43.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 1.29, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var44.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 1.33, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var45.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 1.36, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var46.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 1.39, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var47.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 1.42, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var48.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 1.45, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var49.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 1.48, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var50.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.013, 1.09, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var51.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.016, 1.09, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var52.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.019, 1.09, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var53.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.022, 1.09, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var54.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.025, 1.09, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var55.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.028, 1.09, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var56.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.031, 1.09, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var57.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.034, 1.09, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var58.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.037, 1.09, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var59.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.040, 1.09, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var60.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.043, 1.09, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var61.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.046, 1.09, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var62.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.049, 1.09, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var63.ugx")
    BuildDendrite({0.45, 0.11, 10.0, 5.0, 0.052, 1.09, 0.0, 0.0, 0.07925, 0.716, 0.21025, 0.579}, {true, true, true, false}, "var64.ugx")
    
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
    
    

