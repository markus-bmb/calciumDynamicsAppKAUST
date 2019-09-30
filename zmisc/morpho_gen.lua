--------------------------------------------------------------------------------
-- Test script for development of morphology generator                        --
--                                                                            --
-- Author: Markus Breit                                                       --
-- Date:   2016-08-09                                                         --
--------------------------------------------------------------------------------

-- load pre-implemented lua functions
ug_load_script("ug_util.lua")

-- choice of grid
gridName = util.GetParam("-grid", "../grids/modelDendrite.ugx")

dendLength = util.GetParamNumber("-dendLength", 50.0)
dendRadius = util.GetParamNumber("-dendRadius", 0.5)
erRadius = util.GetParamNumber("-erRadius", 0.158)
synArea = util.GetParamNumber("-synArea", math.pi*0.5)
nSeg = util.GetParamNumber("-nSeg", 100)
bobbel = util.HasParamOption("-bob")

gen = DendriteGenerator()
gen:set_dendrite_length(dendLength)
gen:set_dendrite_radius(dendRadius)
gen:set_er_radius(erRadius)
gen:set_synapse_area(synArea)
gen:set_num_segments(nSeg)
if bobbel then gen:set_bobbel_er(1, 1) end

gen:create_dendrite_middle_influx(gridName)
