------------------------------------------------------------
-- Control script for Ashley's 1d/3d simulations          --
--                                                        --
-- author: mbreit                                         --
-- date:   2017-05-22                                     --
------------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("ashley_util.lua")

-- get cell ID name to simulate on
cellID = util.GetParam("-cell", "SCR_15_9cell1_2rec20150506_tracing")

-- parse ball regions and soma name for this cell
layers = parse_params(cellname, "coordinatesForModeling_final.txt")

-- get output path from command line
-- it needs to contain at least the subdirectory "grid"
-- as well as the subdirectory "vtk" if vtk output is enabled
outPath = util.GetParam("-outName", cellID)


-- simulation setup
simulation_setup = {
	-- cell ID
	cellName = cellID,
	
	-- number of synapses
	numSyn = 100,
	
	-- stimulation of synapses
	stimSyn = {
		-- bi-exponential synapses
		biexp = {
			onset_time = 0.0,
			onset_time_dev = 0.0,
			tau1_mean = 7e-4,
			tau1_dev = 0.0,
			tau2_mean = 5e-3,
			tau2_dev = 0.0,
			peak_conductance = 1.2e-9,
			peak_cond_dev = 0.0
		},

		-- stimulation region (e.g. basal oriens, apical radiatum or apical lm)
		region = {
			x = -5.274e-05,
			y = -0.00029732,
			z = -5.8e-05,
			w = 50, -- diameter of each layer is 50 µm
		}
	},
  
	-- folder for paraview output and csv data write out
	outputPath = ""
}



-- main batch loop (loops layers)
for i = 1, #layers do
	-- only simulate if layer coordinates are not "nottraced";
	if not (layers[i][1] == "nottraced" or layers[i][2] == "nottraced" or layers[i][3] == "nottraced") then
		simulation_setup.outputPath = outPath .. "__" .. layers[i][4]:gsub("%s+", "_")
		
		simulation_setup.stimSyn.region.x = layers[i][1] * 1e-6
		simulation_setup.stimSyn.region.y = layers[i][2] * 1e-6
		simulation_setup.stimSyn.region.z = layers[i][3] * 1e-6
		simulation_setup.stimSyn.region.w = 50 * 1e-6
		
		dofile 'ashley-1d-3d.lua'
	end
end
