------------------------------------------------------------
-- Control script for Ashley's 1d/3d simulations          --
--                                                        --
-- author: mbreit                                         --
-- date:   2017-05-22                                     --
------------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("ashley_util.lua")


---------------------------------------------
---------------------------------------------
-- The following parameters need to be set --
-- before any simulation.                  --
---------------------------------------------
-- cell ID name to simulate on
cellID = util.GetParam("-cell", "SCR_140702_tracing")

-- stimulation layer (1: basal o., 2: apical rm., 3: apical lm.) 
layer = util.GetParamNumber("-layer", 1)

-- number of synapses to distribute within stimulation ball region
nSyn = util.GetParamNumber("-nsyn", 50)

-- peak conductance for synapses
peakCond = util.GetParamNumber("-cond", 1.2e-09)

-- output directory
outPath = util.GetParam("-outPath", cellID .. "_peakCond" .. peakCond)
----------------------------------------------
----------------------------------------------


-- parse ball regions and soma name for this cell
layers = parse_params(cellID, "coordinatesForModeling_final.txt")


-- simulation setup
simulation_setup = {
	-- cell ID
	cellName = cellID,
	
	-- number of synapses
	numSyn = nSyn,
	
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
			peak_conductance = peakCond,
			peak_cond_dev = 0.0
		},

		-- stimulation region center (placeholder only)
		region = {
			x = 0.0,
			y = 0.0,
			z = 0.0,
			w = 5e-5, -- diameter of each layer is 50 µm
		}
	},
  
	-- folder for paraview output (placeholder only)
	outputPath = ""
}


-- only simulate if layer center coordinate is not "nottraced";
if not (layers[layer][1] == "nottraced" or
        layers[layer][2] == "nottraced" or
        layers[layer][3] == "nottraced")
then
	layerName = layers[layer][4]:gsub("%s+", "_")

	-- create subdirectory for stimulation region in outPath (only root proc)
	if ProcRank() == 0 then
		os.execute("mkdir -p " .. outPath .. "/" .. layerName)
		
		-- create subdirectories for simulation output
		os.execute("mkdir " .. outPath .. "/" .. layerName .. "/grid")
		os.execute("mkdir " .. outPath .. "/" .. layerName .. "/meas")
		os.execute("mkdir " .. outPath .. "/" .. layerName .. "/vtk")
	end
	
	-- append layer name to output path
	simulation_setup.outputPath = outPath .. "/" .. layerName
	
	-- store layer center
	simulation_setup.stimSyn.region.x = layers[layer][1] * 1e-6
	simulation_setup.stimSyn.region.y = layers[layer][2] * 1e-6
	simulation_setup.stimSyn.region.z = layers[layer][3] * 1e-6
	
	-- do the actual simulation
	dofile 'ashley_1d-3d.lua'
end
