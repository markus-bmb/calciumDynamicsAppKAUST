------------------------------------------------------------------
-- Examination of calcium wave prerequisites                    --
--                                                              --
-- This script performs a bisection on either ER radius or RyR  --
-- channel density at a given dendritic radius to find out      --
-- threshold values for the creation of a stable calcium wave.  --
-- The actual simulation is performed by the script             --
-- 'wave_exam.lua' which is called in each step of the          --
-- bisection.                                                   --
--                                                              --
-- Author: mbreit                                               --
-- Date:   2017-08-28                                           --
------------------------------------------------------------------

-- load pre-implemented lua functions
ug_load_script("ug_util.lua")

-- current working directory
simFile = FindFileInStandardPaths("calciumDynamics_app/wave/wave_exam.lua")

-- choice of bisection parameter
bisecParam = util.GetParam("-bisec", "erRad")

if bisecParam == "erRad" then
	dendRadius = util.GetParamNumber("-dendRadius", 0.5)
	erRadius = dendRadius / 16.0
	inc = erRadius	
	accuracy = 0.001
elseif bisecParam == "ryrDens" then
	ryrDens = 1.0
	inc = ryrDens
	accuracy = 0.001
else
	error("Illegal bisection parameter '" .. bisecParam .. "' chosen.\n"
		  .. "Valid values are 'erRad' and 'ryrDens'.")
end

-- phase 1: double parameter until wave is elicited
waveHitEnd = false -- global variable waveHitEnd will be set to true by simulation script
			       -- iff a wave has been detected
while true do
	-- execute simulation with this parameter setting
	print("Executing simulation script '" .. simFile .. "' ...")
	dofile(simFile)
	
	if not waveHitEnd then
		print("")
		print("----------------------------------")
		print("No wave has been elicited.")
		if bisecParam == "erRad" then
			if erRadius >= 0.95*dendRadius then
				print("ERROR: Even with an ER radius of >= 0.95*dendRadius, "
					  .. "no calcium wave could be elicited.")
				exit()
			end
			inc = erRadius
			erRadius = 2.0*erRadius
			if erRadius > 0.95*dendRadius then
				inc = inc - (erRadius - 0.95*dendRadius)
				erRadius = 0.95*dendRadius
			end
			print("Increasing erRadius to " .. erRadius .. ".")
		elseif bisecParam == "ryrDens" then
			inc = ryrDens
			ryrDens = 2.0*ryrDens
			print("Increasing ryrDensity to " .. ryrDens .. ".")
		end
		print("----------------------------------")
		print("")
	else
		break
	end
end

-- phase 2: bisection to find the exact threshold
inc = 0.5*inc
while inc > accuracy do
	-- adjust bisection parameter
	print("")
	print("----------------------------------")
	if waveHitEnd then
		print("A calcium wave has been elicited.")
		if bisecParam == "erRad" then
			erRadius = erRadius - inc
			print("Decreasing erRadius to " .. erRadius .. ".")
		elseif bisecParam == "ryrDens" then
			ryrDens = ryrDens - inc
			print("Decreasing ryrDensity to " .. ryrDens .. ".")
		end
	else
		print("No wave has been elicited.")
		if bisecParam == "erRad" then
			erRadius = erRadius + inc
			print("Increasing erRadius to " .. erRadius .. ".")
		elseif bisecParam == "ryrDens" then
			ryrDens = ryrDens + inc
			print("Increasing ryrDensity to " .. ryrDens .. ".")
		end
	end
	print("----------------------------------")
	print("")
		
	-- execute simulation with this parameter setting
	dofile(simFile)
	
	-- adjust step size
	inc = 0.5*inc
end


-- generate output
print("")
print("")
print("----------------------------------------------")
if bisecParam == "erRad" then
	if waveHitEnd then
		erRadius = erRadius - inc
	else
		erRadius = erRadius + inc
	end
	print("erRadius threshold: " .. erRadius .." +/- " .. inc)
elseif bisecParam == "ryrDens" then
	if waveHitEnd then
		ryrDens = ryrDens - inc
	else
		ryrDens = ryrDens + inc
	end
	print("ryrDensity threshold: " .. ryrDens .." +/- " .. inc)
end
print("----------------------------------------------")
print("")
print("")

