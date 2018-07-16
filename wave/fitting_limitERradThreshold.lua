------------------------------------------------------------------
-- Least squares fit for limit ER radius threshold              --
-- Fitting function is of the form                              --
--   r(rho) = (c + b*rho) / (rho - a).                          --
--                                                              --
-- Author: mbreit                                               --
-- Date:   2017-09-21                                           --
------------------------------------------------------------------

--[[
limitRad = {}
limitRad[1.0] = 0.686
limitRad[1.25] = 0.427
limitRad[1.5] = 0.298
limitRad[2.0] = 0.190
limitRad[2.5] = 0.142
limitRad[3.0] = 0.114
limitRad[3.5] = 0.097
limitRad[4.0] = 0.085
limitRad[4.5] = 0.075
limitRad[5.0] = 0.069
limitRad[5.5] = 0.063
limitRad[6.0] = 0.060
limitRad[10.0] = 0.038
limitRad[20.0] = 0.024

local a = 0.65
local b = 0.015
local c = 0.23
--]]

---[[
-- low buffer
limitRad = {}
limitRad[0.5] = 0.299
limitRad[1.0] = 0.157
limitRad[1.5] = 0.114
limitRad[2.0] = 0.091
limitRad[2.5] = 0.079
limitRad[3.0] = 0.069
limitRad[3.5] = 0.063
limitRad[4.0] = 0.058
limitRad[4.5] = 0.054
limitRad[5.0] = 0.050
limitRad[5.5] = 0.048
--limitRad[6.0] = 0.046
--limitRad[10.0] = 0.034
--limitRad[20.0] = 0.023

local a = 0.055
local b = 0.025
local c = 0.125
--]]

function fittingFct(a, b, c, rho)
	return (c + b*rho) / (rho - a)
end
function fittingFctDeriv(a, b, c, rho)
	return (c + b*rho) / ((rho-a)*(rho-a)), rho/(rho-a), 1.0 / (rho-a)
end
function fittingFctHess(a, b, c, rho)
	return 2.0*(c+b*rho)/((rho-a)*(rho-a)*(rho-a)), rho / ((rho-a)*(rho-a)), 1.0 / ((rho-a)*(rho-a)),
           rho / ((rho-a)*(rho-a)), 0.0, 0.0,
           1.0 / ((rho-a)*(rho-a)), 0.0, 0.0
end


function relErrorSq(a, b, c)
	local sum = 0.0
	for rho,rad in pairs(limitRad) do
		sum = sum + math.pow(1.0 - rad / fittingFct(a, b, c, rho), 2)  
	end
	
	return sum
end
function relErrorSqDeriv(a, b, c)
	local sumA = 0.0
	local sumB = 0.0
	local sumC = 0.0
		
	for rho,rad in pairs(limitRad) do
		local ffdA, ffdB, ffdC = fittingFctDeriv(a, b, c, rho)
		sumA = sumA + 2.0 * (1.0 - rad / fittingFct(a, b, c, rho)) *
				rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdA
		sumB = sumB + 2.0 * (1.0 - rad / fittingFct(a, b, c, rho)) *
				rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdB
		sumC = sumC + 2.0 * (1.0 - rad / fittingFct(a, b, c, rho)) *
				rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdC
	end
	
	return sumA, sumB, sumC
end
function relErrorSqHess(a, b, c)
	local sumAA = 0.0
	local sumAB = 0.0
	local sumAC = 0.0
	local sumBA = 0.0
	local sumBB = 0.0
	local sumBC = 0.0
	local sumCA = 0.0
	local sumCB = 0.0
	local sumCC = 0.0
		
	for rho,rad in pairs(limitRad) do
		local ffdA, ffdB, ffdC = fittingFctDeriv(a, b, c, rho)
		local ffdAA, ffdAB, ffdAC, ffdBA, ffdBB, ffdBC, ffdCA, ffdCB, ffdCC = fittingFctHess(a, b, c, rho)
		sumAA = sumAA + 2.0 * rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdA *
				rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdA
				+ 2.0 * (1.0 - rad / fittingFct(a, b, c, rho)) *
					(-2.0*rad / math.pow(fittingFct(a, b, c, rho), 3) * ffdA*ffdA
					 + rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdAA)
		sumAB = sumAB + 2.0 * rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdB *
				rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdA
				+ 2.0 * (1.0 - rad / fittingFct(a, b, c, rho)) *
					(-2.0*rad / math.pow(fittingFct(a, b, c, rho), 3) * ffdA*ffdB
					 + rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdAB)
		sumAC = sumAC + 2.0 * rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdC *
				rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdA
				+ 2.0 * (1.0 - rad / fittingFct(a, b, c, rho)) *
					(-2.0*rad / math.pow(fittingFct(a, b, c, rho), 3) * ffdA*ffdC
					 + rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdAC)
		sumBA = sumBA + 2.0 * rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdA *
				rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdB
				+ 2.0 * (1.0 - rad / fittingFct(a, b, c, rho)) *
					(-2.0*rad / math.pow(fittingFct(a, b, c, rho), 3) * ffdA*ffdB
					 + rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdBA)
		sumBB = sumBB + 2.0 * rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdB *
				rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdB
				+ 2.0 * (1.0 - rad / fittingFct(a, b, c, rho)) *
					(-2.0*rad / math.pow(fittingFct(a, b, c, rho), 3) * ffdB*ffdB
					 + rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdBB)
		sumBC = sumBC + 2.0 * rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdC *
				rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdB
				+ 2.0 * (1.0 - rad / fittingFct(a, b, c, rho)) *
					(-2.0*rad / math.pow(fittingFct(a, b, c, rho), 3) * ffdB*ffdC
					 + rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdBC)
		sumCA = sumCA + 2.0 * rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdA *
				rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdC
				+ 2.0 * (1.0 - rad / fittingFct(a, b, c, rho)) *
					(-2.0*rad / math.pow(fittingFct(a, b, c, rho), 3) * ffdA*ffdC
					 + rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdCA)
		sumCB = sumCB + 2.0 * rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdB *
				rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdC
				+ 2.0 * (1.0 - rad / fittingFct(a, b, c, rho)) *
					(-2.0*rad / math.pow(fittingFct(a, b, c, rho), 3) * ffdB*ffdC
					 + rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdCB)
		sumCC = sumCC + 2.0 * rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdC *
				rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdC
				+ 2.0 * (1.0 - rad / fittingFct(a, b, c, rho)) *
					(-2.0*rad / math.pow(fittingFct(a, b, c, rho), 3) * ffdC*ffdC
					 + rad / math.pow(fittingFct(a, b, c, rho), 2) * ffdCC)
	end
	
	return sumAA, sumAB, sumAC, sumBA, sumBB, sumBC, sumCA, sumCB, sumCC
end

local defa, defb, defc = relErrorSqDeriv(a, b, c)
local def = math.sqrt(defa*defa + defb*defb + defc*defc)
local def0 = def

print("Step 0: a = " .. a .. ", b = " .. b ..",  c = " .. c ..",  ||def|| = " .. def)

it = 0
while (def > 1e-12 and def > 1e-8*def0 and it < 10) do
	local haa, hab, hac, hba, hbb, hbc, hca, hcb, hcc = relErrorSqHess(a, b, c)
	
	-- solve H*c = d (Gauss elimination)
	corra = defa
	corrb = defb
	corrc = defc
	
	local fac = hba/haa
	hbb = hbb - hab*fac
	hbc = hbc - hac*fac
	corrb = corrb - corra*fac
	fac = hca/haa
	hcb = hcb - hab*fac
	hcc = hcc - hac*fac
	corrc = corrc - corra*fac
	fac = hcb/hbb
	hcc = hcc - hbc*fac
	corrc = corrc - corrb*fac
	
	corrc = corrc / hcc
	corrb = (corrb - corrc*hbc) / hbb
	corra = (corra - corrb*hab - corrc*hac) / haa

	a = a - corra
	b = b - corrb
	c = c - corrc
	
	defa, defb, defc = relErrorSqDeriv(a, b, c)
	def = math.sqrt(defa*defa + defb*defb + defc*defc)
	
	it = it + 1
	print("Step " .. it .. ": a = " .. a .. ", b = " .. b .. ", c = " .. c ..",  ||def|| = " .. def)
end