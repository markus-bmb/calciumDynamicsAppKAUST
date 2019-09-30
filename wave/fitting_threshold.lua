--------------------------------------------------------------------------------
-- Data fitting for the wave ER threshold values.                             --
-- The fitting function is given by                                           --
--     r(rho) = - c*rho + ((c*rho)^2 + 2bc*rho + R^2) [mode 1]                --
-- or by                                                                      --
--     r(rho) = (a*c + b*rho) / (rho - c)             [mode 2]                --
--                                                                            --
-- Author: Markus Breit                                                       --
-- Date:   2017-10-04                                                         --
--------------------------------------------------------------------------------

mode = 1

threshRad = {}

--[[
R = 0.1
threshRad[0.0625] = 0.089
threshRad[0.125] = 0.076
threshRad[0.25] = 0.061
threshRad[0.5] = 0.046
threshRad[0.75] = 0.038
threshRad[1.0] = 0.034
threshRad[1.25] = 0.032
threshRad[1.5] = 0.029
threshRad[1.75] = 0.027
threshRad[2.0] = 0.026
threshRad[2.5] = 0.024
threshRad[3.0] = 0.023
threshRad[3.5] = 0.021
threshRad[4.0] = 0.021
threshRad[4.5] = 0.020
threshRad[5.0] = 0.020
threshRad[5.5] = 0.018
threshRad[6.0] = 0.018

local x = {}
x[1] = 0.008
x[2] = 0.26
--]]

--[[
R = 0.2
--threshRad[0.0625] = 0.182
--threshRad[0.125] = 0.164
threshRad[0.25] = 0.141
threshRad[0.5] = 0.113
threshRad[0.75] = 0.096
threshRad[1.0] = 0.085
threshRad[1.25] = 0.076
threshRad[1.5] = 0.070
threshRad[1.75] = 0.065
threshRad[2.0] = 0.062
threshRad[2.5] = 0.055
threshRad[3.0] = 0.051
threshRad[3.5] = 0.046
threshRad[4.0] = 0.045
threshRad[4.5] = 0.041
threshRad[5.0] = 0.040
threshRad[5.5] = 0.038
threshRad[6.0] = 0.037

local x = {}
x[1] = 0.03
x[2] = 0.4
--]]

--[[
R = 0.3
--threshRad[0.0625] = 0.280
--threshRad[0.125] = 0.256
threshRad[0.25] = 0.224
threshRad[0.5] = 0.187
threshRad[0.75] = 0.162
threshRad[1.0] = 0.142
threshRad[1.25] = 0.127
threshRad[1.5] = 0.115
threshRad[1.75] = 0.106
threshRad[2.0] = 0.098
threshRad[2.5] = 0.086
threshRad[3.0] = 0.077
threshRad[3.5] = 0.070
threshRad[4.0] = 0.064
threshRad[4.5] = 0.059
threshRad[5.0] = 0.056
threshRad[5.5] = 0.052
threshRad[6.0] = 0.050

local x = {}
if mode == 1 then
	x[1] = 0.01
	x[2] = 0.2
elseif mode == 2 then
	x[1] = -0.95
	x[2] = 0.26
	x[3] = 0.014
end
--]]

--[[
R = 0.4
--threshRad[0.0625] = 0.379
--threshRad[0.125] = 0.351
threshRad[0.25] = 0.312
threshRad[0.5] = 0.265
threshRad[0.75] = 0.232
threshRad[1.0] = 0.204
threshRad[1.25] = 0.180
threshRad[1.5] = 0.162
threshRad[1.75] = 0.145
threshRad[2.0] = 0.132
threshRad[2.5] = 0.112
threshRad[3.0] = 0.096
threshRad[3.5] = 0.085
threshRad[4.0] = 0.076
threshRad[4.5] = 0.070
threshRad[5.0] = 0.063
threshRad[5.5] = 0.059
threshRad[6.0] = 0.055

local x = {}
if mode == 1 then
	x[1] = 0.013
	x[2] = 0.35
elseif mode == 2 then
	x[1] = -0.5
	x[2] = 0.3
	x[3] = 0.014
end
--]]

---[[
R = 0.5
--threshRad[0.0625] = 0.448
--threshRad[0.125] = 0.402
threshRad[0.25] = 0.402
threshRad[0.5] = 0.346
threshRad[0.75] = 0.304
threshRad[1.0] = 0.265
threshRad[1.25] = 0.231
threshRad[1.5] = 0.204
threshRad[1.75] = 0.179
threshRad[2.0] = 0.159
threshRad[2.5] = 0.128
threshRad[3.0] = 0.108
threshRad[3.5] = 0.093
threshRad[4.0] = 0.081
threshRad[4.5] = 0.073
threshRad[5.0] = 0.067
threshRad[5.5] = 0.062
threshRad[6.0] = 0.058

local x = {}
if mode == 1 then
	x[1] = 0.005
	x[2] = 0.38
elseif mode == 2 then
	x[1] = 4.65
	x[2] = 0.0103
	x[3] = 0.0608
end
--]]

--[[
R = 0.6
--threshRad[0.125] = 0.545
--threshRad[0.25] = 0.495
--threshRad[0.5] = 0.425
--threshRad[0.75] = 0.378
--threshRad[1.0] = 0.328
--threshRad[1.25] = 0.282
threshRad[1.5] = 0.240
threshRad[1.75] = 0.204
threshRad[2.0] = 0.176
threshRad[2.5] = 0.138
threshRad[3.0] = 0.113
threshRad[3.5] = 0.097
threshRad[4.0] = 0.084
threshRad[4.5] = 0.076
threshRad[5.0] = 0.069
threshRad[5.5] = 0.063
threshRad[6.0] = 0.058

local x = {}
if mode == 1 then
	x[1] = 0.005
	x[2] = 0.38
elseif mode == 2 then
	x[1] = 1.0
	x[2] = 0.007
	x[3] = 0.25
end
--]]

--[[
R = 0.7
--threshRad[0.125] = 0.643
--threshRad[0.25] = 0.588
--threshRad[0.5] = 0.509
--threshRad[0.75] = 0.452
--threshRad[1.0] = 0.389
--threshRad[1.25] = 0.325
threshRad[1.5] = 0.266
threshRad[1.75] = 0.219
threshRad[2.0] = 0.185
threshRad[2.5] = 0.142
threshRad[3.0] = 0.114
threshRad[3.5] = 0.098
threshRad[4.0] = 0.085
threshRad[4.5] = 0.076
threshRad[5.0] = 0.069
threshRad[5.5] = 0.062
threshRad[6.0] = 0.058

local x = {}
if mode == 1 then
	x[1] = 0.005
	x[2] = 0.38
elseif mode == 2 then
	x[1] = 0.7
	x[2] = 0.006
	x[3] = 0.4
end
--]]

--[[
R = 0.8
--threshRad[0.125] = 0.742
--threshRad[0.25] = 0.683
--threshRad[0.5] = 0.595
--threshRad[0.75] = 0.524
--threshRad[1.0] = 0.447
--threshRad[1.25] = 0.360
threshRad[1.5] = 0.284
threshRad[1.75] = 0.227
threshRad[2.0] = 0.190
threshRad[2.5] = 0.143
threshRad[3.0] = 0.115
threshRad[3.5] = 0.098
threshRad[4.0] = 0.085
threshRad[4.5] = 0.076
threshRad[5.0] = 0.068
threshRad[5.5] = 0.063
threshRad[6.0] = 0.059

local x = {}
if mode == 1 then
	x[1] = 0.005
	x[2] = 0.38
elseif mode == 2 then
	x[1] = 0.45
	x[2] = 0.011
	x[3] = 0.56
end
--]]

--[[
R = 0.9
--threshRad[0.125] = 0.840
--threshRad[0.25] = 0.778
--threshRad[0.5] = 0.680
--threshRad[0.75] = 0.598
--threshRad[1.0] = 0.501
--threshRad[1.25] = 0.386
threshRad[1.5] = 0.293
threshRad[1.75] = 0.231
threshRad[2.0] = 0.191
threshRad[2.5] = 0.143
threshRad[3.0] = 0.115
threshRad[3.5] = 0.098
threshRad[4.0] = 0.085
threshRad[4.5] = 0.076
threshRad[5.0] = 0.069
threshRad[5.5] = 0.062
threshRad[6.0] = 0.059

local x = {}
if mode == 1 then
	x[1] = 0.005
	x[2] = 0.38
elseif mode == 2 then
	x[1] = 0.37
	x[2] = 0.013
	x[3] = 0.63
end
--]]

--[[
R = 1.0
--threshRad[0.125] = 0.940
--threshRad[0.25] = 0.875
--threshRad[0.5] = 0.768
--threshRad[0.75] = 0.671
--threshRad[1.0] = 0.550
--threshRad[1.25] = 0.401
threshRad[1.5] = 0.296
threshRad[1.75] = 0.233
threshRad[2.0] = 0.192
threshRad[2.5] = 0.144
threshRad[3.0] = 0.116
threshRad[3.5] = 0.097
threshRad[4.0] = 0.085
threshRad[4.5] = 0.075
threshRad[5.0] = 0.069
threshRad[5.5] = 0.063
threshRad[6.0] = 0.058

local x = {}
if mode == 1 then
	x[1] = 0.005
	x[2] = 0.38
elseif mode == 2 then
	x[1] = 0.37
	x[2] = 0.013
	x[3] = 0.63
end
--]]




function numData(table)
	count = 0
	for i,j in pairs(table) do
		count = count + 1
	end
	return count
end


function fitting(rho, arglist)
	if mode == 1 then
		local b = arglist[1]
		local c = arglist[2]
		local R = R
		if #arglist > 2 then R = arglist[3] end
		return -c*rho + math.pow(c*c*rho*rho + b*rho + R*R, 0.5)
	end
	if mode == 2 then
		local a = arglist[1]
		local b = arglist[2]
		local c = arglist[3]
		return (a*c + b*rho) / (rho - c)
	end
end

function fittingDeriv(rho, arglist)
	if mode == 1 then
		local b = arglist[1]
		local c = arglist[2]
		local R = R
		if #arglist > 2 then R = arglist[3] end
		local tmp = math.pow(c*c*rho*rho + b*rho + R*R, -0.5)
		
		result = {}
		result[1] = 0.5 * tmp * rho
		result[2] = -rho + tmp * c*rho*rho
		result[3] = tmp * R
	
		return result
	end
	if mode == 2 then
		local a = arglist[1]
		local b = arglist[2]
		local c = arglist[3]
		
		result = {}
		result[1] = c / (rho - c)
		result[2] = rho / (rho - c)
		result[3] = (a+b)*rho / ((rho - c)*(rho - c))
	
		return result
	end
end

function fittingHess(rho, arglist)
	if mode == 1 then
		local b = arglist[1]
		local c = arglist[2]
		local R = R
		if #arglist > 2 then R = arglist[3] end
		local tmp = math.pow(c*c*rho*rho + b*rho + R*R, -1.5)
		local tmp2 = math.pow(c*c*rho*rho + b*rho + R*R, -0.5)
		
		result = {}
		
		if #arglist > 2 then
			result[1] = -0.25 * tmp * rho*rho
			result[2] = -0.5 * tmp * c*rho*rho*rho
			result[3] = -0.5 * tmp * rho*R
			result[4] = -0.5 * tmp * c*rho*rho*rho
			result[5] = -tmp * c*c*rho*rho*rho*rho + tmp2 * rho*rho
			result[6] = -tmp * c*rho*rho*R
			result[7] = -0.5 * tmp * rho*R
			result[8] = -tmp * c*rho*rho*R
			result[9] = -tmp * R*R + tmp2
		else
			result[1] = -0.25 * tmp * rho*rho
			result[2] = -0.5 * tmp * c*rho*rho*rho
			result[3] = -0.5 * tmp * c*rho*rho*rho
			result[4] = -tmp * c*c*rho*rho*rho*rho + tmp2 * rho*rho
		end
		
		return result
	end
	if mode == 2 then
		local a = arglist[1]
		local b = arglist[2]
		local c = arglist[3]
		
		result = {}
		result[1] = 0.0
		result[2] = 0.0
		result[3] = rho / ((rho - c)*(rho - c))
		result[4] = 0.0
		result[5] = 0.0
		result[6] = rho / ((rho - c)*(rho - c))
		result[7] = rho / ((rho - c)*(rho - c))
		result[8] = rho / ((rho - c)*(rho - c))
		result[9] = 2.0*(a+b)*rho / ((rho - c)*(rho - c)*(rho - c))
		
		return result
	end
end



function objFct(arglist)
	local sum = 0.0
	for rho, rad in pairs(threshRad) do
		local tmp = fitting(rho, arglist) - rad
		tmp = tmp*tmp
		sum = sum + tmp
	end

	return 0.5*sum
end

function objFctDeriv(arglist)
	local sum = {}
	for i = 1, #arglist do
		sum[i] = 0.0
	end
	
	for rho, rad in pairs(threshRad) do
		local fit = fitting(rho, arglist)
		local fitDeriv = fittingDeriv(rho, arglist)
		for i = 1, #arglist do
			sum[i] = sum[i] + (fit - rad) * fitDeriv[i]
		end
	end

	return sum
end

function objFctHess(arglist)
	local sum = {}
	local n = #arglist
	for i = 1, n do
		for j = 1, n do
			sum[(i-1)*n+j] = 0.0
		end
	end
	
	for rho, rad in pairs(threshRad) do
		local fit = fitting(rho, arglist)
		local fitDeriv = fittingDeriv(rho, arglist)
		local fitHess = fittingHess(rho, arglist)
		for i = 1, n do
			for j = 1, n do
				sum[(i-1)*n+j] = sum[(i-1)*n+j] + fitDeriv[i]*fitDeriv[j] +
				                 (fit - rad)*fitHess[(i-1)*n+j]
			end
		end
	end

	return sum
end

--[[
-- test derivs
rhoTest = 2.0
d = fittingDeriv(rhoTest, x)
dbh = (fitting(rhoTest, {x[1]+1e-8, x[2]}) - fitting(rhoTest, x)) * 1e8
dch = (fitting(rhoTest, {x[1], x[2]+1e-8}) - fitting(rhoTest, x)) * 1e8
print("|(db - dbh) / dbh| = " .. math.abs((d[1]-dbh) / dbh))
print("|(dc - dch) / dch| = " .. math.abs((d[2]-dch)/dch))

dd = fittingHess(rhoTest, x)
ddxbh1 = fittingDeriv(rhoTest, {x[1]+1e-8, x[2]})
ddxch1 = fittingDeriv(rhoTest, {x[1], x[2]+1e-8})
ddxxh2 = fittingDeriv(rhoTest, {x[1], x[2]})
ddbbh = (ddxbh1[1] - ddxxh2[1]) * 1e8
ddbch = (ddxch1[1] - ddxxh2[1]) * 1e8
ddcbh = (ddxbh1[2] - ddxxh2[2]) * 1e8
ddcch = (ddxch1[2] - ddxxh2[2]) * 1e8
print("|(ddbb - ddbbh) / ddbbh| = " .. math.abs((dd[1]-ddbbh) / ddbbh))
print("|(ddbc - ddbch) / ddbch| = " .. math.abs((dd[2]-ddbch)/ddbch))
print("|(ddbc - ddbch) / ddbch| = " .. math.abs((dd[3]-ddbch) / ddbch))
print("|(ddcc - ddcch) / ddcch| = " .. math.abs((dd[4]-ddcch)/ddcch))
--]]

function defNorm(def)
	local sum = 0.0
	for i = 1, #def do
		sum = sum + def[i]*def[i]
	end
	return math.sqrt(sum)
end

function output(it, x, defN)
	local n = #x
	if mode == 1 then
		if n <= 2 then
			print("step " .. it .. ": b = " .. x[1]/(2*x[2]) .. ",  c = " .. x[2] .. 
			      ",  ||def|| = " .. defN .. ",  obj = " .. objFct(x)/numData(threshRad))
		else
			print("step " .. it .. ": b = " .. x[1]/(2*x[2]) .. ",  c = " .. x[2] .. ",  R = " .. x[3] .. 
			      ",  ||def|| = " .. defN .. ",  msq = " .. objFct(x)/numData(threshRad))
		end
		return
	end
	
	if mode == 2 then
		print("step " .. it .. ": a = " .. x[1] .. ",  b = " .. x[2] .. ",  c = " .. x[3] .. 
			      ",  ||def|| = " .. defN .. ",  msq = " .. objFct(x)/numData(threshRad))
	end
end


-- Newton iteration
function newton(x)
	local n = #x

	local def = objFctDeriv(x)
	local defN = defNorm(def)
	local defN0 = defN
	
	output(0, x, defN)
	
	local it = 1
	local maxIt = 20
	local minDef = 1e-12
	local minRed = 1e-08
	
	while defN > minRed*defN0 and defN > minDef and it <= maxIt do
		local H = objFctHess(x)
		
		-- solve H*c = d (Gauss elimination)
		local c = def
		for i = 2, n do
			for j = i, n do
				local fac = H[(j-1)*n+i-1] / H[(i-2)*n+i-1]
				for k = i, n do
					H[(j-1)*n+k] = H[(j-1)*n+k] - fac * H[(i-2)*n+k]
				end
				c[j] = c[j] - fac * c[i-1]
			end
		end
		for i = n, 1, -1 do
			for j = i+1, n do
				c[i] = c[i] - H[(i-1)*n+j] * c[j]
			end
			c[i] = c[i] / H[(i-1)*n+i]
		end
		
		-- apply correction and update
		for i = 1, n do
			x[i] = x[i] - c[i]
		end
		
		def = objFctDeriv(x)
		defN = defNorm(def)
		
		output(it, x, defN)
		
		it = it + 1
	end
end

print("Fitting with two parameters b, c:")
newton(x)

---[[
print("\nFitting with three parameters b, c, R:")
x[3] = R
newton(x)
--]]
