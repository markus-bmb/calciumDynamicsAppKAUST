-- This file converts xlsx data (V_m or fluorescence) from Anna
-- to files that can be read using Vm2Ug.

io.input("o_table.csv")

-- first line: time points
local first_line = io.read()
if first_line == nil then exit("File is empty!") end

local count = 0
time_points = {}
values = {}

for value in string.gmatch(first_line, '([^,]+)') do
	count = count+1
	--print("time point "..count..": "..value)
	time_points[count] = value
	values[count] = {}
end

--print("count: "..count)

-- read fluorescence values
local point = 0
while true do
	point = point+1
	--print("point:"..point)
	local line = io.read()
	if line == nil then break end
	
	-- values are comma-separated
	local time = 0
	for value in string.gmatch(line, '([^,]+)') do
    	time = time+1
   		--print("time:"..time..": "..value)
    	values[time][point] = value
	end
end

-- output files
for i=2,count do
	io.output(string.format("dFoF_%.4f.dat", 0.001*time_points[i]))
	
	for j=1,point-1 do
		io.write(string.format("%.4f %.4f %.4f \t%.4f", 0.0, 0.0, values[1][j], values[i][j]), "\n")
	end
	--print("output for time point "..i.." done")
end




