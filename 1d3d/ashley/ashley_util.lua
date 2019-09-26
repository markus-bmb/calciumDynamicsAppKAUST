------------------------------------------------------
-- Utility functions for Ashley's 1d/3d simulations --
-- \author stephan grein <grein@temple.edu>         --
------------------------------------------------------

-- parse activation ball region centers for a cell --
function parse_params(cellname, coordFile)
   cellname = cellname:gsub("%+", "%%+")
   cellname = cellname:gsub("%-", "%%-")

   local f1 = assert(io.open(coordFile, "r"))
 
   local line = true
   local xCoord = {}
   local yCoord = {}
   local zCoord = {}
   local soma = ""

   while line do
     line = f1:read("*line")
     if (line) then
        if (string.find(line, cellname)) then
         local x = f1:read("*line")
         local y = f1:read("*line")
         local z = f1:read("*line")

         local first = true
         for i in string.gmatch(x, "%S+") do
            if (not first) then table.insert(xCoord, i) end
            first = false
         end

         first = true
         for i in string.gmatch(y, "%S+") do
            if (not first) then table.insert(yCoord, i) end
            first = false
         end

         first = true
         for i in string.gmatch(z, "%S+") do
            if (not first) then table.insert(zCoord, i) end
            first = false
         end

       line = false
       end
     end
   end
   
  local layers = {}
  table.insert(layers, {xCoord[1], yCoord[1], zCoord[1], "basal oriens"})
  table.insert(layers, {xCoord[2], yCoord[2], zCoord[2], "apical radiatum"})
  table.insert(layers, {xCoord[3], yCoord[3], zCoord[3], "apical lm"})
  return layers
end
