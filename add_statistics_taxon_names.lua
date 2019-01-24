local file_input_path = nil
local file_input = nil
local taxonomy_directory = nil
local file_output_path = nil
local file_output = nil
local file_size = nil
local file_position = nil
local line = nil
local next_line = nil
local taxons_id = {}
local names = {}
local total_taxons_id = 0

file_input_path = arg[1]
taxonomy_directory = arg[2]
file_output_path = arg[3]

-- Store hits

file_input = io.open(file_input_path, "r")
file_size = file_input:seek("end")
file_input:seek("set", 0)
next_line = file_input:lines()
line = next_line()
if line ~= nil then
	if string.sub(line, 1, 1) == "#" then
		line = next_line()
	end
end
while line ~= nil do
	local rank = {}
	rank.name, rank.taxon_id, rank.abundance = string.match(line, "([%w-]+)\t([%d-]+)\t(%d+)")
	taxons_id[rank.taxon_id] = rank
	if rank.name == "-" then
		rank.taxon_name = "-"
	end
	total_taxons_id = total_taxons_id + 1
	line = next_line()
end
io.close(file_input)

-- Read names

print("Adding names to statistics.")
file_input = io.open(taxonomy_directory .. "/names.dmp", "r")
file_size = file_input:seek("end")
file_input:seek("set", 0)
next_line = file_input:lines()
line = next_line()
while line ~= nil and total_taxons_id > 1 do
	local taxon_id = nil
	local taxon_name = nil
	local class = nil
	taxon_id, taxon_name, class = string.match(line, "(%d+)\t|\t([%w%p ]+)\t|[%w%p%s]+|\t([%w%p ]+)")
	if taxon_id ~= nil then
		if taxons_id[taxon_id] ~= nil then
			if taxons_id[taxon_id].taxon_name == nil then
				taxons_id[taxon_id].taxon_name = taxon_name
				total_taxons_id = total_taxons_id - 1
			elseif class == "scientific name" then
				taxons_id[taxon_id].taxon_name = taxon_name
			end
		end
	end
	line = next_line()
end
io.close(file_input)

-- Write output

file_output = io.open(file_output_path, "w")
file_output:write("#rank\ttaxon id\ttaxon name\tabundance\n")
for taxon_id, rank in pairs(taxons_id) do
	file_output:write(rank.name .. "\t" ..
		rank.taxon_id .. "\t" ..
		rank.taxon_name .. "\t" ..
		rank.abundance .. "\n")
end
io.flush()
io.close(file_output)
