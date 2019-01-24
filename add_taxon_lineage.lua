local file_name_input = nil
local file_input = nil
local taxonomy_directory = nil
local file_name_output = nil
local file_output = nil
local file_position = nil
local file_size = nil
local ranks_text = nil
local get_rank = nil
local line = nil
local next_line = nil
local hits = {}
local nodes = {}
local ranks = {}

file_name_input = arg[1]
taxonomy_directory = arg[2]
ranks_text = arg[3]
file_name_output = arg[4]
file_input = io.open(file_name_input, "r")
file_output = io.open(file_name_output, "w")
if file_input == nil or file_output == nil then
	abort[1] = true
end

-- Ranks of interest 

for rank in string.gmatch(ranks_text, "%a+") do
	ranks[#ranks + 1] = rank
end

-- Store hits

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
	local hit = {}
	hit.sample, hit.otu, hit.size, hit.accession, hit.coverage, hit.similarity, hit.taxon =
		string.match(line, "([%w%p]+)\t(%d+)\t(%d+)\t([%w%p]+)\t([%d%p]+)\t([%d%p]+)\t(%d+)")
	hits[#hits + 1] = hit
	line = next_line()
	file_position = file_input:seek()
	if file_position % 10000 == 0 or file_position == file_size then
		io.write("\rReading input at " .. string.format("%0.2f", (file_position / file_size) * 100.0) .. "%.")
		io.flush()
	end
end
io.write("\n")
io.flush()
io.close(file_input)

-- Read nodes

file_input = io.open(taxonomy_directory .. "/nodes.dmp", "r")
file_size = file_input:seek("end")
file_input:seek("set", 0)
file_input_size = file_input:seek("end")
file_input:seek("set", 0)
next_line = file_input:lines()
line = next_line()
while line ~= nil do
	local taxon = nil
	local rank = nil
	local parent = nil
	local node = {}
	taxon, parent, rank = string.match(line, "(%d+)\t|\t(%d+)\t|\t(%w+)")
	node.parent = parent
	node.rank = rank
	nodes[taxon] = node
	line = next_line()
	file_position = file_input:seek()
	if file_position % 10000 == 0 or file_position == file_size then
		io.write("\rReading input at " .. string.format("%0.2f", (file_position / file_size) * 100.0) .. "%.")
		io.flush()
	end
end
io.write("\n")
io.flush()
io.close(file_input)

-- Save to file output

get_rank = function (rank, taxon_id)
	local result_taxon_id = nil
	local current_taxon_id = taxon_id
	while result_taxon_id == nil and current_taxon_id ~= nil and current_taxon_id ~= "1" do
		if nodes[current_taxon_id] ~= nil then
			if nodes[current_taxon_id].rank == rank then
				result_taxon_id = current_taxon_id
			else
				current_taxon_id = nodes[current_taxon_id].parent
			end
		else
			current_taxon_id = nil
		end
	end
	return result_taxon_id
end

file_output:write("#sample\totu\tsize\taccession\tcoverage\tsimilarity\ttaxon")
for rank_index, rank in ipairs(ranks) do
	file_output:write("\trank(" .. rank .. ")")
end
file_output:write("\n")

for hit_index, hit in ipairs(hits) do
	file_output:write(hit.sample .. "\t" ..
		hit.otu .. "\t" ..
		hit.size .. "\t" ..
		hit.accession .. "\t" ..
		hit.coverage .. "\t" ..
		hit.similarity .. "\t" ..
		hit.taxon)
	for rank_index, rank in ipairs(ranks) do
		local parent = nil
		parent = get_rank(rank, hit.taxon)
		if parent ~= nil then
			file_output:write("\t" .. parent)
		else
			file_output:write("\t-")
		end
	end
	file_output:write("\n")
	if hit_index % 10000 == 0 or hit_index == #hits then
		io.write("\rWriting output at " .. string.format("%0.2f", (hit_index / #hits) * 100.0) .. "%.")
		io.flush()
	end
end
io.write("\n")
io.flush()
io.close(file_output)
