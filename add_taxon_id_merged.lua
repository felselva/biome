local file_input_path = nil
local file_input = nil
local taxonomy_directory = nil
local file_output_path = nil
local file_output = nil
local file_position = nil
local file_size = nil
local line = nil
local next_line = nil
local hits = {}
local taxons_id_merged = {}

file_input_path = arg[1]
taxonomy_directory = arg[2]
file_output_path = arg[3]
file_input = io.open(file_input_path, "r")
file_output = io.open(file_output_path, "w")
if file_input == nil or file_output == nil then
	abort[1] = true
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
	hit.sample, hit.otu_identification, hit.size, hit.accession, hit.coverage, hit.similarity, hit.taxon =
		string.match(line, "([%w%p]+)\t(%d+)\t(%d+)\t([%w%p]+)\t([%d%p]+)\t([%d%p]+)\t(%d+)")
	hits[#hits + 1] = hit
	taxons_id_merged[hit.taxon] = true
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

-- Read merged taxon IDs

file_input = io.open(taxonomy_directory .. "/merged.dmp", "r")
file_size = file_input:seek("end")
file_input:seek("set", 0)
next_line = file_input:lines()
line = next_line()
while line ~= nil do
	local taxon = nil
	local new_taxon = nil
	taxon, new_taxon = string.match(line, "(%d+)\t|\t(%d+)")
	if taxons_id_merged[taxon] == true then
		taxons_id_merged[taxon] = new_taxon
	end
	file_position = file_input:seek()
	if file_position % 10000 == 0 or file_position == file_size then
		io.write("\rReading input at " .. string.format("%0.2f", (file_position / file_size) * 100.0) .. "%.")
		io.flush()
	end
	line = next_line()
end
io.write("\n")
io.flush()
io.close(file_input)

-- Save to file output

file_output:write("#sample\totu\tsize\taccession\tcoverage\tsimilarity\ttaxon\n")
for hit_index, hit in ipairs(hits) do
	local taxon_id = nil
	if taxons_id_merged[hit.taxon] == nil or taxons_id_merged[hit.taxon] == true then
		taxon_id = hit.taxon
	else
		taxon_id = taxons_id_merged[hit.taxon]
	end
	file_output:write(hit.sample .. "\t" ..
			hit.otu_identification .. "\t" ..
			hit.size .. "\t" ..
			hit.accession .. "\t" ..
			hit.coverage .. "\t" ..
			hit.similarity .. "\t" ..
			taxon_id .. "\n")
	if hit_index % 10000 == 0 or hit_index == #hits then
		io.write("\rWriting output at " .. string.format("%0.2f", (hit_index / #hits) * 100.0) .. "%.")
		io.flush()
	end
end
io.write("\n")
io.flush()
io.close(file_output)
