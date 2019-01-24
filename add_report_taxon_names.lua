local file_input_path = nil
local file_input = nil
local taxonomy_directory = nil
local file_output_path = nil
local file_output = nil
local file_size = nil
local file_position = nil
local line = nil
local next_line = nil
local otus = {}
local ranks = {}
local taxons_id_name = {}

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
		-- Get the ranks
		for word in string.gmatch(line, "rank%((%w+)%)") do
			ranks[#ranks + 1] = word
		end
		line = next_line()
	end
end
while line ~= nil do
	local ranks_values = nil
	local otu = {}
	otu.sample, otu.identification, otu.size, otu.sequence, otu.accession, otu.coverage, otu.similarity, otu.taxon, otu.rank_resolved, ranks_values =
		string.match(line, "([%w%p]+)\t([%d]+)\t([%d-]+)\t(%w+)\t([%w%p-]+)\t([%d%p-]+)\t([%d%p-]+)\t([%d-]+)\t([%w%p-]+)\t(.+)")
	taxons_id_name[otu.taxon] = true
	otu.ranks_index_taxon_id = {}
	for word in string.gmatch(ranks_values, "[%d-%(%)|]+") do
		if string.find(word, "%(") == nil then
			if word ~= "-" then
				taxons_id_name[word] = true
			end
		else
			for subword_identification, subword_identification_count in string.gmatch(word, "([%d-]+)%((%d+)%)") do
				if subword_identification ~= "-" then
					taxons_id_name[subword_identification] = true
				end
			end
		end
		otu.ranks_index_taxon_id[#otu.ranks_index_taxon_id + 1] = word
	end
	otus[#otus + 1] = otu
	line = next_line()
	file_position = file_input:seek()
	if file_position % 100 == 0 or file_position == file_size then
		io.write("\rReading input at " .. string.format("%0.2f", (file_position / file_size) * 100.0) .. "%.")
		io.flush()
	end
end
io.write("\n")
io.flush()
io.close(file_input)

-- Read names

file_input = io.open(taxonomy_directory .. "/names.dmp", "r")
file_size = file_input:seek("end")
file_input:seek("set", 0)
next_line = file_input:lines()
line = next_line()
while line ~= nil do
	local taxon_id = nil
	local taxon_name = nil
	local taxon_class = nil
	taxon_id, taxon_name, taxon_class = string.match(line, "(%d+)\t|\t([%w%p ]+)\t|[%w%p%s]+|\t([%w%p ]+)")
	if taxon_id ~= nil then
		if taxons_id_name[taxon_id] ~= nil then
			if taxons_id_name[taxon_id] == true then
				taxons_id_name[taxon_id] = name
			elseif class == "scientific name" then
				taxons_id_name[taxon_id] = name
			end
		end
	end
	line = next_line()
	file_position = file_input:seek()
	if file_position % 10000 == 0 or file_position == file_size then
		io.write("\rReading names at " .. string.format("%0.2f", (file_position / file_size) * 100.0) .. "%.")
		io.flush()
	end
end
io.write("\n");
io.flush()
io.close(file_input)

-- Write output

file_output = io.open(file_output_path, "w")
file_output:write("#sample\totu\tsize\tsequence\taccession\tcoverage\tsimilarity\ttaxon\trank_resolved")
for rank_index, rank in ipairs(ranks) do
	file_output:write("\trank(" .. rank .. ")")
end
for rank_index, rank in ipairs(ranks) do
	file_output:write("\tname(" .. rank .. ")")
end
file_output:write("\n")
for otu_index, otu in ipairs(otus) do
	local lineage = ""
	file_output:write(otu.sample .. "\t" ..
		otu.identification .. "\t" ..
		otu.size .. "\t" ..
		otu.sequence .. "\t" ..
		otu.accession .. "\t" ..
		otu.coverage .. "\t" ..
		otu.similarity .. "\t" ..
		otu.taxon .. "\t" ..
		otu.rank_resolved)
	for rank_index, rank in ipairs(ranks) do
		file_output:write("\t" .. otu.ranks_index_taxon_id[rank_index])
	end
	for rank_index, rank in ipairs(ranks) do
		if string.find(otu.ranks_index_taxon_id[rank_index], "%(") == nil then
			if taxons_id_name[otu.ranks_index_taxon_id[rank_index]] ~= true and taxons_id_name[otu.ranks_index_taxon_id[rank_index]] ~= nil then
				file_output:write("\t" .. taxons_id_name[otu.ranks_index_taxon_id[rank_index]])
			else
				file_output:write("\t-")
			end
		else
			file_output:write("\t")
			for rank_name, rank_count in string.gmatch(otu.ranks_index_taxon_id[rank_index], "([%d-]+)%((%d+)%)") do
				if taxons_id_name[rank_name] ~= true and taxons_id_name[rank_name] ~= nil then
					file_output:write(taxons_id_name[rank_name] .. "(" .. rank_count .. ")|")
				else
					file_output:write("-(" .. rank_count .. ")|")
				end
			end
		end
	end
	file_output:write("\n")
	if otu_index % 200 == 0 or otu_index == #otus then
		io.write("\rWriting output at " .. string.format("%0.2f", (otu_index / #otus) * 100.0) .. "%.")
		io.flush()
	end
end
io.write("\n")
io.flush()
io.close(file_output)
