local file_input_path = nil
local file_input = nil
local file_input_fasta_path = nil
local file_output_path = nil
local file_output = nil
local file_position = nil
local file_size = nil
local line = nil
local next_line = nil
local otus = {}
local ranks = {}
local otu_index = 0
local otus_total = 0
local check_otu = nil
local get_rank_index = nil
local statistics = {}

file_input_path = arg[1]
file_input_fasta_path = arg[2]
file_output_path = arg[3]

-- Store OTUs

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
			statistics[word] = {}
			statistics[word].names = {}
		end
		line = next_line()
	end
end
statistics["-"] = {}
statistics["-"].names = {}
while line ~= nil do
	local ranks_values = nil
	local otu = nil
	local hit = {}
	hit.sample, hit.otu_identification, hit.size, hit.accession, hit.coverage, hit.similarity, hit.taxon, ranks_values =
		string.match(line, "([%w%p]+)\t(%d+)\t(%d+)\t([%w%p]+)\t([%d%p]+)\t([%d%p]+)\t(%d+)(.+)")
	hit.ranks = {}
	for word in string.gmatch(ranks_values, "[%d-]+") do
		hit.ranks[#hit.ranks + 1] = word
	end
	otu = otus[hit.sample .. "_" .. hit.otu_identification]
	if otu == nil then
		otu = {}
		otu.hits = {}
		otu.identification = hit.otu_identification
		otu.sample = hit.sample
		otu.size = hit.size
		otu.ranks = {}
		otu.ranks_resolved = {}
		otu.ranks_best_candidate = {}
		otu.ranks_best_candidate_total = {}
		otus[hit.sample .. "_" .. hit.otu_identification] = otu
		otus_total = otus_total + 1
	end
	otu.hits[#otu.hits + 1] = hit
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
print("Total of " .. otus_total .. " OTUs.")

-- Read FASTA

file_input = io.open(file_input_fasta_path, "r")
file_size = file_input:seek("end")
file_input:seek("set", 0)
next_line = file_input:lines()
line = next_line()
while line ~= nil do
	local otu = nil
	local sample_name = nil
	local otu_name = nil
	local sequence = nil
	sample_name = string.match(line, "sample=([%w%.]+)")
	otu_name = string.match(line, "otu=(%d+)")
	sequence = next_line()
	if otus[sample_name .. "_" .. otu_name] ~= nil and sequence ~= nil then
		otus[sample_name .. "_" .. otu_name].sequence = sequence
	end
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
print("Total of " .. otus_total .. " OTUs.")

-- Check OTUs

get_rank_index = function(rank)
	local index = 0
	for rank_index, rank_name in ipairs(ranks) do
		if rank == rank_name then
			index = rank_index
		end
	end
	return index
end

check_resolved_rank = function(otu, rank_index)
	local rank_list = {}
	local rank_total = 0
	local best_hit = otu.hits[1]
	local contain_undefined_rank = false
	for hit_index, hit in ipairs(otu.hits) do
		if string.find(hit.ranks[rank_index], "-") ~= nil then
			contain_undefined_rank = true
		end
		if rank_list[hit.ranks[rank_index]] == nil then
			rank_list[hit.ranks[rank_index]] = 1
			rank_total = rank_total + 1
		else
			rank_list[hit.ranks[rank_index]] = rank_list[hit.ranks[rank_index]] + 1
		end
		if tonumber(hit.similarity) > tonumber(best_hit.similarity) then
			best_hit = hit
		elseif tonumber(hit.similarity) == tonumber(best_hit.similarity) then
			if tonumber(hit.coverage) > tonumber(best_hit.coverage) then
				best_hit = hit
			end
		end
	end
	otu.best_hit = best_hit
	otu.ranks[rank_index] = rank_list
	if rank_total == 1 and contain_undefined_rank == false then
		otu.resolved = true
		if otu.best_rank_index_resolved ~= nil then
			if rank_index < otu.best_rank_index_resolved then
				otu.best_rank_index_resolved = rank_index
			end
		else
			otu.best_rank_index_resolved = rank_index
		end
		otu.ranks_resolved[rank_index] = true
		otu.ranks_best_candidate[rank_index] = otu.hits[1].ranks[rank_index]
	else
		local rank_best_candidate_total = 0
		for otu_rank_name, otu_rank_count in pairs(otu.ranks[rank_index]) do
			if otu_rank_count > rank_best_candidate_total then
				rank_best_candidate_total = otu_rank_count
			end
		end
		otu.ranks_best_candidate_total[rank_index] = rank_best_candidate_total
	end
end

for sample_otu, otu in pairs(otus) do
	for rank_index, rank in ipairs(ranks) do
		check_resolved_rank(otu, rank_index)
	end
	-- Get rank resolved
	for rank_index, rank in ipairs(ranks) do
		if otu.rank_resolved == nil and otu.ranks_resolved[rank_index] == true then
			otu.rank_resolved = rank
			otu.rank_index_resolved = rank_index
		end
	end
end

-- Fill statistics

for sample_otu, otu in pairs(otus) do
	for rank_index, rank in ipairs(ranks) do
		check_resolved_rank(otu, rank_index)
	end
	-- Get rank resolved
	for rank_index, rank_name in ipairs(ranks) do
		if otu.rank_resolved == rank_name then
			if statistics[rank_name].names[otu.ranks_best_candidate[rank_index]] == nil then
				statistics[rank_name].names[otu.ranks_best_candidate[rank_index]] = 0
			end
			statistics[rank_name].names[otu.ranks_best_candidate[rank_index]] = statistics[rank_name].names[otu.ranks_best_candidate[rank_index]] + tonumber(otu.size)
		end
	end
	if otu.rank_resolved == nil then
		if statistics["-"].names["-"] == nil then
			statistics["-"].names["-"] = 0
		end
		statistics["-"].names["-"] = statistics["-"].names["-"] + tonumber(otu.size)
	end
end

-- Write output

file_output = io.open(file_output_path, "w")
file_output:write("#sample\totu\tsize\tsequence\taccession\tcoverage\tsimilarity\ttaxon\trank_resolved")
for rank_index, rank_name in ipairs(ranks) do
	file_output:write("\trank(" .. rank_name .. ")")
end
for rank_index, rank_name in ipairs(ranks) do
	file_output:write("\tname(" .. rank_name .. ")")
end
file_output:write("\n")
for sample_otu, otu in pairs(otus) do
	otu_index = otu_index + 1
	file_output:write(otu.sample .. "\t" ..
		otu.identification .. "\t" ..
		otu.size .. "\t" ..
		otu.sequence .. "\t" ..
		otu.best_hit.accession .. "\t" ..
		otu.best_hit.coverage .. "\t" ..
		otu.best_hit.similarity .. "\t" ..
		otu.best_hit.taxon .. "\t")
	if otu.rank_resolved ~= nil then
		file_output:write(otu.rank_resolved)
	else
		file_output:write("-")
	end
	-- Write ranks
	for rank_index, rank_name in ipairs(ranks) do
		file_output:write("\t")
		if otu.ranks_resolved[rank_index] == true then
			file_output:write(otu.ranks_best_candidate[rank_index])
		else
			local total = otu.ranks_best_candidate_total[rank_index]
			while total > 0 do
				for otu_rank_name, otu_rank_total in pairs(otu.ranks[rank_index]) do
					if otu_rank_total == total then
						file_output:write(otu_rank_name .. "(" .. otu_rank_total .. ")|")
					end
				end
				total = total - 1
			end
		end
	end
	file_output:write("\n")
	if otu_index % 200 == 0 or otu_index == otus_total then
		io.write("\rWriting output at " .. string.format("%0.2f", (otu_index / otus_total) * 100.0) .. "%.")
		io.flush()
	end
end
io.write("\n")
io.flush()
io.close(file_output)

-- Write output statistics

file_output = io.open(file_output_path .. ".statistics", "w")
file_output:write("#rank\tname\tabundance\n")
for rank_index, rank_name in ipairs(ranks) do
	for taxon_name, abundance in pairs(statistics[rank_name].names) do
		file_output:write(rank_name .. "\t" .. taxon_name .. "\t" .. abundance .. "\n")
	end
end
file_output:write("-\t-\t" .. statistics["-"].names["-"] .. "\n")
io.close(file_output)
