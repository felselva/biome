local file_input_path = nil
local file_input = nil
local taxonomy_directory = nil
local file_output_path = nil
local file_output = nil
local file_position = nil
local file_size = nil
local check_accessions_taxon_id = nil
local read_accessions_taxon_id_from_file = nil
local line = nil
local next_line = nil
local hits = {}
local accessions = {}
local accessions_total = 0
local accessions_with_taxon_id_total = 0

file_input_path = arg[1]
taxonomy_directory = arg[2]
file_output_path = arg[3]
file_input = io.open(file_input_path, "r")
file_output = io.open(file_output_path, "w")
if file_input == nil or file_output == nil then
	abort[1] = true
end

-- Store hits and accessions

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
	hit.sample, hit.otu_identification, hit.size, hit.accession, hit.coverage, hit.similarity =
		string.match(line, "([%w%.]+)\t(%d+)\t(%d+)\t([%w%p]+)\t([%d%p]+)\t([%d%p]+)")
	if string.sub(hit.accession, 1, 4) == "pdb|" then
		hit.accession = string.match(hit.accession, "pdb|(%w+)|")
	end
	if accessions[hit.accession] == nil then
		accessions[hit.accession] = true
		accessions_total = accessions_total + 1
	end
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

-- Store the taxon ID of the accessions

check_accessions_taxon_id = function()
	local complete = nil
	complete = true
	for accession, accession_state in pairs(accessions) do
		if accession_state == true then
			complete = false
		end
	end
	return complete
end

read_accessions_taxon_id_from_file = function(file_input_accession_to_taxon_id_name)
	print("Trying to open the file \"" .. taxonomy_directory .. "/" .. file_input_accession_to_taxon_id_name .. "\".")
	file_input = io.open(taxonomy_directory .. "/" .. file_input_accession_to_taxon_id_name, "r")
	if file_input ~= nil then
		file_size = file_input:seek("end")
		file_input:seek("set", 0)
		next_line = file_input:lines()
		line = next_line()
		line = next_line()
		while line ~= nil and accessions_with_taxon_id_total < accessions_total do
			local accession = nil
			local accession_version = nil
			local tmp = nil
			local taxon_id = nil
			accession, accession_version, taxon_id = string.match(line, "([%w%p_]+)\t([%w%p_]+)\t([%d]+)")
			if accessions[accession] == true then
				accessions[accession] = taxon_id
				accessions_with_taxon_id_total = accessions_with_taxon_id_total + 1
			end
			if accessions[accession_version] == true then
				accessions[accession_version] = taxon_id
				accessions_with_taxon_id_total = accessions_with_taxon_id_total + 1
			end
			-- Try the PDB format
			tmp = string.match(accession, "(%w+)_%w")
			if tmp ~= nil then
				if accessions[tmp] == true then
					accessions[tmp] = taxon_id
					accessions_with_taxon_id_total = accessions_with_taxon_id_total + 1
				end
				accession_tmp = string
			end
			line = next_line()
			file_position = file_input:seek()
			if file_position % 100000 == 0 or file_position == file_size then
				io.write("\rReading input \"" .. file_input_accession_to_taxon_id_name .. "\" at " .. string.format("%0.2f", (file_position / file_size) * 100) .. "%, with " .. tostring(accessions_with_taxon_id_total) .. " of " .. tostring(accessions_total) .. " accessions with taxon ID found.")
				io.flush()
			end
		end
		io.write("\n")
		io.flush()
		io.close(file_input)
	end
end

if check_accessions_taxon_id() == false then
	read_accessions_taxon_id_from_file("pdb.accession2taxid")
end
if check_accessions_taxon_id() == false then
	read_accessions_taxon_id_from_file("nucl_gb.accession2taxid")
end
if check_accessions_taxon_id() == false then
	read_accessions_taxon_id_from_file("nucl_est.accession2taxid")
end
if check_accessions_taxon_id() == false then
	read_accessions_taxon_id_from_file("nucl_gss.accession2taxid")
end
if check_accessions_taxon_id() == false then
	read_accessions_taxon_id_from_file("nucl_wgs.accession2taxid")
end
if check_accessions_taxon_id() == false then
	read_accessions_taxon_id_from_file("prot.accession2taxid")
end

-- Check if all accessions have their taxon ID

if check_accessions_taxon_id() == false then
	print("Failed to obtain the taxon ID of all accessions.")
	print("Accessions without taxon ID:")
	for accession, accession_taxon_id in pairs(accessions) do
		if accession_taxon_id == true then
			print("  " .. accession)
		end
	end
	abort[1] = true
end

-- Save to file output

file_output:write("#sample\totu\tsize\taccession\tcoverage\tsimilarity\ttaxon\n")
for hit_index, hit in ipairs(hits) do
	file_output:write(hit.sample .. "\t" ..
		hit.otu_identification .. "\t" ..
		hit.size .. "\t" ..
		hit.accession .. "\t" ..
		hit.coverage .. "\t" ..
		hit.similarity .. "\t" ..
		accessions[hit.accession] .. "\n")
	if hit_index % 10000 == 0 or hit_index == #hits then
		io.write("\rWriting output at " .. string.format("%0.2f", (hit_index / #hits) * 100.0) .. "%.")
		io.flush()
	end
end
io.write("\n")
io.flush()
io.close(file_output)
