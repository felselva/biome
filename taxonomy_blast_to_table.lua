local file_input_blast_path = nil
local file_input_blast = nil
local file_output_path = nil
local file_output = nil
local file_position = nil
local file_size = nil
local line = nil
local next_line = nil

file_input_blast_path = arg[1]
file_output_path = arg[2]
file_input_blast = io.open(file_input_blast_path, "r")
file_output = io.open(file_output_path, "w")
if file_input_blast == nil or file_output == nil then
	abort[1] = true
end

-- Convert taxonomy BLAST to table

file_output:write("#sample\totu\tsize\taccession\tcoverage\tsimilarity\n")

file_size = file_input_blast:seek("end")
file_input_blast:seek("set", 0)
next_line = file_input_blast:lines()
line = next_line()
while line ~= nil do
	local query = nil
	local target = nil
	local sample = nil
	local otu_number = nil
	local otu_size = nil
	local accession = nil
	local coverage = nil
	local similarity = nil
	query, target, coverage, similarity = string.match(line, "([%w%p]+)\t([%w%p ]+)\t([%d%p]+)\t([%d%p]+)")
	if query == nil or target == nil or coverage == nil or similarity == nil then
		print("Failed at line \"" .. line .. "\".")
		abort[1] = true
	end
	sample = string.match(query, "sample=([%w%.]+);")
	otu_number = string.match(query, "otu=(%d+)")
	otu_size = string.match(query, "size=(%d+)")
	accession = string.match(target, "([%w%p]+)%s")
	if sample == nil or otu_number == nil or otu_size == nil or accession == nil then
		abort[1] = true
	end
	file_output:write(sample .. "\t" ..
		otu_number .. "\t" ..
		otu_size .. "\t" ..
		accession .. "\t" ..
		coverage .. "\t" ..
		similarity .. "\n")
	line = next_line()
	file_position = file_input_blast:seek()
	if file_position % 10000 == 0 or file_position == file_size then
		io.write("\rReading input and writing output at " .. string.format("%0.2f", (file_position / file_size) * 100.0) .. "%.")
		io.flush()
	end
end
io.write("\n")
io.flush()
io.close(file_input_blast)
io.close(file_output)
