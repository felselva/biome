# Biome

This is a in-progress (sketch for a future repository) project to compile tools developed by me for bioinformatics.

- `randomfasta` (https://github.com/felselva/randomfasta) Generates FASTA files with random content based on a provided total of sequences, length of sequences and distribution of the units in the sequences.
- `taxonomy_blast_to_table.lua` converts a `.blast` output from `BLASTN` to a tab separated table.
- `add_taxon_id.lua` add the taxon ID using NCBI taxonomy database.
- `add_taxon_id_merged.lua` converts deprecated taxon IDs to newer IDs using NCBI taxonomy database.
- `add_taxon_lineage.lua` add the taxon ID at each taxonomy level.
- `report.lua` generate a table with the best results considering each taxon level.
- `add_report_taxon_names.lua` add taxon names to each taxon ID.
