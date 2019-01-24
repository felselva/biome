# Pipeline written by Felipe Ferreira da Silva
# Github: @felselva
# Contact: felselva@gmail.com


# SETUP ########################################################################
################################################################################

PROGRAMS_DIR=${1}
SAMPLE_NAMES=${2}
SAMPLES_DIR=${3}
OUTPUT_DIR=${4}
REFERENCE_DB=${5}
TAXONOMY_DIR=${6}
MIN_LENGTH=${7}			# Minimum number of bases per sequence
MAX_EE=${8}			# 0.0-1.0
MAX_NS=${9}			# Maximum number allowed of N bases
CLUSTERING_ID=${10}		# 0.0-1.0
BLAST_QCOV_HSP_PERC=${11}	# 0.0-100.0
BLAST_PERC_IDENTITY=${12}	# 0.0-100.0
THREADS=4

VSEARCH=${PROGRAMS_DIR}/vsearch
MAKEBLASTDB=${PROGRAMS_DIR}/makeblastdb
BLASTN=${PROGRAMS_DIR}/blastn
MAP=${PROGRAMS_DIR}/mapperncbi.lua
SEPARATEBYSAMPLE=${PROGRAMS_DIR}/separate_by_sample.lua
DATABASE_ADD_TAXON=${PROGRAMS_DIR}/database_add_taxon.lua
: '
rm -f "${OUTPUT_DIR}/*"



# LOG PARAMETERS AND PROGRAM VERSIONS ##########################################
################################################################################

echo "# BLASTN VERSION" > "${OUTPUT_DIR}/log"
${BLASTN} -version >> "${OUTPUT_DIR}/log"
echo "" >> "${OUTPUT_DIR}/log"
echo "" >> "${OUTPUT_DIR}/log"
echo "" >> "${OUTPUT_DIR}/log"
echo "# VSEARCH VERSION" >> "${OUTPUT_DIR}/log"
${VSEARCH} --version &>> "${OUTPUT_DIR}/log"
echo "" >> "${OUTPUT_DIR}/log"
echo "" >> "${OUTPUT_DIR}/log"
echo "" >> "${OUTPUT_DIR}/log"
echo "# PARAMETERS" >> "${OUTPUT_DIR}/log"
echo "Taxonomy directory = ${TAXONOMY_DIR}" >> "${OUTPUT_DIR}/log"
echo "Reference database = ${REFERENCE_DB}" >> "${OUTPUT_DIR}/log"
echo "Reference database (search term):" >> "${OUTPUT_DIR}/log"
cat "${REFERENCE_DB}.description" >> "${OUTPUT_DIR}/log"
echo "Filtering minimum length = ${MIN_LENGTH}" >> "${OUTPUT_DIR}/log"
echo "Filtering maximum expected error = ${MAX_EE}" >> "${OUTPUT_DIR}/log"
echo "Filtering maximum N bases = ${MAX_NS}" >> "${OUTPUT_DIR}/log"
echo "Clustering threshold (ID) = ${CLUSTERING_ID}" >> "${OUTPUT_DIR}/log"
echo "Identification minimum coverage = ${BLAST_QCOV_HSP_PERC}" >> "${OUTPUT_DIR}/log"
echo "Identification minimum similarity (ID) = ${BLAST_PERC_IDENTITY}" >> "${OUTPUT_DIR}/log"



# COUNT SEQUENCES ##############################################################
################################################################################

echo ""
echo ""
echo ""
echo "COUNTING SEQUENCES"

COUNTER=0

while read SAMPLE_NAME
do
	COUNTER=$(($(eval "wc -l < ${SAMPLES_DIR}/${SAMPLE_NAME}.fastq | bc") + ${COUNTER}))
done < ${SAMPLE_NAMES}

COUNTER=$((${COUNTER} / 4))
echo "" >> "${OUTPUT_DIR}/log"
echo "" >> "${OUTPUT_DIR}/log"
echo "" >> "${OUTPUT_DIR}/log"
echo "# COUNT SEQUENCES" >> "${OUTPUT_DIR}/log"
echo "Initial total of sequences = ${COUNTER}" >> "${OUTPUT_DIR}/log"



# FILTER #######################################################################
################################################################################

echo ""
echo ""
echo ""
echo "FILTERING"

COUNTER=0

while read SAMPLE_NAME
do
	"${VSEARCH}" --threads ${THREADS} \
		--fastq_filter "${SAMPLES_DIR}/${SAMPLE_NAME}.fastq" \
		--fastq_maxee ${MAX_EE} \
		--fastq_minlen ${MIN_LENGTH} \
		--fastq_maxns ${MAX_NS} \
		--fasta_width 0 \
		--fastaout "${OUTPUT_DIR}/${SAMPLE_NAME}_filtered.fasta"
	COUNTER=$(($(eval "wc -l < ${OUTPUT_DIR}/${SAMPLE_NAME}_filtered.fasta | bc") + ${COUNTER}))
done < "${SAMPLE_NAMES}"

COUNTER=$((${COUNTER} / 2))
echo "" >> "${OUTPUT_DIR}/log"
echo "" >> "${OUTPUT_DIR}/log"
echo "" >> "${OUTPUT_DIR}/log"
echo "# FILTER" >> "${OUTPUT_DIR}/log"
echo "Total of sequences after filtering = ${COUNTER}" >> "${OUTPUT_DIR}/log"



# DEREPLICATE ##################################################################
################################################################################

echo ""
echo ""
echo ""
echo "DEREPLICATING"

while read SAMPLE_NAME
do
	"${VSEARCH}" --threads ${THREADS} \
		--derep_fulllength "${OUTPUT_DIR}/${SAMPLE_NAME}_filtered.fasta" \
		--strand both \
		--sizein \
		--sizeout \
		--fasta_width 0 \
		--relabel "sample=${SAMPLE_NAME};sequence_dereplicated_id=" \
		--uc "${OUTPUT_DIR}/${SAMPLE_NAME}_dereplicated.uc" \
		--output "${OUTPUT_DIR}/${SAMPLE_NAME}_dereplicated.fasta"
	cat "${OUTPUT_DIR}/${SAMPLE_NAME}_dereplicated.fasta" >> "${OUTPUT_DIR}/all_dereplicated.fasta"
done < "${SAMPLE_NAMES}"

COUNTER=$(eval "wc -l < ${OUTPUT_DIR}/all_dereplicated.fasta | bc")
COUNTER=$((${COUNTER} / 2))
echo "" >> "${OUTPUT_DIR}/log"
echo "" >> "${OUTPUT_DIR}/log"
echo "" >> "${OUTPUT_DIR}/log"
echo "# DEREPLICATE" >> "${OUTPUT_DIR}/log"
echo "Total of sequences after dereplication = ${COUNTER}" >> "${OUTPUT_DIR}/log"



# REMOVE CHIMERAS ##############################################################
################################################################################

echo ""
echo ""
echo ""
echo "REMOVING CHIMERAS"

"${VSEARCH}" --uchime_denovo "${OUTPUT_DIR}/all_dereplicated.fasta" \
	--sizein \
	--sizeout \
	--fasta_width 0 \
	--nonchimeras "${OUTPUT_DIR}/all_non_chimera.fasta"

COUNTER=$(eval "wc -l < ${OUTPUT_DIR}/all_non_chimera.fasta | bc")
COUNTER=$((${COUNTER} / 2))
echo "" >> "${OUTPUT_DIR}/log"
echo "" >> "${OUTPUT_DIR}/log"
echo "" >> "${OUTPUT_DIR}/log"
echo "# REMOVE CHIMERAS" >> "${OUTPUT_DIR}/log"
echo "Total of sequences after removing chimeras = ${COUNTER}" >> "${OUTPUT_DIR}/log"



# SEPARATE NON CHIMERAS BY SAMPLE ##############################################
################################################################################

echo ""
echo ""
echo ""
echo "# SEPARATING NON CHIMERAS BY SAMPLE"

lua ${SEPARATEBYSAMPLE} "${OUTPUT_DIR}/all_non_chimera.fasta" "${OUTPUT_DIR}" "_non_chimera.fasta"



# CLUSTER ######################################################################
################################################################################

echo ""
echo ""
echo ""
echo "# CLUSTERING"

while read SAMPLE_NAME
do
	echo "Clustering ${OUTPUT_DIR}/${SAMPLE_NAME}_non_chimera.fasta"
	"${VSEARCH}" --threads ${THREADS} \
		--cluster_fast "${OUTPUT_DIR}/${SAMPLE_NAME}_non_chimera.fasta" \
		--id ${CLUSTERING_ID} \
		--sizein \
		--sizeout \
		--fasta_width 0 \
		--relabel "sample=${SAMPLE_NAME};otu=" \
		--uc "${OUTPUT_DIR}/${SAMPLE_NAME}_clusters.uc" \
		--centroids "${OUTPUT_DIR}/${SAMPLE_NAME}_otus.fasta"
done < "${SAMPLE_NAMES}"



# CONCATENATE OTUS IN SINGLE FILE ##############################################
################################################################################

echo ""
echo ""
echo ""
echo "CONCATENATING OTUS AND CLUSTERS IN SINGLE FILE"

while read SAMPLE_NAME
do
	echo "Concatenating ${OUTPUT_DIR}/${SAMPLE_NAME}_otus.fasta"
	cat "${OUTPUT_DIR}/${SAMPLE_NAME}_otus.fasta" >> "${OUTPUT_DIR}/all_otus.fasta"
done < "${SAMPLE_NAMES}"

COUNTER=$(eval "wc -l < ${OUTPUT_DIR}/all_otus.fasta | bc")
COUNTER=$((${COUNTER} / 2))
echo "" >> "${OUTPUT_DIR}/log"
echo "" >> "${OUTPUT_DIR}/log"
echo "" >> "${OUTPUT_DIR}/log"
echo "# CLUSTER" >> "${OUTPUT_DIR}/log"
echo "Total of clusters = ${COUNTER}" >> "${OUTPUT_DIR}/log"



# MAKE DATABASE ################################################################
################################################################################

echo ""
echo ""
echo ""
echo "MAKING DATABASE"

"${MAKEBLASTDB}" -dbtype nucl \
	-in "${REFERENCE_DB}" \



# TAXONOMY #####################################################################
################################################################################

echo ""
echo ""
echo ""
echo "IDENTIFYING TAXONOMY"

"${BLASTN}" -num_threads ${THREADS} \
	-query "${OUTPUT_DIR}/all_otus.fasta" \
	-db "${REFERENCE_DB}" \
	-qcov_hsp_perc ${BLAST_QCOV_HSP_PERC} \
	-perc_identity ${BLAST_PERC_IDENTITY} \
	-outfmt '6 qseqid stitle qcovhsp pident' \
	-out "${OUTPUT_DIR}/taxonomy.blast"
'


# MAP ##########################################################################
################################################################################

echo ""
echo ""
echo ""
echo "# MAKING TABLE"

lua ${MAP} "${SAMPLES_DIR}/info" "${OUTPUT_DIR}" "${TAXONOMY_DIR}"



# ALL DONE #####################################################################
