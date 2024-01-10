/*
 * Prepares index for bowtie2
 * Input: [FASTA] Reference file 
 * Output: Tuple of [FASTA] reference file and [BT2] index files generated for the reference   
 */
process build_index_bowtie {

	input:
	path(ref)

	output:
	tuple path("${ref}"), path("${ref}.*"), emit: index
	path("${task.process}.version.txt"), 	emit: version

	"""
	bowtie2-build ${ref} ${ref}
	
	echo -e "${task.process}\tbowtie2\t\$(bowtie2-build --version | head -1 | rev | cut -f1 -d' ' | rev)" > ${task.process}.version.txt
	"""
}

/*
 * Alignes reads to a reference via bowtie2, filters out unaligned sequeces
 * and converts output to BAM. Reference needs to be indexed.
 * Input: Tuple of [FASTA] reference file and [BT2] index files generated for the reference
 * 		[FASTA] Read files to be aligned to the reference
 * Params:  params.report_all_alignments    -> Reports all possible alignments and dismisses params.max_alignments
 *          params.max_alignments           -> Amount of alignemnts to report. Only used when params.report_all_alignments is false
 * Output: bam_alignments 	-> [BAM] Aligned sequences
 * 		  report_alignments -> [TXT] Alignment reports
 */
process mapping_bowtie{
	tag {query.simpleName}
	publishDir "${params.output_dir}/alignments", mode: 'copy', pattern: "${query.simpleName}.bam"
	publishDir "${params.output_dir}/statistics", mode: 'copy', pattern: "${query.simpleName}.statistics.txt"

	input:
	tuple path(ref), path(index)
	path(query)

	output:
	path "${query.simpleName}.bam", 			emit: bam_alignments
	path "${query.simpleName}.statistics.txt", 	emit: report
	path("${task.process}.version.txt"), 		emit: version

	script:
	def all_alignments = params.report_all_alignments ? '-a' : ''
	def some_alignments = params.max_alignments && !params.report_all_alignments ? "-k " + params.max_alignments : ''

	"""
	bowtie2 --no-unal \
		-q \
		${all_alignments} \
		${some_alignments} \
		-p ${task.cpus} \
		--seed 0 \
		-U ${query} \
		-x ${ref} \
		2> ${query.simpleName}.statistics.txt | samtools view -bS - > ${query.simpleName}.bam

	echo -e "${task.process}\tbowtie2\t\$(bowtie2 --version | head -1 | rev | cut -f1 -d' ' | rev)" > ${task.process}.version.txt
	echo -e "${task.process}\tsamtools\t\$(samtools --version | head -1 | rev | cut -f1 -d' ' | rev)" >> ${task.process}.version.txt
	"""
}

/*
 * Prepares index for STAR. Includes alignment if one is given
 * Input: [FASTA] Reference sequence
 *		[GTF]|[GFF3] Annotation file - if none is given it should say [NO_FILE]
 * Output: [DIR] Directory with index generated for given reference file
 */
process build_index_STAR {

	input:
	path(referenceGenome)
	path(gtf)

	output:
	path(index), 							emit: index
	path("${task.process}.version.txt"), 	emit: version


	script:
	annotation_file = !params.annotation == 'NO_FILE' ? '--sjdbGTFfile ' + ${gtf} : ''

	"""
	mkdir index
	STAR --runThreadN ${task.cpus} \
		--runMode genomeGenerate \
		--genomeDir ./index \
		--genomeFastaFiles ${referenceGenome} \
		${annotation_file}
	
	echo -e "${task.process}\tSTAR\t\$(STAR --version)" > ${task.process}.version.txt
	"""
}

/*
 * Alignes reads to reference via STAR. Reference needs to be index before.
 * Output is provided as BAM file that is sorted by coordinates
 * Input: Tuple of [FASTQ] Read files to be aligned and [DIR] Directory containing the STAR index
  * Params:  params.report_all_alignments    -> Reports all possible alignments and dismisses params.max_alignments
 *          params.max_alignments           -> Amount of alignemnts to report. Only used when params.report_all_alignments is false
 * Output: bam_alignments -> [BAM] Aligned reads sorted by coordinates 
 *		report_alignments -> [TXT] Alignment report
 */
process mapping_STAR{
	tag {query.simpleName}
	publishDir "${params.output_dir}/alignments", mode: 'copy', pattern: "${query.baseName}.Aligned.sortedByCoord.out.bam"
	publishDir "${params.output_dir}/statistics", mode: 'copy', pattern: "${query.baseName}.Log.*"

	input:
	tuple path(query), path(indexDir)

	output:
	path("${query.baseName}.Aligned.sortedByCoord.out.bam"), 	emit: bam_alignments
	path("${query.baseName}.Log.*"), 							emit: report
	path("${task.process}.version.txt"), 						emit: version

	script:
	def all_alignments = params.report_all_alignments ? '--outSAMmultNmax -1' : ''
	def some_alignments = params.max_alignments && !params.report_all_alignments ? "--outSAMmultNmax " + params.max_alignments : ''

	"""
	STAR --runThreadN ${task.cpus} \
		--genomeDir ${indexDir} \
		--readFilesIn ${query} \
		--outFileNamePrefix ${query.baseName}. \
		${all_alignments} \
		${some_alignments} \
		--outSAMtype BAM SortedByCoordinate

	echo -e "${task.process}\tSTAR\t\$(STAR --version)" > ${task.process}.version.txt
	"""
}