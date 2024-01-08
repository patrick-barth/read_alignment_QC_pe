/*
 * Checks reads for general metrics
 * Input: [FASTQ] Unpreprocessed reads 
 * Output: [HTML] General report for unpreprocessed reads 
 */

process quality_control {
	tag {id}
	publishDir "${params.output_dir}/statistics/qc-preprocessing", mode: 'copy', pattern: "*_fastqc.{html,zip}"

	
	input:
	tuple val(id), path(reads)

	output:
	path "*_fastqc.{html,zip}",				emit: output
	path("${task.process}.version.txt"), 	emit: version

	"""
	fastqc ${reads[0]} ${reads[1]} -o .

	echo -e "${task.process}\tFastQC\t\$(fastqc --version | rev | cut -f 1 -d' ' | rev)" > ${task.process}.version.txt
	"""
}


/*
 * Checks preprocessed reads for general metrics
 * Input: [FASTQ] Preprocessed reads 
 * Output: [HTML] General report for preprocessed reads
 */
process quality_control_2 {
	tag {id}
	
	input:
	tuple val(id), path(reads1), path(reads2)

	output:
	path "*_fastqc.{html,zip}",				emit: output
	path("${task.process}.version.txt"), 	emit: version

	"""
	cat ${reads1} > ${id}_R1.processed.fastq
	cat ${reads2} > ${id}_R2.processed.fastq
	fastqc ${id}_R1.processed.fastq ${id}_R2.processed.fastq -o .

	echo -e "${task.process}\tFastQC\t\$(fastqc --version | rev | cut -f 1 -d' ' | rev)" > ${task.process}.version.txt
	"""
}

/*
 * Removes adapters from reads
 * Input: [FASTQ] Read file
 * Params: params.min_length -> Minimum length reads need after trimming to not be omitted
 * Output:  fastq_trimmed 	-> [FASTQ] Read file with afdapters trimmed
 *			report_trimming -> [TXT] Report adapter trimming
 */
process adapter_removal {
	tag {id}
	publishDir "${params.output_dir}/statistics/adapter_removal", mode: 'copy', pattern: "*_trimming_report.txt"

	input:
	tuple val(id), path(reads)

	output:
	tuple val(id), path("${reads[0].simpleName}_val_1.fq"), path("${reads[1].simpleName}_val_2.fq"), 	emit: reads
	path "*_trimming_report.txt", 																		emit: report
	tuple path("${task.process}.version.txt"), path("${task.process}.version2.txt"), 					emit: version

	"""
	trim_galore --cores ${task.cpus} -o . --length ${params.min_length} --quality 0 --paired ${reads[0]} ${reads[1]}

	echo -e "${task.process}\ttrim_galore\t\$(trim_galore -v | head -4 | tail -1 | sed -e 's/^[ \t]*//' | rev | cut -f 1 -d' ' | rev)" > ${task.process}.version.txt
	echo -e "${task.process}\tcutadapt\t\$(cutadapt --version)" > ${task.process}.version2.txt
	"""
}

/*
 * Removes bases with low quality from reads
 * Input: [FASTQ] Read file
 * Params: 	params.min_qual					-> Bases below this threshold are omitted 
 *			params.min_length 				-> Minimum length reads need after trimming to not be omitted 
 * Output: 	fastq_quality_filtered 	-> [FASTQ] Read file with low quality bases filtered out
 *			report_quality_filter 	-> [TXT] Report of quality filtering
 */
process quality_filter {
	tag {id}
	publishDir "${params.output_dir}/statistics/quality_filter", mode: 'copy', pattern: "${id}_report_quality.txt"
	publishDir "${params.output_dir}/processed-reads", mode: 'copy', pattern: "*.qtrim.fastq"

	input:
	tuple val(id), path(reads1), path(reads2)

	output:
	tuple val(id), path("${id}_R1.qtrim.fastq"), path("${id}_R2.qtrim.fastq"), 	emit: reads 
	path("${id}_report_quality.txt"), 											emit: report 
	path("${task.process}.version.txt"), 										emit: version


	"""
	sickle pe -f ${reads1} -r ${reads2} -t sanger -o ${id}_R1.qtrim.fastq -p ${id}_R2.qtrim.fastq -s ${id}_single.qtrim.fastq \
		-q ${params.min_qual} -l ${params.min_length} > ${id}_report_quality.txt

		
	echo -e "${task.process}\tsickle\t\$(sickle --version | head -1 | rev | cut -f1 -d' ' | rev)" > ${task.process}.version.txt
	"""
}