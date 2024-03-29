#!/usr/bin/env nextflow

import groovy.json.JsonOutput // used for parameter output

nextflow.enable.dsl=2

include{
    collect_metadata
    get_md5sum
    multiqc
    collect_versions
} from './modules/default_processes.nf'

include {
    quality_control
    quality_control_2
    adapter_removal
    quality_filter
} from './modules/read_processing.nf'

if (params.aligner == "bowtie2"){
    include{
        build_index_bowtie
        mapping_bowtie
    } from './modules/alignment.nf'
} else if (params.aligner == "star"){
    include{
        build_index_STAR
        mapping_STAR
    } from './modules/alignment.nf'
}

/*
 * Prints help and exits workflow afterwards when parameter --help is set to true
 */

if ( params.help ) {
    help = """main.nf: A description of your script and maybe some examples of how
                |                to run the script
                |Required arguments:
                |   --reads         Location of the input file file (FASTQ).
                |   --reference     Reference sequence to align reads to (FASTA)
                |
                |Optional arguments:
                |   --min_length    Minimum length for reads after adapter trimming.
                |                   [default: ${params.min_length}]
                |   --min_qual      Minimum base quality.
                |                   [default: ${params.min_qual}]
                |   --min_percent_qual_filter   Minimum percentage of bases within a read that need to
                |                               be above the quality threshold
                |                               [default: ${params.min_percent_qual_filter}]
                |  -w            The NextFlow work directory. Delete the directory once the process
                |                is finished [default: ${workDir}]""".stripMargin()
    // Print the help with the stripped margin and exit
    println(help)
    exit(0)
}

//preparation for workflow

/*
 * Welcome log to be displayed before workflow
 */
log.info """\
        ${params.manifest.name} v${params.manifest.version}
        ==========================
        reads           : ${params.reads}
        reference       : ${params.reference}
        output to       : ${params.output_dir}
        --
        aligner         : ${params.aligner}
        --
        run as          : ${workflow.commandLine}
        started at      : ${workflow.start}
        config files    : ${workflow.configFiles}
        """
        .stripIndent()

//essential input files
input_reads     = Channel.fromFilePairs( params.reads )			//FASTQ file(s) containing reads
reference       = Channel.fromPath( params.reference )
//non essential input files
if(params.annotation != 'NO_FILE'){
    annotation_file = Channel.fromPath( params.annotation )
}
annotation = file(params.annotation)

// Collect all input files
input_files = input_reads
                .map{id,reads -> reads}
                .concat(Channel.of(annotation))
                .concat(reference)
                .flatten().toList()

/*
 * Starting subworkflow descriptions
 */

workflow alignment {
    take:
        reference
        annotation
        reads

    main:
        if(params.aligner == "bowtie2"){
            build_index_bowtie(reference)
            mapping_bowtie(build_index_bowtie.out.index.first(),
                            reads)

            alignments_tmp          =   mapping_bowtie.out.bam_alignments
            version_index_tmp       =   build_index_bowtie.out.version
            version_align_tmp       =   mapping_bowtie.out.version
            report_tmp              =   mapping_bowtie.out.report
        } else if (params.aligner == "star"){
            build_index_STAR(reference,
                            annotation)
            mapping_STAR(reads
                        .combine(build_index_STAR.out.index))

            alignments_tmp          =   mapping_STAR.out.bam_alignments
            version_index_tmp       =   build_index_STAR.out.version
            version_align_tmp       =   mapping_STAR.out.version
            report_tmp              =   mapping_STAR.out.report
        } 

    emit:
        version_index   =   version_index_tmp
        version_align   =   version_align_tmp
        reports         =   report_tmp

        alignments      =   alignments_tmp
}

/*
 * Actual workflow connecting subworkflows
 */
workflow {
    quality_control(input_reads)
    adapter_removal(input_reads)
    quality_filter(adapter_removal.out.reads)
    quality_control_2(quality_filter.out.reads)

    alignment(reference,
            annotation,
            quality_filter.out.reads
    )

    // Collect metadata
    collect_metadata()
    get_md5sum(input_files)
    collect_versions(collect_metadata.out.version
                        .concat(get_md5sum.out.version)
                        .concat(quality_control.out.version.first())
                        .concat(adapter_removal.out.version.first())
                        .concat(quality_filter.out.version.first())
                        .concat(quality_control_2.out.version.first())
                        .concat(alignment.out.version_index.first())
                        .concat(alignment.out.version_align.first())
                        .flatten().toList()
    )
}



/*
 * Prints completion status to command line
 */
workflow.onComplete{
	println "Pipeline completed at: $workflow.complete"
	println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
	println "Something went wrong :("
	println "Pipeline execution stopped with following error message: ${workflow.errorMessage}"
}