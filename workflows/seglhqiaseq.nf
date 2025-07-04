/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_seglhqiaseq_pipeline'
include { ISEC_VALIDATION          } from '../modules/local/bcftools/isec/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SEGLHQIASEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()


    // collect resources
    isec_target_file = Channel.fromPath(params.isec_filter_bed).collect()

    // organize input files into channels
    vcfs_ch = ch_samplesheet.map { meta, sample_vcf, truth_vcf, variant_vcf, json, truth_json-> tuple(meta, sample_vcf, truth_vcf, variant_vcf) }.view()

    //
    // MODULE: Run bcftools isec
    //
        ISEC_VALIDATION (
        vcfs_ch, isec_target_file
    )
    //ch_multiqc_files = ch_multiqc_files.mix(ISEC_VALIDATION.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(ISEC_VALIDATION.out.versions.first())
    
    combine_json_ch = (ch_samplesheet).map {meta, sample_vcf, truth_vcf, variant_vcf, json, truth_json -> tuple(json, truth_json)}.collect().view()
    //
    // another module
    //
    //resources_file_ch = file(params.qiagen_adapters)
    //another_ch = newmodule(ch_samplesheet, resources_file)
    //ch_versions = ch_versions.mix(trimmomatic.out.versions.first())


    //
    // umiextract module
    //
    //UMITOOLS_EXTRACT(trimmomatic.out.filtered_fastq)
    //ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions.first())
    //ch_multiqc_files = ch_multiqc_files.mix(UMITOOLS_EXTRACT.out.umiextract_metrics{it[1]})


    // cutadapt module
    //
    //CUTADAPT(UMITOOLS_EXTRACT.out.reads)
    //ch_versions = ch_versions.mix(CUTADAPT.out.versions.first())
    //ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT.out.log.collect{it[1]})

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
