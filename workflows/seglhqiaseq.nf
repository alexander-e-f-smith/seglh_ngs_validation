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
include { BATCH_CONCORDANT_VARIANTS_PLOTTING_VARDICT    } from '../modules/local/vardict_vcf_combine/main2.nf'
include { CONCORDANT_VARIANTS_PLOTTING_VARDICT} from '../modules/local/concordant_variants_plot_vardict/main.nf'

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
    vcfs_ch = ch_samplesheet.map { meta, sample_vcf, truth_vcf, variant_vcf, json, truth_json, exon_cov, truth_exon_cov-> tuple(meta, sample_vcf, truth_vcf, variant_vcf) }
    vcfs_combine_ch = ch_samplesheet.map { meta, sample_vcf, truth_vcf, variant_vcf, json, truth_json, exon_cov, truth_exon_cov-> tuple(sample_vcf) }.collect()
    truth_vcfs_combine_ch = ch_samplesheet.map { meta, sample_vcf, truth_vcf, variant_vcf, json, truth_json, exon_cov, truth_exon_cov-> tuple(truth_vcf) }.flatten().collect()
    json_ch = (ch_samplesheet).map {meta, sample_vcf, truth_vcf, variant_vcf, json, truth_json, exon_cov, truth_exon_cov-> tuple(meta, json, truth_json, exon_cov, truth_exon_cov) }
    json_combine_ch = (ch_samplesheet).map {meta, sample_vcf, truth_vcf, variant_vcf, json, truth_json, exon_cov, truth_exon_cov-> tuple(json, truth_json, exon_cov, truth_exon_cov) }.flatten().collect()
    vcfs_combine_ch.view()
    truth_vcfs_combine_ch
    json_combine_ch


    //
    // Module: Run COMBINE_QC
    //
    //resources_file = file(params.qiagen_adapters)
    //extract_qc = COMBINE_QC(json_ch)
    //ch_versions = ch_versions.mix(QC_COMBINE.out.versions.first())

    //
    // MODULE: Run bcftools isec
    //
        ISEC_VALIDATION (
        vcfs_ch, isec_target_file
    )
    //ch_multiqc_files = ch_multiqc_files.mix(ISEC_VALIDATION.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(ISEC_VALIDATION.out.versions.first())
    ISEC_VALIDATION.out.isec_concordant_sample.view()
    ISEC_VALIDATION.out.isec_unique_variants.view()
    
    //
    // Module: Run CONCORDANT_VARIANTS_PLOTTING
    //
    //resources_file = file(params.qiagen_adapters)
    vcfs_plot_ch = CONCORDANT_VARIANTS_PLOTTING_VARDICT(ISEC_VALIDATION.out.isec_concordant_sample, ISEC_VALIDATION.out.isec_concordant_truth)
    //ch_versions = ch_versions.mix(QC_COMBINE.out.versions.first())


    //
    // Module: Run MERGE_VCFS_PLOTTING 
    //
    //resources_file = file(params.qiagen_adapters)
    combine_sample_concordant = (ISEC_VALIDATION.out.isec_concordant_sample).map {meta, vcf-> tuple(vcf) }.collect()
    combine_truth_concordant = (ISEC_VALIDATION.out.isec_concordant_truth).map {meta, vcf-> tuple(vcf) }.collect()
    merge_vcfs_plot_ch = BATCH_CONCORDANT_VARIANTS_PLOTTING_VARDICT(combine_sample_concordant, combine_truth_concordant)
    combine_sample_concordant.view()
    //ch_versions = ch_versions.mix(QC_COMBINE.out.versions.first())


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
