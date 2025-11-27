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
include { EXTRACT_KCH_QC         } from '../modules/local/kch_qc_extract/main.nf'
include { COMBINE_KCH_QC         } from '../modules/local/kch_qc_combine/main.nf'
include { ISEC_VALIDATION        } from '../modules/local/bcftools/isec/main.nf'
include { BATCH_CONCORDANT_VARIANTS_PLOTTING_VARDICT    } from '../modules/local/batch_concordant_variants_plot_vardict/main.nf'
include { CONCORDANT_VARIANTS_PLOTTING_VARDICT} from '../modules/local/concordant_variants_plot_vardict/main.nf'
include { BATCH_SENSITIVITY_SPECIFICITY    } from '../modules/local/batch_unique_variants_merge/main.nf'
include { COMPARE_CNVKIT_BED        } from '../modules/local/bedtools/cnvkit_bed_compare/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SEGLHQIASEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input and processed in pipeline initiation task
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()


    // collect resources
    isec_target_file = Channel.fromPath(params.isec_filter_bed).collect()
    githead_file = Channel.fromPath(params.githead).collect()

    // organize input files into channels
    vcfs_ch = ch_samplesheet.map { meta, sample_vcf, truth_vcf, variant_vcf, json, truth_json, exon_cov, truth_exon_cov, cnv_bed, truth_cnv_bed-> tuple(meta, sample_vcf, truth_vcf, variant_vcf) }
    vcfs_combine_ch = ch_samplesheet.map { meta, sample_vcf, truth_vcf, variant_vcf, json, truth_json, exon_cov, truth_exon_cov, cnv_bed, truth_cnv_bed-> tuple(sample_vcf) }.collect()
    truth_vcfs_combine_ch = ch_samplesheet.map { meta, sample_vcf, truth_vcf, variant_vcf, json, truth_json, exon_cov, truth_exon_cov, cnv_bed, truth_cnv_bed-> tuple(truth_vcf) }.flatten().collect()
    json_ch = (ch_samplesheet).map {meta, sample_vcf, truth_vcf, variant_vcf, json, truth_json, exon_cov, truth_exon_cov, cnv_bed, truth_cnv_bed-> tuple(meta, json, truth_json, exon_cov, truth_exon_cov) }
    //json_ch.view()
    cnv_bed_ch = (ch_samplesheet).map {meta, sample_vcf, truth_vcf, variant_vcf, json, truth_json, exon_cov, truth_exon_cov, cnv_bed, truth_cnv_bed-> tuple(meta, cnv_bed, truth_cnv_bed) }
    
    
    //json_combine_ch = (ch_samplesheet).map {meta, sample_vcf, truth_vcf, variant_vcf, json, truth_json, exon_cov, truth_exon_cov-> tuple(json, truth_json, exon_cov, truth_exon_cov) }.flatten().collect()

    if (params.snappy_qc_run == "yes") {
    //
    // Module: Run EXTRACT_QC
    //
    //resources_file = file(params.qiagen_adapters)
        extract_qc_ch = EXTRACT_KCH_QC(json_ch, githead_file)
        ch_versions = ch_versions.mix(EXTRACT_KCH_QC.out.versions.first())
    //   EXTRACT_KCH_QC.out.norm_index.view()
    //
    // Module: Run COMBINE_QC
    //
        combine_kch_qc_files_ch = EXTRACT_KCH_QC.out.sample_kch_qc.map { meta, qc_files -> tuple( qc_files ) }.collect(sort: {it.baseName})
        combine_qc_ch = COMBINE_KCH_QC(combine_kch_qc_files_ch) 
        ch_versions = ch_versions.mix(COMBINE_KCH_QC.out.versions.first()) }
    else {}

    //
    // MODULE: Run bcftools isec
    //
    ISEC_VALIDATION (
      vcfs_ch, isec_target_file
    )
    //ch_multiqc_files = ch_multiqc_files.mix(ISEC_VALIDATION.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(ISEC_VALIDATION.out.versions.first())
    //ISEC_VALIDATION.out.isec_concordant_sample_A.view()
    //ISEC_VALIDATION.out.isec_unique_variants.view()
    
    //
    // Module: Run CONCORDANT_VARIANTS_PLOTTING
    //
    //resources_file = file(params.qiagen_adapters)
    collect_sample_concordant_input = Channel.empty()
    collect_sample_concordant_input = collect_sample_concordant_input.concat(
            ISEC_VALIDATION.out.isec_concordant_sample_A,
            ISEC_VALIDATION.out.isec_concordant_truth_A, 
            ISEC_VALIDATION.out.isec_concordant_sample_B, 
            ISEC_VALIDATION.out.isec_concordant_truth_B, 
            ISEC_VALIDATION.out.isec_vcf_indexes,
            EXTRACT_KCH_QC.out.norm_index
        )
        .groupTuple().map {rg, files -> [rg, *files] }
    collect_sample_concordant_input.view()

    vcfs_plot_ch = CONCORDANT_VARIANTS_PLOTTING_VARDICT(collect_sample_concordant_input)
    ch_versions = ch_versions.mix(CONCORDANT_VARIANTS_PLOTTING_VARDICT.out.versions.first())


    //
    // Module: Run BATCH_CONCORDANT_VARIANTS_PLOTTING_VARDICT 
    //
    //resources_file = file(params.qiagen_adapters)
    combine_sample_concordant_A = (ISEC_VALIDATION.out.isec_concordant_sample_A).map {meta, vcf-> tuple(vcf) }.collect(sort: {it.baseName})
    combine_truth_concordant_A = (ISEC_VALIDATION.out.isec_concordant_truth_A).map {meta, vcf-> tuple(vcf) }.collect(sort: {it.baseName})
    combine_sample_concordant_B = (ISEC_VALIDATION.out.isec_concordant_sample_B).map {meta, vcf-> tuple(vcf) }.collect(sort: {it.baseName})
    combine_truth_concordant_B = (ISEC_VALIDATION.out.isec_concordant_truth_B).map {meta, vcf-> tuple(vcf) }.collect(sort: {it.baseName})
    combine_vcf_indexes = (ISEC_VALIDATION.out.isec_vcf_indexes).map {meta, indexes-> tuple(indexes) }.collect(sort: {it.baseName})
    merge_vcfs_plot_ch = BATCH_CONCORDANT_VARIANTS_PLOTTING_VARDICT(combine_sample_concordant_A, combine_truth_concordant_A, combine_sample_concordant_B, combine_truth_concordant_B, combine_vcf_indexes)
    ch_versions = ch_versions.mix(BATCH_CONCORDANT_VARIANTS_PLOTTING_VARDICT.out.versions.first())

    //
    // Module: Run BATCH_SENSITIVITY_SPECIFICITY
    //
    //resources_file = file(params.qiagen_adapters)
    batch_sample_unique = (ISEC_VALIDATION.out.isec_unique_sample_variants).map {meta, vcf-> tuple(vcf) }.collect(sort: {it.baseName})
    batch_truth_unique = (ISEC_VALIDATION.out.isec_unique_truth_variants).map {meta, vcf-> tuple(vcf) }.collect(sort: {it.baseName})
    batch_variant_stats = (ISEC_VALIDATION.out.isec_variant_comparison_stats).map {meta, stats-> tuple(stats) }.collect(sort: {it.baseName})
    batch_detection_performance_ch = BATCH_SENSITIVITY_SPECIFICITY(batch_sample_unique, batch_truth_unique, batch_variant_stats, combine_vcf_indexes)
    ch_versions = ch_versions.mix(BATCH_SENSITIVITY_SPECIFICITY.out.versions.first())    
    //batch_sample_unique.view()
 
    if (params.cnvkit_compare_run == "yes") {
    //
    // Module: Run COMPARE_CNVKIT_BED (with collect target bed); conditional running
    //
        cnv_lociofinterest = Channel.fromPath(params.cnv_loci_of_interest_bed).collect()
        cnvkit_compare_ch = COMPARE_CNVKIT_BED(cnv_bed_ch, cnv_lociofinterest) 
        ch_versions = ch_versions.mix(COMPARE_CNVKIT_BED.out.versions.first()) }
    else {}

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
