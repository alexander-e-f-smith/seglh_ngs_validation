<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-seglhqiaseq_logo_dark.png">
    <img alt="nf-core/seglhqiaseq" src="docs/images/nf-core-seglhqiaseq_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/seglhqiaseq/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/seglhqiaseq/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/seglhqiaseq/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/seglhqiaseq/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/seglhqiaseq/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/seglhqiaseq)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23seglhqiaseq-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/seglhqiaseq)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

# Introduction

**seglh_ngs_validation** is a bioinformatics pipeline for validation of NGS assays:
- currently the pipleine is configured around Synnovis/SEGLH(KCH) NGS pipelines and the specific outputs therein.

## Input

1. The pipeline requires at least VCF files to compare as a minimum; one VCF from the pipleine to be validated and another in a compatible format to use a the 'truth' to compare against (usually from the same variant caller)
2. Optional input  files  include (and should be paired for variant vcf and cnvkit bed file if cnv assessment is selected. pairing of other files is optional):
   - QC metrics file in json format (currently only supported for KCH/SEGLH snappy/DX based pipelines; can be supplied as single test file or with paired truth for comparison. This QC operation can be switched off via yaml config file sup[plied at run time (see belpW. Also whether a single of paired comparison is required is also specified via the yaml config file.
   - Bed file output for cnv calls (currently using that outputted by cnvkit)
   - Coverage file (exoncoverage) as outputted by kch/DX pipelines
   - CNVkit output files (in bed format; see cnvkit manual for details of bed conversion)
3. A yaml file for custom run configuration parameters (see yaml file configuration below). This contains options to run with paired or unpaired inputs if optional (see above)
4. `Samplesheet.csv` detailing input files including PATHS (see below for details).

## Yaml configuration file 
example of yaml format:
  ```
  pipeline_name : "generic" 
  pipeline_base : "generic"
  isec_filter_bed_rp : "resources/dna_mm/dna_mm_roi_305genes.bed" 
  isec_exclude_filter_1 : "-e 'FORMAT/AF<0.05'" 
  isec_exclude_filter_2 : "-e 'FORMAT/AF<0.02'" 
  cnv_loci_of_interest_bed : "resources/dna_mm/dna_mm_solid_cnvkit_loci_of_interest.bed" 
  expected_data_source : "truth_data" 
  observed_data_source : "test_data" 
  snappy_qc_run : "yes"
  cnvkit_compare_run : "yes" 
  info_tags_to_merge : "DP:join,AF:join" 
  info_tags_to_print : "" 
  githead_rp : ".git/ORIG_HEAD" 
  qc_standalone : "no" 
  ```
Yaml options:

| Option | Detail | Example |
| --- | --- | --- | 
| isec_filter_bed_rp | target regions (bed) being assessed  | resources/test/test.bed |

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->


## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):
-->
First, prepare a samplesheet with your input data that looks as follows (for example purposes):

`samplesheet_example.csv`:

```csv
sample,sample_vcf,truth_vcf,variant_vcf,json,truth_json,exon_cov,truth_exon_cov,cnv_bed,truth_cnv_bed
CONTROL_REP1,/path/to/Sample1_test.filtered.vcf.gz ,/path/to/Sample1_truth_filtered.vcf.gz,,,,,,,
```

Each row represents a sample and associted pipeline output files, where there is an option for a test and paired truth associated file for each file type (paired tructh and test required for variant vcf file input) .


Now, you can run the pipeline using this example (:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run main.nf  
   -profile docker,qiaseq 
   --outdir example_test1  
   --input assets/samplesheet_example.csv 
   -params-file pipeline_params/generic.yaml  
   --pipeline_resources "/home/seglh_ngs_validation/"
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/seglhqiaseq/usage) and the [parameter documentation](https://nf-co.re/seglhqiaseq/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/seglhqiaseq/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/seglhqiaseq/output).

## Credits

nf-core/seglhqiaseq was originally written by Alexander Smith.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#seglhqiaseq` channel](https://nfcore.slack.com/channels/seglhqiaseq) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/seglhqiaseq for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
