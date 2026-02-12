<h1>
  SEGLH NGS Validation Workflow
</h1>

# Introduction

**seglh_ngs_validation** is a bioinformatics pipeline for validation of NGS assays:
- currently the pipleine is configured around Synnovis/SEGLH(KCH) NGS pipelines and the specific outputs therein.
- currently this workflow will only successfully run using Nextflow 24.04.4

## Input

1. The pipeline requires at least VCF files to compare as a minimum; one VCF from the pipleine to be validated and another in a compatible format to use as the 'truth' to compare against (usually from the same variant caller)
2. Optional input files include (and should be paired for variant vcf and cnvkit bed file if cnv assessment is selected. pairing of other files is optional):
   - QC metrics file in json format (currently only supported for KCH/SEGLH snappy/DX based pipelines; can be supplied as single test file or with paired truth for comparison. This QC operation can be switched off via yaml config file supplied at run time see below.) Also whether a single or paired comparison is required is also specified via the yaml config file.
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
| isec_exclude_filter_1 | AF thresholds to exclude on sample vcf  | "-e 'FORMAT/AF<0.05'" |
| isec_exclude_filter_2 | AF thresholds to exclude on truth vcf  | "-e 'FORMAT/AF<0.05'" |
| cnv_loci_of_interest_bed | target regions (bed) being assessed for CNVs  | resources/test/cnvkit_loci_of_interest.bed |
| cnvkit_compare_run | Boolean for applying cnv analysis  | "yes" |
| expected_data_source | Truth Data Name (For Graph Visuals)  | "Truth Data" |
| observed_data_source | Validation Data Name (For Graph Visuals) | "Validation Data" |


<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

## Generating Samplesheet

The Samplesheet needed for this workflow can be generated using the `assets/samplesheet_generatory.py` script as follows:

```python
python3 assets/samplesheet_generator.py \
    -truth_set {filepath to truth dataset} \
    -exp_set {filepath to experimental dataset} \
    -output {complete filepath & name to output samplesheet}
```
The truth and experimental datasets structure should be a parent folder of the dataset and then a folder inside for each filetype to provide. For example:
```
Truth_Data_Set Folder
├── truth_json
│   └── *.json
├── truth_cnvs
│   └── *cnv_call_baf_postFilter.bed
├── truth_exoncov
│   └── *exoncoverage_metrics
└── truth_vcfs
    └── *.vcf.gz
```

This generator currently does not automatically import the variant vcf option, this will need to be added manually.

## Usage

Currently, this workflow will only work using **Nextflow 24.04.4**. Any later versions will not work.


<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):
-->
First, prepare a samplesheet with your input data by using the instructions above. If samplesheet needs to be made manually, make sure that it looks as follows (for example purposes):

`samplesheet_example.csv`:

```csv
sample,sample_vcf,truth_vcf,variant_vcf,json,truth_json,exon_cov,truth_exon_cov,cnv_bed,truth_cnv_bed
CONTROL_REP1,/path/to/Sample1_test.filtered.vcf.gz ,/path/to/Sample1_truth_filtered.vcf.gz,,,,,,,
```

Each row represents a sample and associated pipeline output files, where there is an option for a test and paired truth associated file for each file type (paired truth and test required for variant vcf file input).


Now, you can run the pipeline using this example:

```bash
nextflow run main.nf  
   -profile docker,qiaseq 
   --outdir example_test1  
   --input assets/samplesheet_example.csv 
   -params-file pipeline_params/generic.yaml  
   --pipeline_resources "/home/seglh_ngs_validation/"
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Multiple pipeline parameters files exist at `pipeline_params/`.
Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/seglhqiaseq/usage) and the [parameter documentation](https://nf-co.re/seglhqiaseq/parameters).

## Pipeline output

This pipeline produces results that analyse concordance and specificity between the two data sets.
### Concordance results 
  - contains all variants seen in the same sample in both data sets alongside other key metrics (AF, Depth, norm Depth)
  - contains visual graphs to compare variant loci depth and VAF frequency
  - merged vcfs containing all concordant variants
### Specificity results
  - list of unique variants from either data set and other key metrics (AF, Depth, MSI, Up/Down steam sequence)
  - overall sample statistics on unique variants seen in each sample
  - merged vcfs containing samples with unique variants only
### Combine results
  - Overall QC per sample (obs Q30 & Q40, total reads, duplication rate, insert size etc.)
  - visual graph comparing reads covered at 400X between the two datasets
### Sample specific results
  - Each sample has a breakdown of the results described above for that specific sample


## Credits

SEGLH/seglh_ngs_validation was originally written by Alexander Smith based on nf-core/seglhqiaseq.

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
