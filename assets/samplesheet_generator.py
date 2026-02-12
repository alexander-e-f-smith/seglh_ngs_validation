# This script takes one truth data set and one experimental data set and creates the samplesheet
# required to run through the nextflow validation pipeline

# The data sets parsed to this script need to have a respective folder for each filetype (vcf, cnv, json, exoncov) if applicable.
# Specific folder names are not needed, only that the filetypes are separated in individual folders within the data set folder

import argparse
import os
import glob
import pandas as pd

parser = argparse.ArgumentParser(
    prog = "Validation Samplesheet Generator",
    description = "Generates the samplesheet needed to run the validation pipeline" \
    "based on the two data sets provided to the script"
)
parser.add_argument('-truth_set', help="Filepath to the truth data set")
parser.add_argument('-exp_set', help="Filepath to the experimental data set")
parser.add_argument('-output', help="Complete filepath and name of the samplesheet to be generated")
args = parser.parse_args()

truth_set = args.truth_set
exp_set = args.exp_set
output = args.output

## FUNCTIONS

def get_sample_file_list(set):
    """
    Returns a list of samples based on the data set given
        :param set (str):   Data set filepath
        :return (list):     List of sample names
    """
    set_folder_names = os.listdir(set)
    # Uses filenames in first folder to obtain sample names
    folder1 = set + "/" + set_folder_names[0]
    sample_file_list = os.listdir(folder1)

    sample_return = []
    # Depending on the filetype, uses a different separator to obtain
    #   sample names
    for sample in sample_file_list:
        if "_" in sample:
            sample = (sample.split("_"))[0]
            sample_return.append(sample)
        else:
            sample = (sample.split("."))[0]
            sample_return.append(sample)

    return sample_return

def pair_up_samples(truth_samples, experiment_samples):
    """
    Returns a list of truth and experimental samples paired up, deliminated
    by a :
        :param truth_samples (list)             List of Truth Samples
        :param experiment_samples (list)        List of Experiment Samples
        :return paired_samples (list)           List of Samples Pairs
    """

    # Sort samples to order the lists to then iterate through
    truth_samples_sorted = sorted(truth_samples)
    experiment_samples_sorted = sorted(experiment_samples)

    paired_samples = []

    # Iterate through each lists and extract the sample id part of the name
    for t_sample in truth_samples_sorted:
        for e_sample in experiment_samples_sorted:
            t_sample_id = t_sample.split("-")[0]
            e_sample_id = e_sample.split("-")[0]

            # Check that the sample IDs match
            if t_sample_id in e_sample_id or e_sample_id in t_sample_id or t_sample_id == e_sample_id:
                pair = t_sample + ":" + e_sample
                paired_samples.append(pair)

    return paired_samples

def get_filepaths(sample_set):
    """
    Returns a dictionary of the folder filepaths based on the filetype
    inside them
        :param sample_set (str)             Data set filepath
        :return filepath_dict (dict)        Dictionary of filepaths
    """
    filepath_dict = {}
    set_folder_names = os.listdir(sample_set)
    for folder in set_folder_names:
        file_list = os.listdir(f"{sample_set}{folder}")
        if ".json" in file_list[0]:
            json_folder = f"{sample_set}{folder}"
            filepath_dict["json"] = json_folder
        elif "cnv" in file_list[0]:
            cnv_bed_folder = f"{sample_set}{folder}"
            filepath_dict["cnv"] = cnv_bed_folder
        elif "vcf" in file_list[0]:
            vcf_folder = f"{sample_set}{folder}"
            filepath_dict["vcf"] = vcf_folder
        elif "exoncoverage" in file_list[0]:
            exoncov_folder = f"{sample_set}{folder}"
            filepath_dict["exoncov"] = exoncov_folder

    return filepath_dict

def find_sample_file(sample_name, filepath_dict, file_type, extension):
    """
    Returns the specific filename and filepath for samplesheet entry
        :param sample_name (str)        String of sample name
        :param filepath_dict (dict)     Dictionary of filepaths
        :param file_type (str)          String of filepath to search in dict
        :extension (str)                String of file suffix
        :return file path (str)         String of specific filepath of file
    """
    pattern = f"{filepath_dict[file_type]}/{sample_name}{extension}"
    return glob.glob(pattern)

##############################################################################################################################

# Obtain samples from folders
truth_samples = get_sample_file_list(truth_set)
experiment_samples = get_sample_file_list(exp_set)

# Pair up samples based on ID
paired_samples = pair_up_samples(truth_samples, experiment_samples)

# Obtain filepaths for each filetype and create dictionary
truth_filepath_dict = get_filepaths(truth_set)
experimental_filepath_dict = get_filepaths(exp_set)

# Obtain specific file paths for truth and experiment sets and add to dataframe
data = []
for pair in paired_samples:
    truth_sample, exp_sample = pair.split(":")
    sample_id = exp_sample.split("-")[0]

    # VCF files
    exp_sample_vcf = find_sample_file(exp_sample, experimental_filepath_dict, 'vcf', '*.vcf.gz')
    truth_sample_vcf = find_sample_file(truth_sample, truth_filepath_dict, 'vcf', '*.vcf.gz')
    
    # JSON files
    exp_sample_json = find_sample_file(exp_sample, experimental_filepath_dict, 'json', '*.json')
    truth_sample_json = find_sample_file(truth_sample, truth_filepath_dict, 'json', '*.json')
    
    # Exon coverage files
    exp_sample_exoncov = find_sample_file(exp_sample, experimental_filepath_dict, 'exoncov', '*_metrics')
    truth_sample_exoncov = find_sample_file(truth_sample, truth_filepath_dict, 'exoncov', '*_metrics')
    
    # CNV files
    # Conditional added as some runs may not have CNV files
    exp_sample_cnv = None
    truth_sample_cnv = None
    if "cnv" in truth_filepath_dict:
        exp_sample_cnv = find_sample_file(exp_sample, experimental_filepath_dict, 'cnv', '*.bed')
        truth_sample_cnv = find_sample_file(truth_sample, truth_filepath_dict, 'cnv', '*.bed')

    # Add all file entries to a dataframe
    data.append({
        'sample': sample_id,
        'sample_vcf': exp_sample_vcf[0] if exp_sample_vcf else None,
        'truth_vcf': truth_sample_vcf[0] if truth_sample_vcf else None,
        'variant_vcf': None,  # fill with appropriate value
        'json': exp_sample_json[0] if exp_sample_json else None,
        'truth_json': truth_sample_json[0] if truth_sample_json else None,
        'exon_cov': exp_sample_exoncov[0] if exp_sample_exoncov else None,
        'truth_exon_cov': truth_sample_exoncov[0] if truth_sample_exoncov else None,
        'cnv_bed': exp_sample_cnv[0] if exp_sample_cnv else None,
        'truth_cnv_bed': truth_sample_cnv[0] if truth_sample_cnv else None
    })

df = pd.DataFrame(data)

# Output dataframe to a csv based on output given
df.to_csv(output, index=False)
