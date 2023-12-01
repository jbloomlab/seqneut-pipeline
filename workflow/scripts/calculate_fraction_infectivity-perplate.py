"""Converts barcode counts to a fraction infectivity measurement for each barcode"""

import os
import pandas as pd
import sys

# Write the log file to the specified location
sys.stderr = sys.stdout = log = open(snakemake.log[0], "w")

## ====== Inputs ====== ##
variants = snakemake.input.variants
standards = snakemake.input.standards
counts_csv = snakemake.input.counts

## ====== Output ====== ##
fraction_infectivity_csv = snakemake.output.fraction_infectivity
    
##Here, I am reading in a file that contains all of the barcode counts and a file that contains the linkage between the barcodes and the variants, such that we can link this information later on
variant_counts_samples = pd.read_csv(counts_csv)
variants_in_pool = pd.read_csv(variants).drop(columns=['library'])

#Remove counts for samples that are labelled to not be retained
variant_counts_samples = variant_counts_samples.query("retain")

#create a dataframe that merges the barcode counts with the variant names and has the amino acid information
variants_named_counts = pd.merge(variant_counts_samples, variants_in_pool, on="barcode", how="outer", validate="many_to_one")

#Pull the list of neutralization standard barcodes from the vaariant_in_pool file
neut_standard_barcodes = pd.read_csv(standards)
neut_standard_barcodes_list = neut_standard_barcodes['barcode'].tolist()

#Converting replicate to a string
variants_named_counts['replicate'] = "rep" + variants_named_counts['replicate'].astype(str)

#create a dataframe that is the count for the neut standards, summed over all of the barcodes for the neut standard
variants_named_counts_neut_standard = variants_named_counts.loc[variants_named_counts['barcode'].isin(neut_standard_barcodes_list)].groupby(by=['sample','replicate'], as_index=False).sum(numeric_only=True).rename(columns={"count":"neutstandard_count"}).drop(columns=['concentration','replicate','retain'])

#Add the column containing the neut standard summed value for each of the files to the dataframes as the neutstandard_count_pre column, and rename the other columns so that they are now correct for transformations
#This is done for each sample and barcode independently, I am using separate dataframes for the controls and the selections as I will want to group the controls by plate and then merge that information back into the selections dataframe later on
variants_named_counts_withneut = pd.merge(variants_named_counts, variants_named_counts_neut_standard, on="sample", how="inner", validate="many_to_one")

# We are going to drop the lines that are for the neutstandard barcodes from each file, as this will just get confusing
variants_named_counts_withneut = variants_named_counts_withneut.loc[~variants_named_counts_withneut['barcode'].str.contains('|'.join(neut_standard_barcodes_list))]

# Create a dataframe that is the count for the neut standards, summed over all of the barcodes for the neut standard
variants_named_counts_neut_standard = variants_named_counts.loc[variants_named_counts['barcode'].isin(neut_standard_barcodes_list)].groupby(by=['sample','replicate'], as_index=False).sum(numeric_only=True).rename(columns={"count":"neutstandard_count"}).drop(columns=['concentration','replicate','retain'])

# Create dataframes for selections and controls, use names that can be used every time, and rename columns
variants_named_counts_selections_withneut = variants_named_counts_withneut.loc[~variants_named_counts_withneut['sample'].str.contains('Noselection|CellsOnly|Noserum')].rename(columns={"sample": "serum_sample", "count": "postselection_count", "neutstandard_count":"neutstandard_count_post"}).drop(columns=['library'])
variants_named_counts_controls_withneut = variants_named_counts_withneut.loc[variants_named_counts_withneut['sample'].str.contains('Noselection|Noserum')].rename(columns={"sample": "no-serum_sample", "count": "preselection_count", "neutstandard_count":"neutstandard_count_pre"}).drop(columns=['library'])

### As we have multiple control wells per plate, we need to merge these replicate control wells into a single ratio of each barcode to the neut standard barcodes for the whole plate
#### To do this, I am creating a dataframe where I first calculate the ratio of counts of each barcode to the sum of all neutralization standard barcodes for each of the control wells and then take the median value of this for all wells on the plates, I can use the plate name control for normalization

#Divide the preselection_count by the neutralization standard count for each well for each barcode
variants_named_counts_controls_withneut['preselection_count_normalized'] = variants_named_counts_controls_withneut['preselection_count'] / variants_named_counts_controls_withneut['neutstandard_count_pre']

# Remove any samples where neut standard counts is too low
variants_named_counts_controls_withneut = variants_named_counts_controls_withneut.loc[variants_named_counts_controls_withneut['neutstandard_count_pre']>1000]

# I am currently taking the median of the normalized count for all the control wells for each barcode count/the summed counts for the neut-standard barcodes for each plate
variants_named_counts_controls_withneut_median = variants_named_counts_controls_withneut.groupby(by=['barcode','plate'], as_index=False).median(numeric_only=True)

### Now we want to merge the dataframes so that we can link the pre-selection counts for each barcode and neut-standard with the post selection counts for each barcode and neut standards

#We are just dividing the count for the barcode by the summed counts for the neutstandard
variants_named_counts_selections_withneut['postselection_count_normalized'] = variants_named_counts_selections_withneut['postselection_count'] / variants_named_counts_selections_withneut['neutstandard_count_post'] 

### Now, map the preselection counts with the post selection counts:

#Merge selections and controls on barcode and dated_controls
variant_counts_samples_mapped = pd.merge(variants_named_counts_selections_withneut, variants_named_counts_controls_withneut_median, on=["barcode","plate"],how="left", validate="many_to_one")
variant_counts_withpreselection = variant_counts_samples_mapped.loc[variant_counts_samples_mapped['preselection_count'] > 10]

### Calculate the normalized count for each barcoded variant by dividing the count preselection by the count post selection
variant_counts_withpreselection['normalized_count'] = variant_counts_withpreselection["postselection_count_normalized"].div(variant_counts_withpreselection["preselection_count_normalized"].values)

#Remove unnecessary columns to get the simplified dataframe
variant_counts_normalized = variant_counts_withpreselection.drop(columns = ['postselection_count','neutstandard_count_post','postselection_count_normalized','neutstandard_count_pre','preselection_count','preselection_count_normalized','retain_x','retain_y','standard_set'])

#Add columns which correspond to serum concentration from dilution factor
variant_counts_normalized = variant_counts_normalized.rename(columns={"concentration_x": "serum_dilution"})
variant_counts_normalized['serum_dilution'] = pd.to_numeric(variant_counts_normalized['serum_dilution'], errors='ignore')
variant_counts_normalized['serum_concentration'] = 1/variant_counts_normalized['serum_dilution']

#rename columns for loading into neutcurve to fit Hill curve and calculate NT50
fractioninfectivity_df = variant_counts_normalized[['serum','individual','condition','barcode','serum_concentration','normalized_count','strain','serum_sample','replicate']].copy()
fractioninfectivity_df['virus'] = fractioninfectivity_df['barcode']+"_"+fractioninfectivity_df['replicate'].astype(str)
fractioninfectivity_df = fractioninfectivity_df.rename(columns={"serum_concentration":"concentration","normalized_count":"fraction infectivity","serum_sample":"sample"})

#Save fraction infectivity to file
#os.mkdir(config["selection_dir"])
fractioninfectivity_df.to_csv(fraction_infectivity_csv)
