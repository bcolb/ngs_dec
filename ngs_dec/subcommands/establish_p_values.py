"""Calculates P values from a given variants file and baseline.

Takes in a varscan SNP file and a baseline file
and in turn calculates P-values for each variant.
"""
import logging
import os
import sys
import ipdb
import numpy as np
import pandas as pd
from scipy import stats
import vcf

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('sample_variants',
                        help='Variants file for a sample')
    parser.add_argument('error_baseline',
                        help='An error profile baseline file to analyze against')
    parser.add_argument('output_file',
                        help='Output file to write to')
    parser.add_argument('--full_output', action='store_true',
                        help="Output full merged data from error profile and variant file as csv")

def action(args):
    sample_variants = args.sample_variants
    error_baseline = args.error_baseline
    output_file = args.output_file
    full_output = args.full_output
    
    # Check that files exist
    if not os.path.isfile(sample_variants):
        print("Sample variants file not found. Exiting.")
        sys.exit()
    if not os.path.isfile(error_baseline):
        print("Error baseline file not found. Exiting.")
        sys.exit()

    # build pandas dataframe for sample from VCF
    columns = ['position_base','chromosome','position','var_base','ref_base','var_freq_pct','var_freq_flt','var_pval']
    vcf_reader = vcf.Reader(open(sample_variants))
    row_list = []
    for record in vcf_reader:
        call_data = record.samples[0].data
        # Need following data: Chromosome, Position, Var Base, Ref Base, Var Freq
        row_dict = {}
        row_dict['chromosome'] = record.CHROM
        row_dict['position'] = record.POS
        if len(record.ALT) == 1:
            row_dict['var_base'] = record.ALT[0]
        else:
            print("Multiple variant bases...not supported - passing")
            continue
        row_dict['ref_base'] = record.REF
        row_dict['var_freq_pct'] = call_data.FREQ
        row_dict['var_freq_flt'] = percent_to_float(call_data.FREQ)
        row_dict['var_pval'] = call_data.PVAL
        row_dict['position_base'] = str(row_dict['chromosome']) + ':' + str(row_dict['position']) + ':' + str(row_dict['var_base'])
        row_list.append(row_dict)
    sample_df = pd.DataFrame(row_list, columns=columns)

    # create another dataframe from the error baseline
    dtype={'chrom':str,'position':int,'base':str,'mean':float,'std':float,'alpha':float,'beta':float}
    #error_df = pd.DataFrame.from_csv(error_baseline, dtype=dtype)
    error_df = pd.read_csv(error_baseline, dtype=dtype)
    error_df['position_base'] = error_df.apply(get_error_position_base_key,axis=1)
    error_df.set_index('position_base')

    # do left join with sample_df as left
    #joined_df = sample_df.join(error_df, on='position_base', how='left', lsuffix='sample', rsuffix='error')
    merged_df = sample_df.merge(error_df, on='position_base',how='inner')

    # uw_dec_p_value calculated using beta distribution parameters
    merged_df['uw_dec_p_value'] = merged_df.apply(calculate_p_value,axis=1)

    # used in full output only
    merged_df['pdf'] = merged_df.apply(calculate_pdf,axis=1)
    merged_df['cdf'] = merged_df.apply(calculate_cdf,axis=1)


    # drop rows with null results
    final_df = merged_df[merged_df['uw_dec_p_value'].isnull()==False]

    # Write output file
    if full_output:
        merged_df_to_full_output(final_df, output_file)
    else:
        merged_df_to_munge_ready_output(final_df, output_file)

def get_sample_position_base_key(row):
    ''' Creates a position key consisting of the chromosome, position, and base'''
    return str(row['Chrom']) + ":" + str(row['Position']) + ":" + str(row['Var'])

def get_error_position_base_key(row):
    ''' Creates a position key consisting of the chromosome, position, and base'''
    return str(row['chrom']) + ":" + str(row['position']) + ":" + str(row['base'])

def calculate_pdf(row):
    ''' Calculates the probability density function given beta values.'''
    if row['var_freq_flt'] and row['alpha'] and row['beta']:
        var_freq_flt = row['var_freq_flt']
        if var_freq_flt == 1.0:
            var_freq_flt = var_freq_flt - 0.0001 #0.9999
        return stats.beta.pdf(var_freq_flt, row['alpha'], row['beta'])
    else:
        return None

def calculate_cdf(row):
    ''' Calculates the cumulative density function given beta distribution parameters.'''
    if row['var_freq_flt'] and row['alpha'] and row['beta']:
        var_freq_flt = row['var_freq_flt']
        if var_freq_flt == 1.0:
            var_freq_flt = var_freq_flt - 0.0001 #0.9999
        return stats.beta.cdf(var_freq_flt, row['alpha'], row['beta'])
    else:
        return None

def calculate_p_value(row):
    """ Calculates the cdf and inverts it."""
    if row['var_freq_flt'] and row['alpha'] and row['beta']:
        var_freq_flt = row['var_freq_flt']
        if var_freq_flt == 1.0:
            var_freq_flt = var_freq_flt - 0.0001 #0.9999
        return (1 - stats.beta.cdf(var_freq_flt, row['alpha'], row['beta']))
    else:
        return None

def percent_to_float(val):
    '''turns a percentage to a float.'''
    return float(val.strip('%'))/100

def merged_df_to_full_output(merged_df, output_file):
    '''output full csv'''
    merged_df.to_csv(output_file)

def merged_df_to_munge_ready_output(merged_df, output_file):
    ''' Just need to output:
    generic\tp.value\tchr\tstart\tstop\tref_base\tvar_base\n
    '''
    out_f = open(output_file, 'w')
    # loop through dataframe, write row to file
    for index, row in merged_df.iterrows():
        line = "generic\t" + str(row['uw_DEC_probability']) + "\t" + str(row['chrom']) + "\t" + str(row['position_x']) + "\t" + str(row['position_x']) + "\t" + str(row['ref_base']) + "\t" + str(row['var_base']) + "\n"
        out_f.write(line)
    out_f.close()
