""" Consolidates multiple error profiles into a single file.

Takes multiple error profiles in and computes the mean and 
standard deviation for each base at each position.
"""

import logging
import os
import sys
from collections import OrderedDict
import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('variant_profile_directory',
                        help='directory with variant profile csv files')
    parser.add_argument('output_file',
                        help='File to write profile to')

def action(args):
    vp_dir = args.variant_profile_directory
    output_file = args.output_file

    # find all variant_profile.csv files
    vp_files = []
    for item in os.listdir(vp_dir):
        # descend into sample subdirectories
        item_dir = os.path.join(vp_dir, item)
        if os.path.isdir(item_dir):
            for sub_item in os.listdir(item_dir):
                # check if ends with variant_profile.csv
                sub_item_path = os.path.join(item_dir, sub_item)
                if os.path.isfile(sub_item_path) and sub_item.endswith('variant_profile2.csv'):
                    vp_files.append(sub_item_path)

    # Turn into pandas dataframes and concatenate together
    df_list = []
    header_vals = ['chrom','position','base','quality_depth','base_reads','base_read_fraction']
    dtype = {'chrom':str,'position':np.int64,'base':str,'quality_depth':np.int64,'base_reads':np.int64,'base_read_fraction':np.float64}
    for i in range(0, len(vp_files)):
        f = vp_files[i]
        df = pd.read_csv(f, dtype=dtype) # , dtype=dtype) #,dtype=dtype) #,names=header_vals)
        df_list.append(df)
    
    master_df = pd.concat(df_list).reset_index()

    # in master df group by chrom,position,base
    grouped_df = master_df.groupby(['chrom','position','base'])
    results_df = pd.DataFrame()
    results_df['mean'] = grouped_df['base_read_fraction'].mean()
    results_df['std'] = grouped_df['base_read_fraction'].std()
    
    # write to output
    results_df.reset_index().to_csv(output_file, index=False)

