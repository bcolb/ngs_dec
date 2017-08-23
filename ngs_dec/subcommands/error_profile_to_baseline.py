"""Turns a consolidated error profile into an error baseline.

Takes in the output from consolidate error profile and precomputes
data necessary to later establish p-values from experimental results.
"""

import logging
import argparse
import os
import sys
from collections import OrderedDict
import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

def build_parser(parser):
    ''' Parsers input/output arguments.'''
    parser.add_argument("input_file", 
                        help="error model file")
    parser.add_argument("output_file", 
                        help="file to write final error profile/baseline to")

def process_input(error_model, output_file):
    ''' Turns the error_model and snp_input files into dataframes,
    then calculates alpha and beta values for each. '''
    error_df = pd.read_csv(error_model)
    error_df['alpha'] = error_df.apply(calculate_alpha_value,axis=1)
    error_df['beta'] = error_df.apply(calculate_beta_value,axis=1)
    error_df.to_csv(output_file)
    
def calculate_alpha_value(row):
    ''' Calculates the alpha value for later use.'''
    mu = row['mean']
    std = row['std']
    var = std**2
    if var==0 or mu==0:
        return 0
    alpha = ((1-mu) / var - 1 / mu) * mu**2
    return alpha

def calculate_beta_value(row):
    ''' Calculates the beta value for later use.'''
    mu = row['mean']
    alpha = row['alpha']
    if mu == 0:
        return 0
    beta = alpha * (1/mu - 1)
    return beta

def action(args):
    ''' Takes in a baseline error file and outputs a version of the file with a new column for p values.'''
    error_model = args.input_file
    output_file = args.output_file

    if not os.path.isfile(error_model):
        sys.exit()

    process_input(error_model, output_file)
