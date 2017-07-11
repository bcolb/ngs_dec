"""Turns a varscan readcounts file into an error profile.

This command takes the output from varscan readcounts and 
processes it to create an error profile for use later in 
digital error correction.
"""

import logging
import os
import sys
from collections import OrderedDict

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('input_file',
                        help='Varscan readcounts input')
    parser.add_argument('output_file',
                        help='Output file to write to')

def parse_varscan_line(line, header_values):
    ''' Parses each individual line from varscan readcounts output and turns it into a dictionary.'''
    # chrom   position        ref_base        depth   q15_depth       base:reads:strands:avg_qual:map_qual:plus_reads:minus_reads
    line_values = line.strip().split('\t')
    vdict={}
    vdict['chrom'] = line_values[0]
    vdict['position'] = line_values[1]
    vdict['ref_base'] = line_values[2]
    vdict['depth'] = line_values[3]
    vdict['quality_depth'] = line_values[4] # varscan outputs q10_depth, q15_depth, etc. Just use 'quality depth' here to account for different thresholds
    # Line value 5 not used - would be reference base
    bases = OrderedDict()
    bases_list = ['A', 'T', 'C', 'G']
    empty_base_dict = {'reads':'0','strands':'0','avg_qual':'0','map_qual':'0','plus_reads':'0','minus_reads':'0'}
    for base in bases_list:
        bases[base] = dict(empty_base_dict)
    default_base = line_values[6] # should be T
    bases[default_base]['reads'] = line_values[7]
    bases[default_base]['strands'] = line_values[8]
    bases[default_base]['avg_qual'] = line_values[9]
    bases[default_base]['map_qual'] = line_values[10]
    bases[default_base]['plus_reads'] = line_values[11]
    bases[default_base]['minus_reads'] = line_values[12]
    if len(line_values) > 13:
        # parse variant values
        variant_list = line_values[13:]
        for variant in variant_list:
            variant_parts = variant.split(':')
            base = variant_parts[0]
            bases[base]={}
            bases[base]['reads'] = variant_parts[1]
            bases[base]['strands'] = variant_parts[2]
            bases[base]['avg_qual'] = variant_parts[3]
            bases[base]['map_qual'] = variant_parts[4]
            bases[base]['plus_reads'] = variant_parts[5]
            bases[base]['minus_reads'] = variant_parts[6]
    vdict['bases'] = bases
    return vdict

def varscan_dict_to_error_profile(vdict):
    ''' Turns the varscan dictionary into an error profile.'''
    lines = []
    for base in vdict['bases']:
        base_dict = vdict['bases'][base]
        if any(base_dict): # check that base_dict is not empty
            read_fraction = 0
            if int(vdict['quality_depth']) > 0:
                read_fraction = float(base_dict['reads']) / float(vdict['quality_depth'])
            # chrom,position,base,q15_depth,base_reads,base_read_fraction
            new_line = vdict['chrom'] + ',' + vdict['position'] + ',' + base + ',' + vdict['quality_depth'] + ',' + base_dict['reads'] + ',' + str(read_fraction) + '\n'
            lines += [new_line]
    return lines

def action(args):
    ''' Sets up processing of the input.'''
    input_file = args.input_file
    output_file = args.output_file
    if not os.path.isfile(input_file):
        print("Input file not found, exiting...")
        sys.exit()

    varscan_file = open(input_file, 'r')
    csv_output = open(output_file, 'w')
    csv_output.write('chrom,position,base,quality_depth,base_reads,base_read_fraction\n')

    lines = varscan_file.readlines() # may need to change this approach with file size
    header_values = lines[0].strip().split('\t')

    # break out last line of header values
    base_header_values = header_values[-1].split(':')
    header_values = header_values[:-1] + base_header_values

    for i in range(1, len(lines)):
        vdict = parse_varscan_line(lines[i], header_values)
        output_lines = varscan_dict_to_error_profile(vdict)
        for line in output_lines:
            csv_output.write(line)

    csv_output.close()
