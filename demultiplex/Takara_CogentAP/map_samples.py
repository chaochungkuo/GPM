#!/usr/bin/env python

import argparse
import pandas as pd

 
parser = argparse.ArgumentParser(

    add_help=False,

    description='''Use a 96-well sample layout file and a demux counts file to

        create a file that maps barcode to sample name.''',

    formatter_class=argparse.ArgumentDefaultsHelpFormatter,)

 

informational = parser.add_argument_group('informational arguments')

informational.add_argument(

    '-h',

    '--help',

    action='help',

    help='show this help message and exit',)

 
required = parser.add_argument_group('required arguments')

required.add_argument(

    '-c',

    '--counts',

    default=argparse.SUPPRESS,

    dest='counts_file',

    help='Demux counts file',)

required.add_argument(

    '-o',

    '--output',

    default=argparse.SUPPRESS,

    dest='barcode_sample_map_file',

    help='Barcode x Sample mapping output file',)

required.add_argument(

    '-s',

    '--sample',

    default=argparse.SUPPRESS,

    dest='sample_layout_file',

    help='Sample layout file on 96-well plate',)

 

optional = parser.add_argument_group('optional arguments')

optional.add_argument(

    '-b',

    '--bc1',

    default='BC1_layout.csv',

    dest='bc1_layout_file',

    help='BC1 layout file on 96-well plate',)

 

args = parser.parse_args()



sample_layout = pd.read_csv(args.sample_layout_file, index_col=0)


bc1_layout = pd.read_csv(args.bc1_layout_file, index_col=0)
bc1_stacked = bc1_layout.stack()
bc1_stacked = bc1_stacked.rename('BC1').reset_index()
bc1_stacked['Well position'] = bc1_stacked['level_0'] + bc1_stacked['level_1'].astype(str)



counts_column_names = ['Barcode', 'BarcodeID', 'Count']

counts = pd.read_csv(args.counts_file, names=counts_column_names)
 
counts['BC1'] = counts['BarcodeID'].str.split('_', n=1, expand=True).iloc[:, 1]



#sample_info = pd.DataFrame({'BC1':bc1_layout.stack(sort=False), 'Sample': sample_layout.stack(sort=False)})
sample_info = pd.merge(bc1_stacked, sample_layout, on = 'Well position', how ='left')

sample_info = pd.merge(sample_info, counts, on = 'BC1', how ='left')
#print(sample_info.head())



sample_info[['Well position',"Barcode",'BarcodeID','TSO Barcode', 'Sample Name']].to_csv(args.barcode_sample_map_file, index=False)
