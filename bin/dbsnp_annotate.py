#!/usr/local/bin/python3

import pandas as pd
import os

my_parser = argparse.ArgumentParser(description='find the missing data in cnv files')
my_parser.add_argument('-sample',
                       type=str,
                       help='sample')
my_parser.add_argument('-pre_dbsnp_file',
                       type=float,
                       help='the pre-dbsnp annotated input')
my_parser.add_argument('-chr_1',
                       type=str,
                       help='dbsnp per chromosome')
my_parser.add_argument('-chr_2',
                       type=str,
                       help='dbsnp per chromosome')
my_parser.add_argument('-chr_3',
                       type=str,
                       help='dbsnp per chromosome')
my_parser.add_argument('-chr_4',
                       type=str,
                       help='dbsnp per chromosome')
my_parser.add_argument('-chr_5',
                       type=str,
                       help='dbsnp per chromosome')
my_parser.add_argument('-chr_6',
                       type=str,
                       help='dbsnp per chromosome')
my_parser.add_argument('-chr_7',
                       type=str,
                       help='dbsnp per chromosome')
my_parser.add_argument('-chr_8',
                       type=str,
                       help='dbsnp per chromosome')
my_parser.add_argument('-chr_9',
                       type=str,
                       help='dbsnp per chromosome')
my_parser.add_argument('-chr_10',
                       type=str,
                       help='dbsnp per chromosome')
my_parser.add_argument('-chr_11',
                       type=str,
                       help='dbsnp per chromosome')
my_parser.add_argument('-chr_12',
                       type=str,
                       help='dbsnp per chromosome')
my_parser.add_argument('-chr_13',
                       type=str,
                       help='dbsnp per chromosome')
my_parser.add_argument('-chr_14',
                       type=str,
                       help='dbsnp per chromosome')
my_parser.add_argument('-chr_15',
                       type=str,
                       help='dbsnp per chromosome')
my_parser.add_argument('-chr_16',
                       type=str,
                       help='dbsnp per chromosome')
my_parser.add_argument('-chr_17',
                       type=str,
                       help='dbsnp per chromosome')
my_parser.add_argument('-chr_18',
                       type=str,
                       help='dbsnp per chromosome')
my_parser.add_argument('-chr_19',
                       type=str,
                       help='dbsnp per chromosome')
my_parser.add_argument('-chr_20',
                       type=str,
                       help='dbsnp per chromosome')
my_parser.add_argument('-chr_21',
                       type=str,
                       help='dbsnp per chromosome')
my_parser.add_argument('-chr_22',
                       type=str,
                       help='dbsnp per chromosome')
my_parser.add_argument('-chr_X',
                       type=str,
                       help='dbsnp per chromosome')

args = my_parser.parse_args()
sample = args.sample
gistc_input = args.pre_dbsnp_file
chr_1 = args.chr_1
chr_2 = args.chr_2
chr_3 = args.chr_3
chr_4 = args.chr_4
chr_5 = args.chr_5
chr_6 = args.chr_6
chr_7 = args.chr_7
chr_8 = args.chr_8
chr_9 = args.chr_9
chr_10 = args.chr_10
chr_11 = args.chr_11
chr_12 = args.chr_12
chr_13 = args.chr_13
chr_14 = args.chr_14
chr_15 = args.chr_15
chr_16 = args.chr_16
chr_17 = args.chr_17
chr_18 = args.chr_18
chr_19 = args.chr_19
chr_20 = args.chr_20
chr_21 = args.chr_21
chr_22 = args.chr_22
chr_X = args.chr_X

#read in all the per chromosome dbsnps and put the pos of variants into a list

##define the funciton which will get the lists and with the for loop below put them into the snps dictionary
def snp_dict(chr_tab):
    dbsnp_chrom = pd.read_csv(chr_tab,sep='\t',header=None)
    return list(dbsnp_chrom[1])

chr_dict = {'1':chr_1,'2': chr_2, '3': chr_3, '4':chr_4, '5':chr_5, '6':chr_6,'7':chr_7, '8':chr_8, '9':chr_9, '10': chr_10,'11':chr_11, '12': chr_12, '13':chr_13, '14':chr_14, '15':chr_15, '16': chr_16, '17':chr_17, '18':chr_18,'19':chr_19, '20':chr_20, '21':chr_21, '22':chr_22, 'X':chr_X}
snps = {}
for key in chr_dict.keys():
    snps[key] = snp_dict(chr_dict[key])
    

##annotate the existing gistic inpput file with the number of dbsnp common snps occur in each segment
gistc_input.index = gistc_input.reset_index(drop=True)
num_markers = list()
for row in range(len(gistc_input.index)):
    min = gistc_input[2][row]
    max = gistc_input[3][row]
    chr = gistc_input[1][row]
    tmp = snps[chr]
    num_markers.append(len(set([i for i in tmp if (i >= min) and (i <= max)])))
    
gistc_input['Num Markers'] = num_markers    

#output
gistc_input.to_csv(sample + '_gistic_input.txt')
