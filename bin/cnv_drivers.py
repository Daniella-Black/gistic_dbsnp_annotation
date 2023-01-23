#!/usr/local/bin/python3

import pandas as pd
from dataclasses import dataclass
import argparse
#from os.path import exists

@dataclass
class SequenceRange:
    """Class for the start and end of a range."""
    name: str
    start: int
    end: int
    chrom: str
    def overlaps(self, other: "SequenceRange") -> bool:
        if self.chrom != other.chrom:
            return False
        return (other.start <= self.start <= other.end) or (other.start <= self.end <= other.end) or (self.start <= other.start <= self.end) or (self.start <= other.end <= self.end)

my_parser = argparse.ArgumentParser(description='find the missing data in cnv files')
my_parser.add_argument('-sample',
                       type=str,
                       help='sample')
my_parser.add_argument('-ploidy',
                       type=float,
                       help='ploidy')
my_parser.add_argument('-gene_df',
                       type=str,
                       help='path to the df of genes and cooridnates of canonical transcript')
my_parser.add_argument('-somatic_cnv_vcf',
                       type=str,
                       help='path to the mtr input cnv of the sample')
####data input

args = my_parser.parse_args()
sample = args.sample
ploidy = args.ploidy
gene_df_path = args.gene_df
cnv_path = args.somatic_cnv_vcf
    
    
amps = list()
missing_gene_data_sample= list()
gene_df = pd.read_csv(gene_df_path) 
gene_df=gene_df.dropna()

#set threshold for amplifications
if ploidy <2.5:
    amp_threshold = 5
elif ploidy >= 2.5:
    amp_threshold = 9

    
#file_exists = exists(cnv_path):
try:    
    cnv = pd.read_csv(cnv_path, '\t')
    cnv['total_cn'] = cnv['major_cn'] + cnv['minor_cn']
    total_cn = list(cnv['total_cn'])
    cnv['id'] = cnv['seqnames'].astype(str) + '_' + cnv['start'].astype(str) + '_' + cnv['end'].astype(str) + '_' + cnv['total_cn'].astype(str) +'_' + sample
    id_list = list(cnv['id'])
    for contig in range(len(total_cn)):
        if total_cn[contig] >= amp_threshold: 
            amps.append(id_list[contig]) 
    #take the list of amps obtained in for loop above and convert to a table
    if len(amps) >0:
        amps_df = pd.DataFrame(amps)
        amps_df[[ 'chr', 'start', 'end','total_cn', 'sample']] = amps_df[0].str.split('_', 4, expand=True)
        amps_df.drop(columns=[0])
    else:
        amps_df = pd.DataFrame(columns=[0])
    ##for each contig (no matter if it is amplified) report whether any contig overlaps with gene of interest - this is to identify samples with no data for the gene of interest for i in range(len(gene_df.index)): ##need file with all coding gene name chromosome coordinates #for gene in gene_df:
    genes_in_amps = [[] for _ in range(len(amps_df.index))]
    for i in range(len(gene_df.index)):  
        gene  = SequenceRange(gene_df['gene_name'][i], gene_df['start'][i], gene_df['end'][i], gene_df['chr'][i])
        ##find the genes with missing data
        contig_overlapping_gene = list()
        cnv_chr = cnv.loc[cnv['seqnames'].astype('str') == gene.chrom]
        cnv_chr.index = pd.RangeIndex(len(cnv_chr.index))
        for contig in range(len(list(cnv_chr['total_cn']))):
            contig_range= SequenceRange('place_holder', int(cnv_chr['start'][contig]), int(cnv_chr['end'][contig]), str(cnv_chr['seqnames'][contig]))
            if contig_range.overlaps(gene):
                contig_overlapping_gene.append(id_list[contig])
        if len(contig_overlapping_gene)== 0:
            missing_gene_data_sample.append(gene.name + '_' + sample)
        ##report the genes in each amp
        if len(amps) > 0:    
            for amp in range(len(amps_df.index)):
                amp_range= SequenceRange('place_holder', int(amps_df['start'][amp]), int(amps_df['end'][amp]), str(amps_df['chr'][amp]))
                if amp_range.overlaps(gene) and amp_range.chrom == gene.chrom:
                    genes_in_amps[amp].append(gene.name)
    amps_df['genes_in_amps'] = genes_in_amps
        
    #print(sample)
    with open(sample + '_mtr_format_cnv_missing.txt', 'w') as f:
        f.write(sample+' complete')
except FileNotFoundError as e:
    with open(sample + '_mtr_format_cnv_missing.txt', 'w') as f:
        f.write(sample+' no mtr format cnv file')
        
if len(missing_gene_data_sample) >0:
    missing_data_samples_gene_df = pd.DataFrame(missing_gene_data_sample)
    #missing_data_samples_gene_df = missing_data_samples_mdm2_df.rename(columns={0: 'missing_mdm2_samples'})
    missing_data_samples_gene_df[['sample', 'gene']] = missing_data_samples_gene_df[0].str.split('_', 1, expand=True)
    missing_data_samples_gene_df.drop(columns=[0])
else:
    missing_data_samples_gene_df = pd.DataFrame(columns=[0])

#output table of genes with missing data                                                                                                 
missing_data_samples_gene_df.to_csv(sample + '_genes_with_missing_data.csv')

#output amps_df
amps_df.to_csv(sample + '_amplifications.csv')
