#!/usr/local/bin/python3

import pandas as pd
from dataclasses import dataclass

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

sample = args[1]
ploidy = args[2]
gene_df_path = args[3]
cnv_path= args[4]

ploidy_samples = {sample:ploidy}



amps = list()
missing_gene_data_sample= list()
gene_df = pd.read_csv(gene_df_path) 
gene_df=gene_df.dropna()

#set threshold for amplifications
if ploidy_samples [sample]['Ploidy'] <2.5:
    amp_threshold = 5
elif ploidy_samples [sample]['Ploidy'] >= 2.5:
    amp_threshold = 9

cnv = pd.read_csv(cnv_path, '\t')
cnv['total_cn'] = cnv['major_cn'] + cnv['minor_cn']
total_cn = list(cnv['total_cn'])
cnv['id'] = cnv['seqnames'].astype(str) + '_' + cnv['start'].astype(str) + '_' + cnv['end'].astype(str) + '_' + cnv['total_cn'].astype(str) +'_' + sample
id_list = list(cnv['id'])
for contig in range(len(total_cn)):
    if total_cn[contig] >= amp_threshold: 
        amps.append(id_list[contig]) 
##for each contig (no matter if it is amplified) report whether any contig overlaps with gene of interest - this is to identify samples with no data for the gene of interest for i in range(len(gene_df.index)): ##need file with all coding gene name chromosome coordinates #for gene in gene_df:
for i in range(len(gene_df.index)):  
    gene  = SequenceRange(gene_df['gene_name'][i], gene_df['start'][i], gene_df['end'][i], gene_df['chr'][i])
    contig_overlapping_gene = list()
    cnv_chr = cnv.loc[cnv['seqnames'].astype('str') == gene.chrom]
    cnv_chr.index = pd.RangeIndex(len(cnv_chr.index))
    for contig in range(len(list (cnv_chr['total_cn']))):
        contig_range= SequenceRange('place_holder', int(cnv_chr['start'][contig]), int(cnv_chr['end'][contig]), str(cnv_chr['seqnames'][contig]))
        if contig_range.overlaps(gene):
            contig_overlapping_gene.append(id_list[contig])
    if len(contig_overlapping_gene)== 0:
        missing_gene_data_sample.append(gene.name + '_' + sample)
#print(sample)
missing_data_samples_gene_df = pd.DataFrame(missing_gene_data_sample)
#missing_data_samples_gene_df = missing_data_samples_mdm2_df.rename (columns={0: 'missing_mdm2_samples'}) missing_data_samples_gene_df [['sample', 'gene']] = missing_data_samples_gene_df[0].str.split(', 1, expand=True)
missing_data_samples_gene_df.drop(columns=[0])
missing_data_samples_gene_df.to_csv('genes_with_missing_data.csv')
#take the list of amps obtained in for loop above and convert to a table
amps_df = pd.DataFrame(amps)
amps_df[[ 'chr', 'start', 'end','total_cn', 'sample']] = amps_df[0].str.split('_', 4, expand=True)
amps_df.drop(columns=[0])
amps_df.to_csv('amplifications.csv')
