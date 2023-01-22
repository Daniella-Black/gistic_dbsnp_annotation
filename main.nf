#! /usr/bin/env nextflow
//define channels from input file
Channel 
    .fromPath(params.inputlist)
    .ifEmpty {exit 1, "Cannot find input file : ${params.inputlist}"}
    .splitCsv(skip:1)
    .map{tumour_sample_platekey,somatic_cnv_vcf, ploidy, gene_df -> [tumour_sample_platekey, file(somatic_cnv_vcf), ploidy, file(gene_df)]}
    .set{ ch_input }


//run the script to make MTR input on above file paths
process  CloudOS_MTR_input{
    container = 'dockeraccountdani/pydocker:latest' 
    tag"$tumour_sample_platekey"
    publishDir "${params.outdir}/$tumour_sample_platekey", mode: 'copy'
    
    input:
    set val(tumour_sample_platekey), file(somatic_cnv_vcf), val(ploidy), file(gene_df) from ch_input

    output:
    file "*_genes_with_missing_data.csv"
    file "*_amplifications.csv"
    //file "readme.txt"
 
    script:
    """
    cnv_drivers.py -sample '$tumour_sample_platekey' -somatic_cnv_vcf '$somatic_cnv_vcf' -ploidy '$ploidy' -gene_df '$gene_df'
    """ 
}
