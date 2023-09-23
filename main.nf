#! /usr/bin/env nextflow
//define channels from input file
Channel 
    .fromPath(params.inputlist)
    .ifEmpty {exit 1, "Cannot find input file : ${params.inputlist}"}
    .splitCsv(skip:1)
    .map{sample,pre_dbsnp_file, chr_1,chr_2, chr_3, chr_4, chr_5, chr_6,chr_7,chr_8,chr_9,chr_10,chr_11,chr_12,chr_13,chr_14,chr_15,chr_16,chr_17, chr_18,chr_19,chr_20,chr_21,chr_22, chr_X -> [sample, file(pre_dbsnp_file), file(chr_1),file(chr_2), file(chr_3), file(chr_4), file(chr_5), file(chr_6),file(chr_7),file(chr_8),file(chr_9),file(chr_10),file(chr_11),file(chr_12),file(chr_13),file(chr_14),file(chr_15),file(chr_16),file(chr_17), file(chr_18),file(chr_19),file(chr_20),file(chr_21),file(chr_22), file(chr_X)]}
    .set{ ch_input }


//run the script to make MTR input on above file paths
process  CloudOS_MTR_input{
    tag"$sample"
    publishDir "${params.outdir}/$sample", mode: 'copy'
    maxForks 900
    errorStrategy 'ignore'
    maxRetries 9
    
    input:
    set val(sample), file(pre_dbsnp_file), file(chr_1),file(chr_2), file(chr_3), file(chr_4), file(chr_5), file(chr_6),file(chr_7),file(chr_8),file(chr_9),file(chr_10),file(chr_11),file(chr_12),file(chr_13),file(chr_14),file(chr_15),file(chr_16),file(chr_17), file(chr_18),file(chr_19),file(chr_20),file(chr_21),file(chr_22), file(chr_X) from ch_input

    output:
    file "*_gistic_input.txt"
 
    script:
    """
    dbsnp_annotate.py -sample '$sample' -pre_dbsnp_file '$pre_dbsnp_file' -chr_1 '$chr_1' -chr_2 '$chr_2' -chr_3 '$chr_3' -chr_4 '$chr_4' -chr_5 '$chr_5' -chr_6 '$chr_6' -chr_7 '$chr_7' -chr_8 '$chr_8' -chr_9 '$chr_9' -chr_10 '$chr_10' -chr_11 '$chr_11' -chr_12 '$chr_12' -chr_13 '$chr_13' -chr_14 '$chr_14' -chr_15 '$chr_15' -chr_16 '$chr_16' -chr_17 '$chr_17' -chr_18 '$chr_18' -chr_19 '$chr_19' -chr_20 '$chr_20' -chr_21 '$chr_21' -chr_22 '$chr_22' -chr_X  '$chr_X'
    """ 
}
