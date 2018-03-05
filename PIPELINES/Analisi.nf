/*
 * Copyright (c) 2018-2025, .
 *
 *   This file is part of 'variant calling'.
 *
 *   Variant calling using frebays software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 */

/* 
 * Main CALL-nf pipeline script
 *
 * @authors
 * Maurizio Polano <mauriziopolano@blu.it>
 */


params.fasta = "$baseDir/reference/human_g1k_v37_MT.fasta"
params.name          = "CALL-NF"
params.reads         = "$baseDir/fastq/*.gz"
params.output        = "results/"


log.info "C A L L - N F -  E X A M P L E  ~  version 1.0"
log.info "====================================="
log.info "name                   : ${params.name}"
log.info "reads                  : ${params.reads}"
log.info "fasta           : ${params.fasta}"
log.info "output                 : ${params.output}"
log.info "\n"


/*
 * Input parameters validation
 */

fasta_file     = file(prams.fasta)

/*
 * validate input files
 */
if( !fasta_file.exists() ) exit 1, "Missing reference file: ${fasta_file}"


/*
 * Create a channel for read files 
 */
 
Channel
    .fromFilePairs( params.reads, size: -1 )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_files } 


process index {
    input:
    file fasta_file
    
    output:
    file "fasta.index" into fasta_index
    file "fasta.dict" into fasta_dict
      
    script:
    //
    // bwa  index
    //
    """
    samtools faidx  ${fasta_file} && picard CreateSequenceDictionary R=${fasta_file} o=${fasta_dict} && bwa index ${fasta_file}"
    """
}


process mapping {
    tag "reads: $name"

    input:
    file index from fasta_index
    set val(name), file(reads) from read_files

    output:
    file "align_${name}" into align_out_dirs 

    script:
    //
    // bwa align
    //
    def single = reads instanceof Path
    if( !single ) {
        """
        mkdir align_${name}
        bwa mem -M  -b ${params.bootstrap} -i ${index} -t ${task.cpus} -o kallisto_${name} ${reads}
        """
    }  
    else {
        """
        mkdir align${name}
        bwa mem -M -R '@RG\tID:ind1\tSM:ind1\tLB:ind1\tPU:ind1\tPL:Illumina' reference/human_g1k_v37_MT.fasta fastq/ind1_1.fastq.gz fastq/ind1_2.fastq.gz
        """
    }

}



