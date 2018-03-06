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
params.reads         = "$baseDir/fastq/*_{1,2}.fastq.gz"
params.output        = "results/"


log.info "C A L L - N F -  E X A M P L E  ~  version 1.0"
log.info "====================================="
log.info "name                   : ${params.name}"
log.info "reads                  : ${params.reads}"
log.info "fasta                  : ${params.fasta}"
log.info "output                 : ${params.output}"
log.info "\n"


/*
 * Input parameters validation
 */

fasta_file = file(params.fasta)

/*
 * validate input files
 */
if( !fasta_file.exists() ) exit 1, "Missing reference file: ${fasta_file}"


/*
 * Create a channel for read files 
 */
 
Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_files } 


process index {
     tag "$fasta_file.simpleName"
    input:
    file fasta_file
    
    output:
    file outfile
    file("*") into bwa_index
    file 'fasta.dict' into fasta_dict
      
    script:
    //
    // bwa  index
    //
    """
    samtools faidx  ${fasta_file} && picard CreateSequenceDictionary R=${fasta_file} o=fasta.dict && bwa index ${fasta_file} > outfile
    """
}


process mapping {
    tag "$name"
     publishDir  params.output, mode:'copy' 
    input:
    file("*") from bwa_index
    set val(name), file(reads) from read_files


    output:
    file("${name}.bam") into aling_dir
    file("stats/${name}.stats") into aling_stats


    script:
    //
    // bwa align
    //
    def single = reads instanceof Path
    if( !single ) {
        """
       mkdir -p stats; bwa mem -M   ${fasta_file} -t ${task.cpus}  ${reads[0]} ${reads[1]}|samblaster -M | samtools fixmate - - | samtools  sort -O bam -o ${name}.bam -;samtools index  ${name}.bam;samtools stats ${name}.bam | grep ^SN | cut -f 2- > stats/${name}.stats 
        """
    }  
    else {
        """
        
        bwa mem -M  -b -i ${fasta_file} -t ${task.cpus} -o ${name} ${reads[0]} 
        """
    }

}




process vcf {
       
    publishDir  params.output, mode:'copy' 
    input:
    file('sew') from aling_dir.collect()
    file fasta_file
    output:
    file('all.raw.vcf') 
   

    script:
    //
    // vcf
    //
    """
    freebayes -f ${fasta_file} ${sew}  > all.raw.vcf

    """
}










workflow.onComplete {
        println ( workflow.success ? "\nDone! Open the following vcf file --> $params.output/all.raw.vcf\n" : "Oops .. something went wrong" )
}

