#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * ========================================================================================
 * CNV ANALYSIS PIPELINE
 * ========================================================================================
 */


/*
 * ========================================================================================
 * Define workflow parameters
 * ========================================================================================
 */
params.reads = "$baseDir/data/raw/*_{1,2}.fastq.gz"
params.reference_dir = "$baseDir/data/reference"
params.outdir = "$baseDir/results"
params.help = false

// Help message
if (params.help){
    log.info """
    CNV Analysis Pipeline
    =====================
    Usage:
    nextflow run main.nf --reads '<path to your reads>' --reference_dir '<path to reference dir>' --outdir '<path to output dir>'

    Optional arguments:
      --reads             : Path to input FASTQ files (glob pattern) (default: ${params.reads})
      --reference_dir     : Path to directory containing reference genome (must be a .fasta file)
      --outdir            : Path to output results directory (default: ${params.outdir})
      --help              : Display this help message
    """.stripIndent()
    exit 0
}


/*
 * ========================================================================================
 * Define Channels
 * ========================================================================================
 */
Channel
    .fromFilePairs(params.reads)
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!" }
    .set { ch_reads }

Channel
    .fromPath(params.reference_dir)
    .set { ch_ref_dir }


/*
 * ========================================================================================
 * W O R K F L O W
 * ========================================================================================
 */
workflow {
    // Step 1: Quality Control
    FASTQC ( ch_reads )

    // Step 2: Trim Reads
    TRIM_GALORE ( ch_reads )

    // Step 3: Index Reference Genome (runs once)
    BWA_INDEX ( ch_ref_dir )

    // Step 4: Align Reads to Reference
    BWA_ALIGN ( TRIM_GALORE.out.trimmed_reads, BWA_INDEX.out.index )

    // Step 5: Process SAM/BAM files
    SAMTOOLS_CALL ( BWA_ALIGN.out.sam )

    // Step 6: Run CNVnator for CNV analysis
    CNV_ANALYSIS ( SAMTOOLS_CALL.out.bam, ch_ref_dir )
}


/*
 * ========================================================================================
 * P R O C E S S E S
 * ========================================================================================
 */

process FASTQC {
    tag "$pair_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(pair_id), path(reads)

    output:
    path "fastqc_out/*"

    script:
    """
    mkdir fastqc_out
    fastqc -o fastqc_out -t $task.cpus $reads
    """
}

process TRIM_GALORE {
    tag "$pair_id"
    publishDir "${params.outdir}/trimmed_reads", mode: 'copy'

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("${pair_id}*") , emit: trimmed_reads

    script:
    def (read1, read2) = reads
    """
    trim_galore --paired --cores $task.cpus -o . "$read1" "$read2"
    """
}

process BWA_INDEX {
    // This process will only run once due to the nature of the input channel
    publishDir "${params.outdir}/reference_index", mode: 'copy'

    input:
    path ref_dir

    output:
    path "indexed_ref/*" , emit: index

    script:
    """
    mkdir indexed_ref
    cp ${ref_dir}/*.fasta indexed_ref/
    bwa index indexed_ref/*.fasta
    """
}

process BWA_ALIGN {
    tag "$pair_id"
    publishDir "${params.outdir}/aligned", mode: 'copy'

    input:
    tuple val(pair_id), path(reads)
    path index_dir

    output:
    tuple val(pair_id), path("*.sam"), emit: sam

    script:
    def (read1, read2) = reads
    def reference = file(index_dir.list().find { it.name.endsWith('.fasta') })
    """
    bwa mem -t $task.cpus $reference $read1 $read2 > ${pair_id}.sam
    """
}


process SAMTOOLS_CALL {
    tag "$pair_id"
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    tuple val(pair_id), path(sam_file)

    output:
    tuple val(pair_id), path("*.sorted.bam"), emit: bam

    script:
    """
    samtools view -S -b $sam_file | samtools sort -o ${pair_id}.sorted.bam
    samtools index ${pair_id}.sorted.bam
    """
}

process CNV_ANALYSIS {
    tag "$pair_id"
    publishDir "${params.outdir}/cnv", mode: 'copy'

    input:
    tuple val(pair_id), path(bam_file)
    path ref_dir

    output:
    path "*.cnv"

    script:
    def reference = file(ref_dir.list().find { it.name.endsWith('.fasta') })
    """
    cnvnator -root ${pair_id}.root -tree $bam_file
    cnvnator -root ${pair_id}.root -his 10 -d $reference
    cnvnator -root ${pair_id}.root -stat 10
    cnvnator -root ${pair_id}.root -partition 10
    cnvnator -root ${pair_id}.root -call 10 > ${pair_id}.cnv
    """
}
