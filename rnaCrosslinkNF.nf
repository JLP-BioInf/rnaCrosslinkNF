#!/usr/bin/env nextflow

/*
 * Description: This pipeline processes raw RNA-Seq data and performs quality control,
 *              adapter trimming, read alignment, and quantification.
 * Steps:
 *   1. Quality control using FastQC
 *   2. Adapter trimming using Cutadapt
 *   3. Read alignment using STAR
 *   4. Quantification using featureCounts
 *   5. MultiQC report generation
 *
 *   Author - Jonathan price
 */


// enable bianaries for extra functionality
nextflow.enable.dsl=2
nextflow.enable.moduleBinaries = true


// inncludes single
include { pear; gunzip; cutadapt; fastqc; multiqc ; starIndex; tstk} from './rnaCrossPre'
include {samtoolsView; align; samtools_idxstats; processSAM;cutadaptr4} from './rnaCrossPre'


// inncludes repeats
include { fastqc as fastqc2 } from './rnaCrossPre'
include { fastqc as fastqc3 } from './rnaCrossPre'


// Parameters
params.reads = "./data/*neb*_R{1,2}.fastq.gz"
params.transcriptome = "./Transcripts.fasta"
params.outdir = "$baseDir/results_myc_neb"
params.sampleTable= "$baseDir/data/sampleTable.txt"










workflow {

    // Read channel and STAR Index --------
    read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true )
    starIndexCh = starIndex(params.transcriptome)


    // Fastqc of Raw          --------
    f1ch = fastqc(read_pairs_ch,Channel.value("raw"))


    // Trim Adapters                 --------
    trimmedreads =   cutadaptr4(read_pairs_ch)
    trimmedreads2 =   cutadapt(read_pairs_ch)


    // Fastqc of Trimmed      --------
    f3ch = fastqc2(trimmedreads2,Channel.value("trim"))


    // Assemble Reads             --------
    assembledReads = pear( gunzip(trimmedreads2))
    f4ch = fastqc3(assembledReads,Channel.value("assembled"))


    // assemble reads
    cleanReads = tstk( assembledReads)



    // Align the clean reads                   --------
    alignedreads1 = align(starIndexCh,cleanReads,Channel.value(3))


    // Make sam files         --------
    view = samtoolsView(alignedreads1,assembledReads)
    processSAM(view)


    // Make stats with samtools  --------
    samtools = samtools_idxstats(alignedreads1)


    // MultiQC         --------
    multiqc(f1ch.mix(f3ch,f4ch,alignedreads1).collect()) //fix the input channel
}
