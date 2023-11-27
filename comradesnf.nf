#!/usr/bin/env nextflow

/*
 * Description: This pipeline processes raw RNA-Seq data and performs quality control,
 *              adapter trimming, read alignment and parsing for final output for comradesOO
 * Steps:
 *   1. Quality control using FastQC
 *   2. Adapter trimming using Cutadapt
 *   3. Combining paired-end reads
 *   4. Collapsing unique reads based on UMI
 *   2. UMI trimming using Cutadapt
 *   3. Read alignment using STAR
 *   4. SAM alignment parsing for final output
 *   5. MultiQC report generation
 *
 *   Author - Jonathan price
 */


// Input and Output parameters
params.reads         = "$baseDir/data/*_R{1,2}.fastq.gz"
params.transcriptome = "$baseDir/data/18Sref.fasta"
params.outdir        = "$baseDir/results"
params.sampleTable   = "$baseDir/data/sampleTable.txt"


workflow {
    read_pairs_ch  = channel.fromFilePairs( params.reads,
                                            checkIfExists: true )
    starIndexCh    = STARINDEX( params.transcriptome )
    f1ch           = FASTQC( read_pairs_ch )
    trimmedreads   = CUTADAPTONE( read_pairs_ch )
    f2ch           = FASTQC2( trimmedreads )
    assembledReads = PEAR( GUNZIP( trimmedreads ) )
    f3ch           = FASTQC3(assembledReads )
    collapsedReads = TSTK( assembledReads )
    cleanReads     = CUTADAPTWO( collapsedReads )
    alignedreads   = ALIGN( starIndexCh,
                            cleanReads)
    sam            = SAMTOOLSVIEW ( alignedreads,
                                    collapsedReads)
    PROCESSSAM( sam )
    multiqc( f2ch.mix( f1ch,
                       f3ch,
                     alignedreads ).collect( ) )
}



process GUNZIP {
    cache 'deep'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id) , path("${sample_id}*.fastq")

    script:

    reads_1_fq = reads[0].name.split("\\.")[0] + '.fastq'
    reads_2_fq = reads[1].name.split("\\.")[0] + '.fastq'

    """

    gzip -dc ${reads[0]} > ${reads_1_fq}
    gzip -dc ${reads[1]} > ${reads_2_fq}
    """
}

process PEAR {
    cache 'deep'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id) , path("${sample_id}*assembled.fastq")

    script:
    """

    pear  \
        --min-assembly-length 20 \
        -f ${reads[0]}  \
        -r ${reads[1]}   \
        -o $sample_id \
        -v 7 \
        -j 5
    """
}


process TSTK {
    cache 'deep'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id) , path("${sample_id}*assembled.tstk.fasta")

    script:
        output = sample_id + 'assembled.tstk.fasta'
    """

    python3 ~/.local/lib/python3.6/site-packages/tstk/collapse.py \
     --minreads 1 \
    ${reads} \
     ${output}

    """
}







process CUTADAPTONE {
    cache 'deep'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id) , path("${sample_id}_*_T.fastq.gz")

    script:
    """

    cutadapt \
        -b AGATCGGAAGAGCACACGTCTGAACTCCAGTC \
        -b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -B AGATCGGAAGAGCACACGTCTGAACTCCAGTC \
        -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -o ${sample_id}_R1_T.fastq.gz  \
        -p ${sample_id}_R2_T.fastq.gz  \
        -j 10 \
        --minimum-length 10 \
        ${reads[0]} \
        ${reads[1]}
    """
}



process CUTADAPTWO {
    cache 'deep'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id) , path("${sample_id}.clean.fasta")

    script:
    output = sample_id + ".clean.fasta"
    """

    cutadapt \
        -m 20 \
        -u -6 \
        -o ${output} \
        -j 10 \
        --minimum-length 10 \
        ${reads}
    """
}



process SAMTOOLSVIEW{
    cache 'deep'

    publishDir params.outdir, mode : "copy"


    input:
    path reads
    path log
    tuple val(sample_id) , path(collapsed)

    output:
    tuple val(sample_id) , path("${sample_id}.sam")

    script:
    output = sample_id + ".sam"


    """

    samtools view ${reads} > ${output}


    """
}





process PROCESSSAM{
    cache 'deep'

    publishDir params.outdir, mode : "copy"


    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id) , path("${sample_id}_gapped.txt")

    script:
    output = sample_id + "_gapped.txt"
    """

    grep  'SA:' $reads > Chimera.sam

    grep -v 'SA:' $reads > nonChimera.sam



    Rscript $baseDir/bin/sam2hyb.R \
        nonChimera.sam \
        nochim.txt \
        0

    Rscript $baseDir/bin/samChimera2hyb.R \
            Chimera.sam \
            chim.txt \
            0

    cat nochim.txt  chim.txt > ${output}
    """
}


process STARINDEX {
    cache 'deep'

    publishDir params.outdir, mode : "copy"

    input:
    path transcriptome

    output:
    path "index"

    script:
    """
    mkdir index
    STAR --runMode genomeGenerate --genomeDir index --genomeFastaFiles $transcriptome --runThreadN 20 --limitGenomeGenerateRAM 190936728608
    """
}


process FASTQC {
    cache 'deep'

    tag "FASTQC on $sample_id"
    publishDir params.outdir, mode : "copy"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "${sample_id}_fastqc_out"

    script:
    """
    mkdir ${sample_id}_fastqc_out
    fastqc $reads -o ${sample_id}_fastqc_out -t 1
    """
}

process FASTQC2 {
    cache 'deep'
    tag "FASTQC on $sample_id"
    publishDir params.outdir, mode : "copy"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "${sample_id}_T_fastqc_out"

    script:
    """
    mkdir ${sample_id}_T_fastqc_out
    fastqc $reads -o ${sample_id}_T_fastqc_out -t 1
    """
}


process FASTQC3 {
    cache 'deep'
    tag "FASTQC on $sample_id"
    publishDir params.outdir, mode : "copy"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "${sample_id}_assembled_fastqc_out"

    script:
    """
    mkdir ${sample_id}_assembled_fastqc_out
    fastqc $reads -o ${sample_id}_assembled_fastqc_out -t 1
    """
}

process multiqc {
    cache 'deep'
    publishDir params.outdir, mode : "copy"

    input:
    path outdir

    output:
    file "multiqc_report.html"
    file "multiqc_data"

    script:
    """
    multiqc ${outdir}
    """
}

process ALIGN {
    cache 'deep'

    publishDir params.outdir, mode : "copy"

    input:
    path starIndexCh
    tuple val(pair_id), path(reads)

    output:
    path "${pair_id}Aligned.sortedByCoord.out.bam"
    path "${pair_id}Log.final.out"


    script:
    """
    mkdir $pair_id
    STAR \
        --runMode alignReads \
        --genomeDir $starIndexCh \
        --readFilesIn  ${reads}  \
        --outFileNamePrefix $pair_id  \
         --runThreadN 5 --genomeLoad NoSharedMemory \
         --outReadsUnmapped Fastx  \
         --outFilterMultimapNmax 100000 \
         --outSAMattributes All \
         --outSAMtype BAM SortedByCoordinate \
         --alignIntronMin 4 \
         --alignIntronMax 3900000000 \
         --outFilterScoreMinOverLread 0 \
         --scoreGap 0 \
         --scoreGapNoncan 0 \
         --scoreGapGCAG 0 \
         --scoreGapATAC 0 \
        --scoreGenomicLengthLog2scale -1 \
         --chimFilter None \
         --chimOutType WithinBAM HardClip --chimSegmentMin 15 \
         --chimJunctionOverhangMin 15 --chimScoreJunctionNonGTAG 0 \
         --chimScoreDropMax 80 \
         --chimMainSegmentMultNmax  100000
    """
}
