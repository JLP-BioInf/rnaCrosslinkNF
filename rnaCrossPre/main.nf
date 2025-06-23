/*
 * Description: This pipeline processes raw RNA-Seq data and performs quality control,
 *              adapter trimming, read alignment, and quantification
 *        Mature modules kept here
 * Steps:
 *   1. Quality control using FastQC
 *   2. Adapter trimming using Cutadapt
 *   3. Read alignment using STAR
 *   4. Quantification using featureCounts
 *   5. MultiQC report generation
 *
 *   Author - Jonathan price
 */





process gunzip {
    cache 'lenient'

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

process pear {
    cache 'lenient'
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id) , path("${sample_id}*assembled.fastq")

    script:
    """

    pear  \
        --min-assembly-length 25 \
        -f ${reads[0]}  \
        -r ${reads[1]}   \
        -o $sample_id \
        -v 7 \
        -j 5
    """
}


process fastqc {
    cache 'lenient'
    tag "$sample_id"

    publishDir(
    path: "${params.outdir}/fastQC",
    mode: 'copy'
)

    input:
    tuple val(sample_id), path(reads)
    val suffix

    output:
    path "${sample_id}_fastqc_out_${suffix}"

    script:
    """
    mkdir ${sample_id}_fastqc_out_${suffix}
    fastqc $reads -o ${sample_id}_fastqc_out_${suffix}
    """
}


process samtoolsView{
    cache 'lenient'

    publishDir(
        path: "${params.outdir}/star",
        mode: 'copy'
    )


    input:
    path reads
    path log

    output:
    tuple val(sample_id) , path("${reads.baseName}.sam")
    tuple val(sample_id) , path("${reads.baseName}.tNamesCigars")

    script:


    """
    samtools view ${reads} > ${reads.baseName}.sam

    cut -f3,6 ${reads.baseName}.sam > ${reads.baseName}.tNamesCigars

    """
}



process align {
    cache 'lenient'
    publishDir(
        path: "${params.outdir}/star",
        mode: 'copy'
    )

    input:
    path starIndexCh
    tuple val(pair_id), path(reads)
    val insertParam

    output:
    path "${pair_id}_i${insertParam}Aligned.out.bam"
    path "${pair_id}_i${insertParam}Log.final.out"


    script:
    """
    mkdir ${pair_id}_i${insertParam}



    STAR \
         --runMode alignReads \
         --genomeDir $starIndexCh \
         --readFilesIn ${reads}   \
         --outFileNamePrefix  ${pair_id}_i${insertParam}  \
         --runThreadN 12 --genomeLoad NoSharedMemory \
         --outReadsUnmapped Fastx  \
         --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 \
	       --outFilterMismatchNmax 1 \
         --alignSJoverhangMin 11 \
         --outFilterMultimapNmax 50 \
         --outSAMattributes All \
         --outSAMtype BAM Unsorted \
         --alignIntronMin 4 \
         --alignIntronMax 3900000000 \
         --scoreGap 0 \
         --scoreGapNoncan 0 \
         --scoreGapGCAG 0 \
         --scoreGapATAC 0 \
	       --limitBAMsortRAM 54317759228 \
         --scoreGenomicLengthLog2scale -1 \
         --chimFilter None \
         --chimOutType WithinBAM HardClip --chimSegmentMin 12 \
         --chimJunctionOverhangMin 12 --chimScoreJunctionNonGTAG 0 \
         --chimScoreDropMax 80 --chimNonchimScoreDropMin 20 \
         --chimMainSegmentMultNmax  10


    """
}






process cutadaptr4{
    cache 'lenient'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id) , path("${sample_id}_*_T.fastq.gz")

    script:
    """
    cutadapt --interleaved -a NNNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
      -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -O 1 ${reads[0]} <(cutadapt -l -144 ${reads[1]}) > t1.fastq

    cutadapt --interleaved -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        -b GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -b TCGGAGAATCTCGTATGCCGTCTTCTGCTTG -b ATCTCGTATGCCGTCTTCTGCTTG \
        -B GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC -B CAAGCAGAAGACGGCATACTAGATTCTCCGGA -B CAAGCAGAAGACGGCATACGAGAT \
        -O 24 --discard-trimmed  t1.fastq > t2.fastq

      cutadapt --interleaved -a AGATCGGAAGAGC -A AGATCGGAAGAGC -O 8 t2.fastq >t3.fastq

      cutadapt --interleaved --nextseq-trim=30 t3.fastq > t4.fastq

      cutadapt --interleaved -m 20:20 -o ${sample_id}_R1_T.fastq.gz -p ${sample_id}_R2_T.fastq.gz t4.fastq


    """
}

process cutadapt{
    cache 'lenient'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id) , path("${sample_id}_*_T.fastq.gz")

    script:
    """
    cutadapt \
        -q 10 \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -n 2 \
        -o ${sample_id}_t2_R1_T.fastq.gz  \
        -p ${sample_id}_t2_R2_T.fastq.gz  \
        -j 10 \
        --minimum-length 25 \
        ${reads[0]} \
        ${reads[1]}


    """
}



process starIndex {
    cache 'lenient'


    input:
    path transcriptome

    output:
    path "index"

    script:
    """
    mkdir index
    STAR  --runMode genomeGenerate \
          --genomeDir index \
          --genomeSAindexNbases 7 \
          --genomeFastaFiles $transcriptome \
          --runThreadN 20 \
          --limitGenomeGenerateRAM 190936728608
    """
}







process samtools_idxstats {

    // Publish the outputs to a specified directory
    publishDir(
        path: "${params.outdir}/sorted_bam_and_idxstats",
        mode: 'copy'
    )

    input:
    path bam_file // Input BAM file
    path log

    output:
    path "${bam_file.baseName}_sorted.bam" // Sorted BAM file
    path "${bam_file.baseName}_sorted.bam.bai" // Index file
    path "${bam_file.baseName}_idxstats.txt" // idxstats file

    script:
    """
    # Sort the BAM file
    samtools sort  ${bam_file}  ${bam_file.baseName}_sorted

    # Index the sorted BAM file
    samtools index ${bam_file.baseName}_sorted.bam

    # Generate idxstats
    samtools idxstats ${bam_file.baseName}_sorted.bam > ${bam_file.baseName}_idxstats.txt
    """
}


process processSAM{
    cache 'lenient'
    errorStrategy 'ignore'

    publishDir(
        path: "${params.outdir}/output_gapped_and_cigars",
        mode: 'copy'
    )



    input:
    tuple val(sample_id), path(reads)
    tuple val(sample_id2) , path(cigar)

    output:
    tuple val(sample_id) , path("${reads}_gapped.txt")


    script:
    output = reads + "_gapped.txt"
    """

    grep  'SA:' $reads > Chimera.sam

    grep -v 'SA:' $reads > nonChimera.sam



    sam2hyb.R \
        nonChimera.sam \
        nochim.txt \
        0

    samChimera2hyb.R \
            Chimera.sam \
            chim.txt \
            0

    cat nochim.txt  chim.txt > ${output}



    """
}


process tstk {
    cache 'lenient'

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




process multiqc {
    cache 'lenient'
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
