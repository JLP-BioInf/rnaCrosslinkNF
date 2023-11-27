################################################################################
# test the mapping with STAR
################################################################################
options(warn=-1)
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyverse))
################################################################################
# Mapped tt he transcriptome: 
# 
# sbatch -p 1804 --array 1 ../sbatchSTARCOMRADES.sh ~/referenceGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/STARtranscriptsIndex/ ../5_star_transcripts
# 
#
# # START OF SCRIPT  BELOW
#
#
# #!/bin/bash
# #SBATCH -n 5                         # Number of cores
# #SBATCH -N 1                         # Ensure that all cores are on one machine
# #SBATCH --mem=8G                     # Memory pool for all cores (see also --mem-per-cpu)
# #SBATCH --mail-type=END,FAIL         # Type of email notification- BEGIN,END,FAIL,ALL
# #SBATCH --mail-user=jlp76@cam.ac.uk  # Email to which notifications will be sent
# 
# #find all files in the current directory and grab them 1 by 1 depeding on which task is being run
# arrayfile=`ls ./*.assembled.tstk.fasta | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
# 
# echo $arrayfile
# 
# STAR \
# --runMode alignReads \
# --genomeDir $1 \
# --readFilesIn $arrayfile \
# --outFileNamePrefix ${2}/${arrayfile%.assembled.tstk.fastq}  \
# --runThreadN 5 --genomeLoad NoSharedMemory \
# --outReadsUnmapped Fastx  --outFilterMultimapNmax 200 \
# --outFilterScoreMinOverLread 0 \
# --outFilterMatchNminOverLread 0 \
# --outSAMattributes All \
# --outSAMtype BAM SortedByCoordinate \
# --alignIntronMin 1 --scoreGap 0 --scoreGapNoncan 0 --scoreGapGCAG 0 \
# --scoreGapATAC 0 --scoreGenomicLengthLog2scale -1 --chimFilter None \
# --chimOutType WithinBAM HardClip --chimSegmentMin 5 \
# --chimJunctionOverhangMin 5 --chimScoreJunctionNonGTAG 0 \
# --chimScoreDropMax 80 \
# --twopassMode Basic

# # END OF SCRIPT  





#can reads part of bam file with IF BAM IS TOO BIG  :
#yieldSize(bamFile) <- 1
#open(bamFile)
#scanBam(bamFile)[[1]]$seq # repeat this
#close(bamFile)
#yieldSize(bamFile) <- NA

args <- commandArgs(trailingOnly = TRUE)




#read in  the sam files:
gappedAln = read.table(args[1], header = F, skip = args[3], fill = T)
gappedAln = read.table("/Users/jp/projects/COMRADES/nextflow/nf-corePipeline/testChimeras/nonChimeraH.sam", header = F, skip = 0, fill = T)
colnames(gappedAln) = c("readName", "flag", "geneID", "start", "mapQ", "cigar", "strand", "x","y", "seq", 
                        "flag1","flag2","flag3","flag4","flag5","flag6","flag7","flag8","flag9")




# split the cigar into M - N - M
cigarSplit1 = str_split(gappedAln$cigar, "[A-Z]")
cigarSplit2 = str_split(gappedAln$cigar, "[0-9]")

# get a more sensible format
#cigarSplit2 = lapply(cigarSplit2, paste,collapse = "" )

# get a more sensible format
for(i in 1:length(cigarSplit1)){
     #cigarSplit1[[i]] = print(paste(cigarSplit1[[i]], collapse = ""))
    cigarSplit2[[i]] = paste(cigarSplit2[[i]], collapse = "")
}

# Now chek for MNM CIGARS
cigarSplit2 = unlist(cigarSplit2)
goodIndexes = grep('^S{0,1}MNMS{0,1}$', cigarSplit2)

#subset the gappen alignments
cigarSplit2 = cigarSplit2[goodIndexes]
gappedAln = gappedAln[goodIndexes,]

gappedAln = gappedAln[,1:19]

# make solumns to score the number of BPs for each segment
gappedAln$S1 = 0
gappedAln$M1 = 0
gappedAln$N1 = 0
gappedAln$M2 = 0
gappedAln$S2 = 0


# now get those that start with S
Sgapped = gappedAln[grep('^S', cigarSplit2),]
cigarSplitS = str_split(Sgapped$cigar, "[A-Z]")
Sgapped[,c(20:24)] = as.data.frame(do.call(rbind, cigarSplitS))[,1:5]


# now get the ones that dont start with S
Mgapped = gappedAln[grep('^M', cigarSplit2),]
cigarSplitS = str_split(Mgapped$cigar, "[A-Z]")
Mgapped[,c(21:24)] = as.data.frame(do.call(rbind, cigarSplitS))[,1:4]


#put back together
gappedAln = rbind.data.frame(Sgapped,Mgapped)




# fill in the blank S with 0
gappedAln = gappedAln %>%
    mutate(S1 = replace(S1, S1 == "",0))
gappedAln = gappedAln %>%
    mutate(S2 = replace(S2, S2 == "",0))




#dim(gappedAln)
gappedAln$S1 = as.numeric(gappedAln$S1)
gappedAln$M1 = as.numeric(gappedAln$M1)
gappedAln$N1 = as.numeric(gappedAln$N1)
gappedAln$M2 = as.numeric(gappedAln$M2)
gappedAln$S2 = as.numeric(gappedAln$S2)
gappedAln$start = as.numeric(gappedAln$start)








# make the hyb style output
hybStyledf = cbind.data.frame(gappedAln[,c("readName","seq")],
                              rep(".",nrow(gappedAln)), 
                              gappedAln[,"geneID"],
                              gappedAln$S1+1,
                              gappedAln$S1+gappedAln$M1,
                              gappedAln$start,
                              gappedAln$start+gappedAln$M1,
                              rep(".",nrow(gappedAln)),
                              gappedAln$geneID,
                              gappedAln$S1+gappedAln$M1,
                              gappedAln$S1+gappedAln$M1+gappedAln$M2,
                              gappedAln$start+gappedAln$M1+gappedAln$N1,
                              gappedAln$start+gappedAln$M1+gappedAln$N1+gappedAln$M2,
                              rep(".",nrow(gappedAln)))

hybStyledf = hybStyledf[(hybStyledf[,14] - hybStyledf[,13]) >10 &
                            (hybStyledf[,8] - hybStyledf[,7] ) > 10, ]

hybStyledf = hybStyledf[complete.cases(hybStyledf),]



write.table(hybStyledf,file = args[2], quote = F, row.names = F, col.names = F, sep ="\t")
