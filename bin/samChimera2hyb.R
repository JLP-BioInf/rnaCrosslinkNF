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




args <- commandArgs(trailingOnly = TRUE)


gappedAln = read.table(args[1], header = F, skip = args[3], fill = T)
gappedAln = gappedAln[ nchar(as.character(gappedAln$V21)) > 55 , ]
#read in  the sam files:

colnames(gappedAln) = c("readName", "flag", "geneID", "start", "mapQ", "cigar", "strand", "x","y", "seq", 
                        "flag1","flag2","flag3","flag4","flag5","flag6","flag7","flag8","flag9", "flasg10", "chimera")



# Now we need to group the two sets of alignments




# split the cigar into M - N - M
cigarSplit1 = str_split(gappedAln$cigar, "[A-Z]")
cigarSplit2 = str_split(gappedAln$cigar, "[0-9]")



# get a more sensible format
for(i in 1:length(cigarSplit1)){
    # cigarSplit1[[i]] = print(paste(cigarSplit1[[i]], collapse = ""))
    cigarSplit2[[i]] = paste(cigarSplit2[[i]], collapse = "")
}


# Now chek for MNM CIGARS
cigarSplit2 = unlist(cigarSplit2)
goodIndexes = cigarSplit2 %in% c("HM", "MH", "SMH", "HMS", "SM","MS")

#subset the gappen alignments
cigarSplit2 = cigarSplit2[goodIndexes]
gappedAln = gappedAln[goodIndexes,]

gappedAln = gappedAln[,1:21]


gappedAln = gappedAln[order(gappedAln$readName),]

# get thos parts of the table that have the same readID next door
repIndexes = c()
for(i in seq(from = 1,to = nrow(gappedAln)-1)){
    if(gappedAln$readName[i] == gappedAln$readName[i+1]){
        repIndexes = c(repIndexes,i)
    }else{
        #print("no")
    }
}




# now subset the gappedAln to get the lines with two read names next to eachother
repIndexes2 = repIndexes +1
repIndexes = c(repIndexes,repIndexes2)
repIndexes = repIndexes[order(repIndexes)]
gappedAln = gappedAln[repIndexes,]


# make the cigar split again
# split the cigar into M - N - M
cigarSplit1 = str_split(gappedAln$cigar, "[A-Z]")
cigarSplit2 = str_split(gappedAln$cigar, "[0-9]")



# get a more sensible format
for(i in 1:length(cigarSplit1)){
    # cigarSplit1[[i]] = print(paste(cigarSplit1[[i]], collapse = ""))
    cigarSplit2[[i]] = paste(cigarSplit2[[i]], collapse = "")
}






# make solumns to score the number of BPs for each segment
gappedAln$H1 = 0
gappedAln$S1 = 0
gappedAln$M1 = 0
gappedAln$S2 = 0
gappedAln$H2 = 0




# now get those that start with MH
Agapped = gappedAln[ cigarSplit2 == "MH",]
cigarSplitS = str_split(Agapped$cigar, "[A-Z]")
Agapped[,c(24)] = as.data.frame(do.call(rbind, cigarSplitS))[,1]
Agapped[,c(26)] = as.data.frame(do.call(rbind, cigarSplitS))[,2]


# now get those that start with HM
Bgapped = gappedAln[ cigarSplit2 == "HM",]
cigarSplitS = str_split(Bgapped$cigar, "[A-Z]")
Bgapped[,c(22)] = as.data.frame(do.call(rbind, cigarSplitS))[,1]
Bgapped[,c(24)] = as.data.frame(do.call(rbind, cigarSplitS))[,2]




# now get those that start with MS
Cgapped = gappedAln[ cigarSplit2 == "MS",]
cigarSplitS = str_split(Cgapped$cigar, "[A-Z]")
Cgapped[,c(24:25)] = as.data.frame(do.call(rbind, cigarSplitS))[,1:2]

# now get those that start with SM
Dgapped = gappedAln[ cigarSplit2 == "SM",]
cigarSplitS = str_split(Dgapped$cigar, "[A-Z]")
Dgapped[,c(23:24)] = as.data.frame(do.call(rbind, cigarSplitS))[,1:2]

# now get those that start with SMH
Egapped = gappedAln[ cigarSplit2 == "SMH",]
cigarSplitS = str_split(Egapped$cigar, "[A-Z]")
Egapped[,c(23:24)] = as.data.frame(do.call(rbind, cigarSplitS))[,1:2]
Egapped[,c(26)] = as.data.frame(do.call(rbind, cigarSplitS))[,3]


# now get the ones that dont start with HMS
Fgapped = gappedAln[ cigarSplit2 == "HMS",]
cigarSplitS = str_split(Fgapped$cigar, "[A-Z]")
Fgapped[,c(23:24)] = as.data.frame(do.call(rbind, cigarSplitS))[,2:3]
Fgapped[,c(22)] = as.data.frame(do.call(rbind, cigarSplitS))[,1]


#put back together
gappedAln = rbind.data.frame(Agapped,Bgapped,
                             Cgapped,Dgapped,
                             Egapped,Fgapped)


gappedAln = gappedAln[order(gappedAln$readName,gappedAln$flag ),]





# fill in the blank S with 0
#gappedAln = gappedAln %>%
#    mutate(S1 = replace(S1, S1 == "",0))
#gappedAln = gappedAln %>%
#    mutate(S2 = replace(S2, S2 == "",0))

# nwo we need to take eachof the 


#dim(gappedAln)
gappedAln$S1 = as.numeric(gappedAln$S1)
gappedAln$H1 = as.numeric(gappedAln$H1)
gappedAln$M1 = as.numeric(gappedAln$M1)
gappedAln$H2 = as.numeric(gappedAln$H2)
gappedAln$S2 = as.numeric(gappedAln$S2)
gappedAln$start = as.numeric(gappedAln$start)






# make the hyb style output
hybStyledf = cbind.data.frame(gappedAln[,c("readName","seq")],
                              rep(".",nrow(gappedAln)), 
                              gappedAln[,"geneID"],
                              gappedAln$S1+gappedAln$H1+1,
                              gappedAln$S1+gappedAln$H1+gappedAln$M1,
                              gappedAln$start,
                              gappedAln$start+gappedAln$S1+gappedAln$M1-1,
                              rep(".",nrow(gappedAln)),
                              gappedAln$H1, gappedAln$S1, gappedAln$M1, gappedAln$S2, gappedAln$H2)



# get thos parts of the table that have the same readID next door

# first find out how big the DF will be 

indexes = c()
for(i in seq(from = 1,to = nrow(hybStyledf)-1)){
    if(hybStyledf[i,1] == hybStyledf[i+1,1]){
        if(hybStyledf[i,4] == hybStyledf[i+1,4]){
            if( (hybStyledf[i,7] < hybStyledf[i+1,7] & hybStyledf[i,8] > hybStyledf[i+1,7]) | 
                (hybStyledf[i,7] > hybStyledf[i+1,7] & hybStyledf[i,7] < hybStyledf[i+1,8]) ){
                 indexes = c(indexes,i,i+1)
        }
    }
    }else{
        #print("no")
    }
}



homDF = hybStyledf[indexes,]
hetDF = hybStyledf[-indexes,]





# now add the two lines together to make one table


# get the indexes for evens and the iundexes ofr odds
x = seq(from = 1,to = nrow(hetDF)-1, by = 2)
y = seq(from = 2,to = nrow(hetDF), by = 2)
hetDF2 = cbind.data.frame(hetDF[x,],hetDF[y,])

# now make it into a hybStyle Df
hetDF3 = hetDF2[,c(1,2,3,4,5,6,7,8,9,18,19,20,21,22,23)]


write.table(hetDF3,file = args[2], quote = F, row.names = F, col.names = F, sep ="\t")

# 
# gappedAln[gappedAln$readName =="772625-16",]
# hist(gappedAln$start+gappedAln$M1+gappedAln$N1+gappedAln$M2-1)
# max(gappedAln$start+gappedAln$M1+gappedAln$N1+gappedAln$M2)
# 
# 
# hist(cds@hybFiles$`ENSG000000XXXXX_NR003286-2_RN18S1_rRNA`$original$s1$V14)
# 
# 
# # find those on the streak
# cds@hybFiles$`ENSG000000XXXXX_NR003286-2_RN18S1_rRNA`$original$s1[which(cds@hybFiles$`ENSG000000XXXXX_NR003286-2_RN18S1_rRNA`$original$s1$V14 - cds@hybFiles$`ENSG000000XXXXX_NR003286-2_RN18S1_rRNA`$original$s1$V8 > 1000),]
# 
# 
# gappedAln[gappedAln$readName =="936403-15",]
# cds@hybFiles$`ENSG000000XXXXX_NR003286-2_RN18S1_rRNA`$original$s1[cds@hybFiles$`ENSG000000XXXXX_NR003286-2_RN18S1_rRNA`$original$s1$V1 == "936403-15",]
# 
# hybStyledf[hybStyledf$readName == "936403-15",]
# hist(hybStyledf$`gappedAln$start + gappedAln$M1 + gappedAln$N1 + gappedAln$M2 - ` - hybStyledf$`gappedAln$start + gappedAln$M1 - 1`)
