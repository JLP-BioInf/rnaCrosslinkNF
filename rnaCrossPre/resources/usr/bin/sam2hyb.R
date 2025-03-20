#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))


args <- commandArgs(trailingOnly = TRUE)

#read in  the sam files:
gappedAln = read.table(args[1], header = F, skip = args[3], fill = T)
#gappedAln = read.table("/Users/jp/projects/COMRADES/nextflow/nf-corePipeline/testChimeras/nonChimeraH.sam", header = F, skip = 0, fill = T)

colnames(gappedAln) = c("readName", "flag", "geneID", "start", "mapQ", "cigar", "strand", "x","y", "seq", 
                        "flag1","flag2","flag3","flag4","flag5","flag6","flag7","flag8","flag9")

# split the cigar into M - N - M
cigarSplit1 = str_split(gappedAln$cigar, "[A-Z]")
cigarSplit2 = str_split(gappedAln$cigar, "[0-9]")

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

writeLines(capture.output(sessionInfo()), "sessionInfo.txt") 
write.table(hybStyledf,file = args[2], quote = F, row.names = F, col.names = F, sep ="\t")
