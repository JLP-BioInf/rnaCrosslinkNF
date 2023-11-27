# comradesNF - Pre-processing for COMRADES experimental data


## About

This pre-processing pipeline is designed to take raw reads from a COMRADES experiment through to processed files that can be used as input for the comradesOO R package: [github](https://github.com/JLP-BioInf/comradesOO), [CRAN](https://CRAN.R-project.org/package=comradesOO)

The [comradesOO](https://github.com/JLP-BioInf/comradesOO) pipeline will work with other types of RNA crosslinking data and library reparation protocols differ between methods. The pipeline is straightforward and the only requirement for the R package is the production of a text output file with the details of the gapped alignments to the transcriptome and their location.

## The Output

\
The main output files are the files entitled *X_gapped.txt*. These are the input files for [comradesOO](https://github.com/JLP-BioInf/comradesOO). The columns of the output files are as follows:

1. Read Name
2. Read Sequence
3. Side 1 transcript ID
4. Side 1 Position start in read sequence
5. Side 1 Position end in read sequence
6. Side 1 Coordinate start in transcript
7. Side 1 Coordinate end in transcript
8. NA
9. Side 2 transcript ID
10. Side 2 Position start in read sequence
11. Side 2 Position end in read sequence
12. Side 2 Coordinate start in transcript
13. Side 2 Coordinate end in transcript
14. NA


## The Pipeline

### cutAdapt (Adaptors)

### PEAR

### collapse UMI

### Cutadape (UMI)

### STAR alignment

### Process SAM 
