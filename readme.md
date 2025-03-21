# rnaCrosslinkNF - Pre-processing for [rnaCrosslinkOO](https://github.com/JLP-BioInf/comradesOO) - COMRADES experimental data


## About
\
The COMRADES (crosslinking of Matched RNA and Deep Sequencing) experimental protocol for the prediction of RNA structure in vivo was first published in 2018 (Ziv et al., 2019). The protocol has subsequently been useD to analyse the structure of SARS-CoV-2 (Ziv et al., 2020) :

* COMRADES determines in vivo RNA structures and interactions. (2018). Omer Ziv, Marta Gabryelska, Aaron Lun, Luca Gebert. Jessica Sheu-Gruttadauria and Luke Meredith, Zhong-Yu Liu,  Chun Kit Kwok, Cheng-Feng Qin, Ian MacRae, Ian Goodfellow , John Marioni, Grzegorz Kudla, Eric Miska.  Nature Methods. Volume 15. https://doi.org/10.1038/s41592-018-0121-0   

* The Short- and Long-Range RNA-RNA Interactome of SARS-CoV-2. (2020). Omer Ziv, Jonathan Price, Lyudmila Shalamova, Tsveta Kamenova, Ian Goodfellow, Friedemann Weber, Eric A. Miska. Molecular Cell,
Volume 80. https://doi.org/10.1016/j.molcel.2020.11.004


![](https://github.com/JLP-BioInf/comradesOO/blob/main/vignettes/comradesProtocol.jpg)

Ziv et al., 2020. "Virus-inoculated cells are crosslinked using clickable psoralen. Viral RNA is pulled down from the cell lysate using an array of biotinylated DNA probes, following digestion of the DNA probes and fragmentation of the RNA. Biotin is attached to crosslinked RNA duplexes via click chemistry, enabling pulling down crosslinked RNA using streptavidin beads. Half of the RNA duplexes are proximity-ligated, following reversal of the crosslinking to enable sequencing. The other half serves as a control, in which crosslink reversal proceeds the proximity ligation. "

This pre-processing pipeline is designed to take raw reads from a COMRADES experiment through to processed files that can be used as input for the comradesOO R package: [github](https://github.com/JLP-BioInf/rnaCrosslinkOO), [CRAN](https://CRAN.R-project.org/package=rnaCrosslinkOO)

The [rnaCrosslinkOO](https://github.com/JLP-BioInf/comradesOO) pipeline will work with other types of RNA crosslinking data and library reparation protocols differ between methods. The pipeline is straightforward and the only requirement for the R package is the production of a text output file with the details of the gapped alignments to the transcriptome and their location.

## The Output

\
The main output files are the files entitled *X_gapped.txt*. These are the input files for [comradesOO](https://github.com/JLP-BioInf/rnaCrosslinkOO). The columns of the output files are as follows:

1. Read Name
2. Read Sequence
3. NA
4. Side 1 transcript ID
5. Side 1 Position start in read sequence
6. Side 1 Position end in read sequence
7. Side 1 Coordinate start in transcript
8. Side 1 Coordinate end in transcript
9. NA
10. Side 2 transcript ID
11. Side 2 Position start in read sequence
12. Side 2 Position end in read sequence
13. Side 2 Coordinate start in transcript
14. Side 2 Coordinate end in transcript
15. NA


![](https://github.com/JLP-BioInf/comradesOO/blob/main/vignettes/inputFileSchematic.jpg)

## The Pipeline
\
The pipeline is straightforward, here is a run down of the steps and parameters in the pre-processing. You can use this pipeline or alter it to suit your own library preparation protocol.
