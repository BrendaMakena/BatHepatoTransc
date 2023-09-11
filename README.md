# 3' tagRNAseq
Analysis of _Hepatocystis_ parasites transcripts expressed by the infected _E. labiatus_ bat hosts.

First is the analysis of parasites infection intensities estimates by eg organ, infection status.
Second is the analysis of the host DETs based on categories. age, sex, season and scaled rpmh values. 

FastQC was used to QC the 3' Tag RNAseq raw reads from Illumina sequencing. The SE reads were short averaging between 50~70 bases. Therefore, no trimming was done so as not to have poor mapping downstreaam.

Mapping of the 3' TagSeq reads was done using STAR alligner in SE mode.  
