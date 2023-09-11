# 3' tagRNAseq
Analysis of _Hepatocystis_ parasites transcripts expressed by the infected _E. labiatus_ bat hosts.

First is the analysis of parasites infection intensities estimates by eg organ, infection status.
Second is the analysis of the host DETs based on categories. age, sex, season and scaled rpmh values. 

**FastQC** was used to QC the 3' Tag RNAseq raw reads from Illumina sequencing. The SE reads were short averaging between 50~70 bases. Therefore, no trimming was done so as not to have poor mapping downstreaam.

Mapping of the 3' TagSeq reads was done using **STAR** alligner in **SE** **mode**. The mapping index was generated from the merged  reference genomes for _Rousettus aegyptiacus_ and _Hepatocystis_ (from _Piliocolobus_ - aunin et, al. 2020). Script **1_STAR_mapping_AegyptiacusHepatocystis.sh** in the scripts folder describes this.

The mapped reads were indexed using **samtools index** as shown in script **2_3TagSeq_mapped_samtoolsIndex.sh** in the scripts folder.

**Features counts** in **R** was used to generate the transcripts counts for the mapped reads. Script **3_featurecounts.R** in scripts folder describes this. the resulting **countTable.RDS** counts file was saved in the intermediate file folder.

The metadata file, **tagseqRNA_metadata2023.csv**, in the input folder was loaded and modified accordingly. Script **4_metadata.R** describes the changes.

For the first analysis of _Hepatocystis_ transcriptome intensities in categories organ and infection status by PCR and microscopy of blood smears; script **5_analysis_of_hepatocystis_transcriptome_intensity.R** shows the steps done. The resulting plots are in the folder **plots/Hepatocystis_parasites_infection_intensity_plots/**. 

Annotation of the mapped transcripts was done as per script **6_hepatocystis_and_rousettus_transcripts_annotation.R**. The annotation files were the **hepatocystis_proteins.csv** and the **rousettus_proteins.csv** under the input folder. 

For the second analysis of of the host DETs based on categories. age, sex, season and scaled rpmh values; script **7_DETs_analysis_for_host_liver_transcriptome_data.R** shows the liver transcripts DETs analsis and script **8_DETs_analysis_for_host_spleen_transcriptome_data.R** shows the spleen transcripts DETs analsis. The resulting plots are in folder **plots/Host_DETs_plots/ Host_liver_DETs_plots** and **plots/Host_DETs_plots/Host_spleen_DETs_plots**. 
