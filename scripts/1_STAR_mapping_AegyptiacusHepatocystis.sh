#!usr/bin/bash
# For every name in the file
find /SAN/RNASeqHepatoCy/HostTranscriptome/raw_reads/ -name "*.fastq.gz" | while read SAMPLE; do  

# Get single file name
FILEBASE=$(basename "${SAMPLE%.fastq.gz}") 
echo "Running STAR on" ${SAMPLE%.fastq.gz}  

# Make new directory for every sample
mkdir /SAN/RNASeqHepatoCy/HostTranscriptome/host_merged_Raegypticus_and_HepatocystisAunin/STAR/starMapped/$FILEBASE.STAR

# Enter the new directory
cd /SAN/RNASeqHepatoCy/HostTranscriptome/host_merged_Raegypticus_and_HepatocystisAunin/STAR/starMapped/$FILEBASE.STAR

# Align with STAR  	
/tools/STAR-2.5.4b/bin/Linux_x86_64/STAR --runThreadN 20 --genomeDir /SAN/RNASeqHepatoCy/HostTranscriptome/host_merged_Raegypticus_and_HepatocystisAunin/STAR/starIndex/ --readFilesIn $SAMPLE --readFilesCommand gunzip -c --outFileNamePrefix /SAN/RNASeqHepatoCy/HostTranscriptome/host_merged_Raegypticus_and_HepatocystisAunin/STAR/starMapped/$FILEBASE.STAR/$FILEBASE.STAR_ --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard  

done 

echo "done!"
