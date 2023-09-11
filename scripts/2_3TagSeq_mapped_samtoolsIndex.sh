#!usr/bin/bash 

#make list of BAM files and save in a list
ls starMapped/*.STAR/*.bam > BAM.list ;
paste BAM.list | while read BAM ;

do
	samtools index "${BAM}" > "${BAM}".bai

 done
 echo "done!" 
