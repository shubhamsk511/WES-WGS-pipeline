# WES-WGS-pipeline
we developed a pipeline to process WGS/WES data for variant discovery. The pipeline contain processes from QC to annotatating variants. We used 3 variant callers
1.BCFtools 
2.Freebayes
3.GATK Haplotype caller 
After merging the Variants from 3 variant callers they will undergo through GATK filteration(QUAL>20 HIGH QUAL ; DP> 5 HIGH COVERAGE) 
The Annotations will given by using SnpSift dbsnp and SnpEff. 
