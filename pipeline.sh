
#!/bin/bash;

echo -e "\n\n\n";
echo "######################################";
echo "#variable initialisation#";
echo "######################################";

read1="${1}";
read2="${2}";
ref_fasta="${3}";
known_variants="${4}";
adapter_seq1="${5}";
adapter_seq2="${6}"


#Root directories 
pipeline_root="/mnt/x/linux/project/";
picard_root="/mnt/x/linux/picard/";
gatk_root="/mnt/x/linux/gatk-4.2.0.0/";
snpeff_root="/mnt/x/linux/test/snpeff/";
fastqc_root="/mnt/x/linux/FastQC/";

echo -e "\n\n\n";
echo "###################################";
echo "#Pre-processing raw data#";
echo "###################################";

"${fastqc_root}"fastqc "${1}" "${2}";

cutadapt -a "${5}" -A "${6}"\
                     -o "${1}"_trimmed -p "${2}"_trimmed\
                     "${1}" "${2}";



echo -e "\n\n\n";
echo "###################################";
echo "#Indexing reference genome#";
echo "###################################";

bwa index "${3}";

samtools faidx "${3}";

java -jar "${gatk_root}"gatk.jar CreateSequenceDictionary -R "${3}";

echo Done 

echo -e "\n\n\n";
echo "###########";
echo "#Alignment#";
echo "###########";

bwa mem -M -t 4 "${3}" "${1}"_trimmed "${2}"_trimmed\
         >"${1}"_"${2}"_BWA.sam;

samtools view -bS "${1}"_"${2}"_BWA.sam -o "${1}"_"${2}"_BWA.bam;

samtools index "${1}"_"${2}"_BWA.bam;


echo Done

echo -e "\n\n\n";
echo "##############";
echo "#sorting BAM #";
echo "##############";

java -Xmx10g -jar "${picard_root}"picard.jar SortSam VALIDATION_STRINGENCY=SILENT\
                   I="${1}"_"${2}"_BWA.bam O="${1}"_"${2}"_BWA_sorted.bam SORT_ORDER=coordinate;

java -Xmx10g -jar "${picard_root}"picard.jar MarkDuplicates VALIDATION_STRINGENCY=SILENT\
             I="${1}"_"${2}"_BWA_sorted.bam O="${1}"_"${2}"_BWA_sorted_PCR.bam\
             REMOVE_DUPLICATES=true M="${1}"_"${2}"_BWA_sorted_PCR.metrics;

java -Xmx10g -jar "${picard_root}"picard.jar AddOrReplaceReadGroups VALIDATION_STRINGENCY=SILENT\
           I="${1}"_"${2}"_BWA_sorted_PCR.bam O="${1}"_"${2}"_BWA_sorted_PCR_RG.bam\
          SO=coordinate RGID=HISeq RGLB=NIST RGPL=illumina RGPU=R1 RGSM=NA12878 CREATE_INDEX=true;


echo creating index file
samtools index "${1}"_"${2}"_BWA_sorted_PCR_RG.bam

echo Done

echo -e "\n\n\n";
echo "###################################";
echo "#Base Quality score recaliberation#";
echo "###################################";

tabix "${4}";

echo 1.craeting recal data table

java -jar "${gatk_root}"gatk.jar BaseRecalibrator -R "${3}" -I "${1}"_"${2}"_BWA_sorted_PCR_RG.bam\
              --known-sites "${4}" -O recal_data.table;

echo 2.Applying BQSR 

java -jar "${gatk_root}"gatk.jar ApplyBQSR -R "${3}" -I "${1}"_"${2}"_BWA_sorted_PCR_RG.bam\
                      -bqsr recal_data.table -O "${1}"_"${2}"_BWA_sorted_PCR_RG_BQSR.bam;



echo -e "\n\n\n";
echo "##################";
echo "#Variant Calling #";
echo "##################";
 
echo variant calling by bcftools 

bcftools mpileup -Ou -f "${3}" "${1}"_"${2}"_BWA_sorted_PCR_RG.bam\
            |bcftools call -mv -Ov -o "${1}"_"${2}"_BWA_sorted_PCR_RG_BQSR_bcftools_call.vcf;
echo Done

echo -e "\n\n\n";
echo Variant calling by Freebayes 

echo 1.Deleting non DNA character using seqtk

seqtk randbase "${3}" > "${3}"_freebayes.fa;

echo 2.creating .fai file

samtools faidx "${3}"_freebayes.fa;

echo 3.variant calling

freebayes -f "${3}"_freebayes.fa "${1}"_"${2}"_BWA_sorted_PCR_RG.bam\
                >"${1}"_"${2}"_BWA_PCR_RG_BQSR_freebayes.vcf;


echo -e "\n\n\n";
echo Variant calling by GATK Haplotype caller
 

java -jar "${gatk_root}"gatk.jar HaplotypeCaller -R "${3}" -I "${1}"_"${2}"_BWA_sorted_PCR_RG.bam\
           -O "${1}"_"${2}"_BWA_PCR_RG_BQSR_gatk.vcf;

echo Done

echo -e "\n\n\n";
echo "###########################################";
echo "#Merging all VCF files#";
echo "###########################################";

bgzip "${1}"_"${2}"_BWA_sorted_PCR_RG_BQSR_bcftools_call.vcf;
tabix "${1}"_"${2}"_BWA_sorted_PCR_RG_BQSR_bcftools_call.vcf.gz;

bgzip "${1}"_"${2}"_BWA_PCR_RG_BQSR_freebayes.vcf;
tabix "${1}"_"${2}"_BWA_PCR_RG_BQSR_freebayes.vcf.gz;

bgzip "${1}"_"${2}"_BWA_PCR_RG_BQSR_gatk.vcf;
tabix "${1}"_"${2}"_BWA_PCR_RG_BQSR_gatk.vcf.gz;

bcftools merge --force-sample "${1}"_"${2}"_BWA_sorted_PCR_RG_BQSR_bcftools_call.vcf.gz "${1}"_"${2}"_BWA_PCR_RG_BQSR_freebayes.vcf.gz "${1}"_"${2}"_BWA_PCR_RG_BQSR_gatk.vcf.gz>merged_all.vcf;

bcftools norm -d both merged_all.vcf >  "${1}"_"${2}"_merged_all.vcf;

bcftools annotate -x 'FORMAT/GL' "${1}"_"${2}"_merged_all.vcf>"${1}"_"${2}"_merged_all_GL.vcf;


echo done

echo -e "\n\n\n";
echo "######################";
echo "#GATK Hard Filtering##";
echo "######################";

java -jar "${gatk_root}"gatk.jar VariantFiltration -V "${1}"_"${2}"_merged_all_GL.vcf\
			--filter-expression "QUAL < 20" --filter-name LowQUAL --filter-expression "DP < 5" --filter-name LowCoverage\
			-O "${1}"_"${2}"_merged_all_GL_filter.vcf;
			
bcftools view -f PASS "${1}"_"${2}"_merged_all_GL_filter.vcf > "${1}"_"${2}"_merged_all_GL_filtered.vcf;



echo -e "\n\n\n";
echo "#########################";
echo "#Annotation using snpEff#";
echo "#########################";

echo ID column annotate by SnpSift

java -jar "${snpeff_root}"SnpSift.jar annotate -dbsnp "${1}"_"${2}"_merged_all_GL_filtered.vcf\
                       > "${1}"_"${2}"_merged_all_GL_filtered_SnpSift.vcf;



echo Annotation by snpEff

java -Xmx8G -jar "${snpeff_root}"snpEff.jar -c snpEff.config GRCh37.75 "${1}"_"${2}"_merged_all_GL_filtered_SnpSift.vcf\
 > "${1}"_"${2}"_merged_all_GL_filtered_SnpSift_SnpEff.vcf




echo done 

exit 0;

