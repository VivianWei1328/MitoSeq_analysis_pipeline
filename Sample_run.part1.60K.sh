###### Part 1 for MitoSeq, dumultiplexing and recover UMI 
#!/bin/bash
#!/usr/bin/perl

### before running, replace $PATH_to_basecalled_files and $Run_name to actual value
### Need conda environments for longread_umi_Intel and samtools_env_Intel
###### This pipeline use samtools 1.2 in samtools_env_Intel
###### samtools 1.21 in regular and samtools_env



cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101
pwd

source /opt/miniconda3/etc/profile.d/conda.sh

conda activate longread_umi
conda deactivate

parallel

cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101


samtools fastq $Run_name-13V-run1.101.dup.bam > $Run_name-13V-run1.101.dup.fastq

 
wc -l $Run_name-13V-run1.101.dup.fastq 


cp $PATH_to_barcode/Vbarcores.16.top.21.fasta  .

##################################################
################ cut 1 
################ cutadapt 5.1 on Mac Studio

split -l 1000000 $Run_name-13V-run1.101.dup.fastq chunk_

for chunk in chunk_*; do
  (
    echo "Processing $chunk"
    cutadapt -g file:Vbarcores.16.top.21.fasta -O 20 -e 0.15 -m 150 -o ${chunk}.{name}.fastq "$chunk" > "$Run_name-13V-run1.101.dup.${chunk}.cut1.21.m150.log" 2>&1
    echo "Finished $chunk"
  ) &
done

time wait 

for sample in V01 V02 V03 V04 V05 V13 V07 V08 V09 V10 V11 V12 F01 F02 F03 F04 unknown; do cat chunk_*.${sample}.fastq > $Run_name-13V-run1.101.dup.${sample}.cut1.21.m150.fastq; done

 wc -l $Run_name-13V-run1.101.dup.*.cut1.21.m150.fastq

for sample in V01 V02 V03 V04 V05 V13 V07 V08 V09 V10 V11 V12 F01 F02 F03 F04 unknown; do wc -l chunk_*.${sample}.fastq; done

for f in $(ls $Run_name-13V-run1.101.dup.chunk_*.log); do echo $f; head -20 $f; done
$Run_name-13V-run1.101.dup.chunk_aa.log


########## activate the conda environment for Intel chip, if using an Apple silicon chip  

conda activate longread_umi_Intel
cd /opt/miniconda3/envs/longread_umi_Intel/longread_umi_Intel

file="$Run_name-13V-run1.101.dup.V01.cut1.21.m150.fastq $Run_name-13V-run1.101.dup.V02.cut1.21.m150.fastq $Run_name-13V-run1.101.dup.V03.cut1.21.m150.fastq $Run_name-13V-run1.101.dup.V04.cut1.21.m150.fastq $Run_name-13V-run1.101.dup.V05.cut1.21.m150.fastq $Run_name-13V-run1.101.dup.V13.cut1.21.m150.fastq $Run_name-13V-run1.101.dup.V07.cut1.21.m150.fastq $Run_name-13V-run1.101.dup.V08.cut1.21.m150.fastq $Run_name-13V-run1.101.dup.V09.cut1.21.m150.fastq $Run_name-13V-run1.101.dup.V10.cut1.21.m150.fastq $Run_name-13V-run1.101.dup.V11.cut1.21.m150.fastq $Run_name-13V-run1.101.dup.V12.cut1.21.m150.fastq $Run_name-13V-run1.101.dup.F01.cut1.21.m150.fastq"


time for f in $file; do echo $f; longread_umi nanopore_pipeline -d $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101/$f -v 30 -o $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101/${f%.fastq}.LRUMI_out -s 90 -e 90 -m 5000 -M 20000 -f CAAGCAGAAGACGGCATACGAGATGTTGCAGCGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT -F GATCACAGGTCTATCACCCTATT -r AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -R CATCGTGATGTCTTATTTAAGGG -c 3 -p 1 -q r941_min_high_g330 -t 12 -T 1; wc -l  $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101/${f%.fastq}.LRUMI_out/umi_binning/umi_ref/umi12f.fa; done


cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101 
pwd
echo "Generate the nano.input.file"

cat <<EOF >> nano.input.file
MT-28-13V-run1.101.dup.V01.cut1.21.m150.fastq.gz
MT-28-13V-run1.101.dup.V02.cut1.21.m150.fastq.gz
MT-28-13V-run1.101.dup.V03.cut1.21.m150.fastq.gz
MT-28-13V-run1.101.dup.V04.cut1.21.m150.fastq.gz
MT-28-13V-run1.101.dup.V05.cut1.21.m150.fastq.gz
MT-28-13V-run1.101.dup.V13.cut1.21.m150.fastq.gz
MT-28-13V-run1.101.dup.V07.cut1.21.m150.fastq.gz
MT-28-13V-run1.101.dup.V08.cut1.21.m150.fastq.gz
MT-28-13V-run1.101.dup.V09.cut1.21.m150.fastq.gz
MT-28-13V-run1.101.dup.V10.cut1.21.m150.fastq.gz
MT-28-13V-run1.101.dup.V11.cut1.21.m150.fastq.gz
MT-28-13V-run1.101.dup.V12.cut1.21.m150.fastq.gz
MT-28-13V-run1.101.dup.F01.cut1.21.m150.fastq.gz
EOF


for f in $file; do wc -l  ./${f%.fastq}.LRUMI_out/umi_binning/umi_ref/umi12f.fa; done

conda deactivate

###### map to Hg38 reference genome

cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101 

time for f in $(cat nano.input.file); do echo $f; mkdir ${f%.fastq.gz}; mv ${f%.fastq.gz}.* ${f%.fastq.gz}; cd ${f%.fastq.gz}; pwd; minimap2 -x map-ont -t 24 --frag=yes --secondary=yes -N 18 -p 0.8 -a -o ${f%.fastq.gz}.mini.sam $PATH_to_reference/Ref/Hg38/grch38_1kgmaj.fa ${f%.gz}; cd ..; done 


cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101


time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; cp ./${f%.fastq.gz}.LRUMI_out/umi_binning/umi_ref/umi12f.fa .; sed 's/^>\([^>]*\)>.*/>\1/' umi12f.fa > test.umi12.f.fasta; cd ..; pwd; done


conda activate samtools_env_Intel

time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; samtools view -bS ${f%.fastq.gz}.mini.sam | samtools sort - ${f%.fastq.gz}.mini.sort; samtools index ${f%.fastq.gz}.mini.sort.bam; samtools idxstats ${f%.fastq.gz}.mini.sort.bam; cd ..; done 

conda deactivate


echo "end of part 1, start part 2"
cd $PATH_to_script_files

time sh $Run_name-13V-run1.101.dup.part2.60K.sh

























