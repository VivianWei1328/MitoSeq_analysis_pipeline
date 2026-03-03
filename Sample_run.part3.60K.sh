###### This pipeline use samtools 1.2 in samtools_env_Intel
###### samtools 1.21 in regular and samtools_env
###### Try to use samtools v1.21 consensus
#!/bin/bash
#!/usr/bin/perl

### before running, replace $PATHs to actual value
### Need conda environments for longread_umi_Intel and samtools_env_Intel

cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101 
pwd
source /opt/miniconda3/etc/profile.d/conda.sh

#conda create -n samtools_env

#conda install -c bioconda samtools=1.21 bcftools htslib
samtools --version

#samtools consensus -a in.bam -o ref.fa -config r10.4_sup

#samtools consensus -a -X r10.4_sup test.umi12.f.u.c.94.4up.clusters1.bam -o test.umi12.f.u.c.94.4up.clusters1.con.fasta
   
#samtools consensus -a -f fastq  -X r10.4_sup test.umi12.f.u.c.94.4up.clusters1.bam -o test.umi12.f.u.c.94.4up.clusters1.con.fastq



cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101 

### consensus -d 4, -c 0.75 -X r10.4_sup

### consensus for 4up cluster first 

time for f in $(cat nano.input.file)
do 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd	
	echo $f; 
	echo "Generate consensus for each cluster"
	

	rm ../test.umi12.f.60K.u.c.94.4up.clusters.con.fix.fastq
	touch ../test.umi12.f.60K.u.c.94.4up.clusters.con.fix.fastq
#### used samtools v1.21 consensus	
	for f in $(cat ../test.umi12.f.60K.u.c.94.4up.clusters.rename); do samtools consensus -r chrM -d 4 -f fastq  -X r10.4_sup $f.bam -o $f.60K.con.fastq; sed "s/^\(@[^ ]*\)/\1|${f}/g" $f.60K.con.fastq > $f.60K.con.rename.fastq; /Users/sw3203/Documents/Research/Script_lib/fix_fastq.sh $f.60K.con.rename.fastq > $f.60K.con.fix.fastq; rm $f.60K.con.rename.fastq; cat $f.60K.con.fix.fastq >> ../test.umi12.f.60K.u.c.94.4up.clusters.con.fix.fastq; done 	
	
	wc -l ../test.umi12.f.60K.u.c.94.4up.clusters.con.fix.fastq

	cd ../..
	pwd
	
done


### keep the consensus on chrM 
time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	
	pwd	
	
	awk 'BEGIN {keep = 0} /^@/ {keep = index($0, "chrM") > 0} {if (keep) print}' test.umi12.f.60K.u.c.94.4up.clusters.con.fix.fastq > test.umi12.f.60K.u.c.94.4up.clusters.con.fix.chrM.fastq
	
	wc -l test.umi12.f.60K.u.c.94.4up.clusters.con.fix.fastq
	wc -l test.umi12.f.60K.u.c.94.4up.clusters.con.fix.chrM.fastq
	
	cd ..
	pwd
	
done

### map the consensus to hg38, cluster with USEARCH

##### samtools 1.2 in samtools_env_Intel
conda activate samtools_env_Intel
cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101 

time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; minimap2 -x map-ont -t 12 --frag=yes --secondary=yes -N 24 -p 0.8 -a -o ${f%.fastq.gz}.60K.4up.con.chrM.sam /Users/sw3203/Documents/Fun/Ref/Hg38/grch38_1kgmaj.fa test.umi12.f.60K.u.c.94.4up.clusters.con.fix.chrM.fastq; samtools view -bS ${f%.fastq.gz}.60K.4up.con.chrM.sam | samtools sort - ${f%.fastq.gz}.60K.4up.con.chrM.sort; samtools index ${f%.fastq.gz}.60K.4up.con.chrM.sort.bam; samtools idxstats ${f%.fastq.gz}.60K.4up.con.chrM.sort.bam; cd ..; done 

conda deavtivate 
conda activate longread_umi_Intel

cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101 

time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; usearch -fastx_uniques test.umi12.f.60K.u.c.94.4up.clusters.con.fix.chrM.fastq -fastaout test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.fasta -sizeout -minuniquesize 1 -relabel Unicon -strand both -tabbedout test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.tabbed -uc test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.unicon.txt; cd ..; pwd;  done



time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; grep ">" test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.fasta | sed 's/^>//g' > test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.unicon.title; cut -d ";" -f2 test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.unicon.title | sort | uniq -c > test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.unicon.title.count; cd ..; pwd; done

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	pwd
	
	awk -F ";|=" '{print $3;}' test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.unicon.title > test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.unicon.hist.txt
	wc -l test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.unicon.hist.txt
 
	cd ..
	pwd

done

conda deactivate
conda activate samtools_env_Intel
cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101 

time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; minimap2 -x map-ont -t 8 --frag=yes --secondary=yes -N 24 -p 0.8 -a -o ${f%.fastq.gz}.60K.4up.clusters.con.f.u.sam /Users/sw3203/Documents/Fun/Ref/Hg38/grch38_1kgmaj.fa test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.fasta; samtools view -bS ${f%.fastq.gz}.60K.4up.clusters.con.f.u.sam | samtools sort - ${f%.fastq.gz}.60K.4up.clusters.con.f.u.sort; samtools index ${f%.fastq.gz}.60K.4up.clusters.con.f.u.sort.bam; samtools idxstats ${f%.fastq.gz}.60K.4up.clusters.con.f.u.sort.bam; cd ..; done 
conda deactivate

conda activate longread_umi_Intel

cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101 


time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; mkdir cluster_unicon_60K_100; usearch -cluster_fast test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.fasta -id 1.0000 -centroids test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.fasta -uc test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.txt -sizein -sizeout -strand both -minsize 1 -sort size -clusters cluster_unicon_60K_100/test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.clusters; grep ">" test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.fasta | sed 's/^>//g' > test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.title; cut -d ";" -f2 test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.title | sort | uniq -c > test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.title.count; cd ..; pwd; done



conda deactivate
conda activate samtools_env_Intel

cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101 

time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; minimap2 -x map-ont -t 12 --frag=yes --secondary=yes -N 24 -p 0.8 -a -o ${f%.fastq.gz}.60K.4up.clusters.con.f.u.c.100.sam /Users/sw3203/Documents/Fun/Ref/Hg38/grch38_1kgmaj.fa test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.fasta; samtools view -bS ${f%.fastq.gz}.60K.4up.clusters.con.f.u.c.100.sam | samtools sort - ${f%.fastq.gz}.60K.4up.clusters.con.f.u.c.100.sort; samtools index ${f%.fastq.gz}.60K.4up.clusters.con.f.u.c.100.sort.bam; samtools idxstats ${f%.fastq.gz}.60K.4up.clusters.con.f.u.c.100.sort.bam; cd ..; done 

conda deactivate
conda activate longread_umi_Intel
cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101 
### consensus for Ave-1 clusters 

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd	
	echo "Generate consensus for each cluster"
	## assume the con.fix.fastq files had been generated at 4up 

	rm ../test.umi12.f.60K.u.c.94.Ave-1.clusters.con.fix.fastq
	touch ../test.umi12.f.60K.u.c.94.Ave-1.clusters.con.fix.fastq
	
	for f in $(cat ../test.umi12.f.60K.u.c.94.Ave-1.clusters.rename); do cat $f.60K.con.fix.fastq >> ../test.umi12.f.60K.u.c.94.Ave-1.clusters.con.fix.fastq; done 	
	
	wc -l ../test.umi12.f.60K.u.c.94.Ave-1.clusters.con.fix.fastq

	cd ../..
	pwd
	
done


### keep the consensus on chrM 
time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	
	pwd	
	
	awk 'BEGIN {keep = 0} /^@/ {keep = index($0, "chrM") > 0} {if (keep) print}' test.umi12.f.60K.u.c.94.Ave-1.clusters.con.fix.fastq > test.umi12.f.60K.u.c.94.Ave-1.clusters.con.fix.chrM.fastq
	
	wc -l test.umi12.f.60K.u.c.94.Ave-1.clusters.con.fix.fastq
	wc -l test.umi12.f.60K.u.c.94.Ave-1.clusters.con.fix.chrM.fastq
	
	cd ..
	pwd
	
done

time for f in $(cat nano.input.file)
do 
	
	cd ${f%.fastq.gz}; 

	
	wc -l test.umi12.f.60K.u.c.94.Ave-1.clusters.con.fix.chrM.fastq
	
	cd ..

	
done


### map the consensus to hg38, cluster with USEARCH

conda deactivate
conda activate samtools_env_Intel
cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101 

time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; minimap2 -x map-ont -t 8 --frag=yes --secondary=yes -N 24 -p 0.8 -a -o ${f%.fastq.gz}.60K.Ave-1.con.chrM.sam /Users/sw3203/Documents/Fun/Ref/Hg38/grch38_1kgmaj.fa test.umi12.f.60K.u.c.94.Ave-1.clusters.con.fix.chrM.fastq; samtools view -bS ${f%.fastq.gz}.60K.Ave-1.con.chrM.sam | samtools sort - ${f%.fastq.gz}.60K.Ave-1.con.chrM.sort; samtools index ${f%.fastq.gz}.60K.Ave-1.con.chrM.sort.bam; samtools idxstats ${f%.fastq.gz}.60K.Ave-1.con.chrM.sort.bam; cd ..; done 
conda deactivate

conda activate longread_umi_Intel

cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101 

time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; usearch -fastx_uniques test.umi12.f.60K.u.c.94.Ave-1.clusters.con.fix.chrM.fastq -fastaout test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.fasta -sizeout -minuniquesize 1 -relabel Unicon -strand both -tabbedout test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.tabbed -uc test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.unicon.txt; cd ..; pwd;  done



time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; grep ">" test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.fasta | sed 's/^>//g' > test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.unicon.title; cut -d ";" -f2 test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.unicon.title | sort | uniq -c > test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.unicon.title.count; cd ..; pwd; done

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	pwd
	
	awk -F ";|=" '{print $3;}' test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.unicon.title > test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.unicon.hist.txt
	wc -l test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.unicon.hist.txt
 
	cd ..
	pwd

done

conda deactivate
conda activate samtools_env_Intel
cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101 

time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; minimap2 -x map-ont -t 8 --frag=yes --secondary=yes -N 24 -p 0.8 -a -o ${f%.fastq.gz}.60K.Ave-1.clusters.con.f.u.sam /Users/sw3203/Documents/Fun/Ref/Hg38/grch38_1kgmaj.fa test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.fasta; samtools view -bS ${f%.fastq.gz}.60K.Ave-1.clusters.con.f.u.sam | samtools sort - ${f%.fastq.gz}.60K.Ave-1.clusters.con.f.u.sort; samtools index ${f%.fastq.gz}.60K.Ave-1.clusters.con.f.u.sort.bam; samtools idxstats ${f%.fastq.gz}.60K.Ave-1.clusters.con.f.u.sort.bam; cd ..; done 
conda deactivate

conda activate longread_umi_Intel

cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101 


time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; mkdir cluster_unicon_60K_100; usearch -cluster_fast test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.fasta -id 1.0000 -centroids test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.fasta -uc test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.txt -sizein -sizeout -strand both -minsize 1 -sort size -clusters cluster_unicon_60K_100/test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.clusters; grep ">" test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.fasta | sed 's/^>//g' > test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.title; cut -d ";" -f2 test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.title | sort | uniq -c > test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.title.count; cd ..; pwd; done


##### map the centroids -id 1.000
conda deactivate
conda activate samtools_env_Intel
cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101 

time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; minimap2 -x map-ont -t 12 --frag=yes --secondary=yes -N 24 -p 0.8 -a -o ${f%.fastq.gz}.60K.Ave-1.clusters.con.f.u.c.100.sam /Users/sw3203/Documents/Fun/Ref/Hg38/grch38_1kgmaj.fa test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.fasta; samtools view -bS ${f%.fastq.gz}.Ave-1.clusters.con.f.u.c.100.sam | samtools sort - ${f%.fastq.gz}.Ave-1.clusters.con.f.u.c.100.sort; samtools index ${f%.fastq.gz}.Ave-1.clusters.con.f.u.c.100.sort.bam; samtools idxstats ${f%.fastq.gz}.Ave-1.clusters.con.f.u.c.100.sort.bam; cd ..; done 
conda deactivate

conda activate longread_umi_Intel

cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101 

time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; mkdir cluster_unicon_60K_1mis; usearch -cluster_fast test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.fasta -id 0.9999 -centroids test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.1mis.fasta -uc test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.1mis.txt -sizein -sizeout -strand both -minsize 1 -sort size -clusters cluster_unicon_60K_1mis/test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.1mis.clusters; grep ">" test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.1mis.fasta | sed 's/^>//g' > test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.1mis.title; cut -d ";" -f2 test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.1mis.title | sort | uniq -c > test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.1mis.title.count; cd ..; pwd; done




### consensus for Ave-2 cluster 



time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd	
	echo "Generate consensus for each cluster"
	## assume the con.fix.fastq files had been generated at 4up 

	rm ../test.umi12.f.60K.u.c.94.Ave-2.clusters.con.fix.fastq
	touch ../test.umi12.f.60K.u.c.94.Ave-2.clusters.con.fix.fastq
	
	for f in $(cat ../test.umi12.f.60K.u.c.94.Ave-2.clusters.rename); do cat $f.60K.con.fix.fastq >> ../test.umi12.f.60K.u.c.94.Ave-2.clusters.con.fix.fastq; done 	
	
	wc -l ../test.umi12.f.60K.u.c.94.Ave-2.clusters.con.fix.fastq

	cd ../..
	pwd
	
done


### keep the consensus on chrM 
time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	
	pwd	
	
	awk 'BEGIN {keep = 0} /^@/ {keep = index($0, "chrM") > 0} {if (keep) print}' test.umi12.f.60K.u.c.94.Ave-2.clusters.con.fix.fastq > test.umi12.f.60K.u.c.94.Ave-2.clusters.con.fix.chrM.fastq
	
	wc -l test.umi12.f.60K.u.c.94.Ave-2.clusters.con.fix.fastq
	wc -l test.umi12.f.60K.u.c.94.Ave-2.clusters.con.fix.chrM.fastq
	
	cd ..
	pwd
	
done

### map the consensus to hg38, cluster with USEARCH

conda deactivate
conda activate samtools_env_Intel
cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101 

time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; minimap2 -x map-ont -t 8 --frag=yes --secondary=yes -N 24 -p 0.8 -a -o ${f%.fastq.gz}.60K.Ave-2.con.chrM.sam /Users/sw3203/Documents/Fun/Ref/Hg38/grch38_1kgmaj.fa test.umi12.f.60K.u.c.94.Ave-2.clusters.con.fix.chrM.fastq; samtools view -bS ${f%.fastq.gz}.60K.Ave-2.con.chrM.sam | samtools sort - ${f%.fastq.gz}.60K.Ave-2.con.chrM.sort; samtools index ${f%.fastq.gz}.60K.Ave-2.con.chrM.sort.bam; samtools idxstats ${f%.fastq.gz}.60K.Ave-2.con.chrM.sort.bam; cd ..; done 
conda deactivate
conda activate longread_umi_Intel

cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101 

time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; usearch -fastx_uniques test.umi12.f.60K.u.c.94.Ave-2.clusters.con.fix.chrM.fastq -fastaout test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.fasta -sizeout -minuniquesize 1 -relabel Unicon -strand both -tabbedout test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.tabbed -uc test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.unicon.txt; cd ..; pwd;  done



time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; grep ">" test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.fasta | sed 's/^>//g' > test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.unicon.title; cut -d ";" -f2 test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.unicon.title | sort | uniq -c > test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.unicon.title.count; cd ..; pwd; done

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	pwd
	
	awk -F ";|=" '{print $3;}' test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.unicon.title > test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.unicon.hist.txt
	wc -l test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.unicon.hist.txt
 
	cd ..
	pwd

done

conda deactivate
conda activate samtools_env_Intel
cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101 

time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; minimap2 -x map-ont -t 8 --frag=yes --secondary=yes -N 24 -p 0.8 -a -o ${f%.fastq.gz}.60K.Ave-2.clusters.con.f.u.sam /Users/sw3203/Documents/Fun/Ref/Hg38/grch38_1kgmaj.fa test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.fasta; samtools view -bS ${f%.fastq.gz}.60K.Ave-2.clusters.con.f.u.sam | samtools sort - ${f%.fastq.gz}.60K.Ave-2.clusters.con.f.u.sort; samtools index ${f%.fastq.gz}.60K.Ave-2.clusters.con.f.u.sort.bam; samtools idxstats ${f%.fastq.gz}.60K.Ave-2.clusters.con.f.u.sort.bam; cd ..; done 
conda deactivate
conda activate longread_umi_Intel

cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101 


time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; mkdir cluster_unicon_60K_100; usearch -cluster_fast test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.fasta -id 1.0000 -centroids test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.fasta -uc test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.txt -sizein -sizeout -strand both -minsize 1 -sort size -clusters cluster_unicon_60K_100/test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.clusters; grep ">" test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.fasta | sed 's/^>//g' > test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.title; cut -d ";" -f2 test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.title | sort | uniq -c > test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.title.count; cd ..; pwd; done


##### map the centroids -id 1.000
conda deactivate
conda activate samtools_env_Intel
cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101 

time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; minimap2 -x map-ont -t 12 --frag=yes --secondary=yes -N 24 -p 0.8 -a -o ${f%.fastq.gz}.60K.Ave-2.clusters.con.f.u.c.100.sam /Users/sw3203/Documents/Fun/Ref/Hg38/grch38_1kgmaj.fa test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.fasta; samtools view -bS ${f%.fastq.gz}.60K.Ave-2.clusters.con.f.u.c.100.sam | samtools sort - ${f%.fastq.gz}.60K.Ave-2.clusters.con.f.u.c.100.sort; samtools index ${f%.fastq.gz}.60K.Ave-2.clusters.con.f.u.c.100.sort.bam; samtools idxstats ${f%.fastq.gz}.60K.Ave-2.clusters.con.f.u.c.100.sort.bam; cd ..; done 
conda deactivate
conda activate longread_umi_Intel

cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101 

time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; mkdir cluster_unicon_60K_1mis; usearch -cluster_fast test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.fasta -id 0.9999 -centroids test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.1mis.fasta -uc test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.1mis.txt -sizein -sizeout -strand both -minsize 1 -sort size -clusters cluster_unicon_60K_1mis/test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.1mis.clusters; grep ">" test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.1mis.fasta | sed 's/^>//g' > test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.1mis.title; cut -d ";" -f2 test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.1mis.title | sort | uniq -c > test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.1mis.title.count; cd ..; pwd; done





####################  variant calling for 4up 
######   4up first 

time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; grep ">" test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.fasta | sed 's/^>//g' > test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.title; cut -d ";" -f1 test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.title > test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.title.ID; cut -d ";" -f2 test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.title | sort | uniq -c > test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.title.count; cd ..; pwd; done


time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; awk -F "\t|;" 'FNR==NR{a[$0];next}($9 in a)' test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.title.ID test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.txt | cut -f2 | sort -k1n | uniq > test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.clusters; wc -l test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.title.ID; wc -l test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.clusters; cd ..; pwd; done

conda deactivate
#### results
#### end of results
conda activate samtools_env_Intel
cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101

################# sam for each cluster, 4up  ######################
#########not the longread_umi_Intel samtool

#parallel
 time for f in $(cat nano.input.file); do  echo $f; cd ${f%.fastq.gz}; gunzip ${f%.fastq.gz}.60K.4up.con.chrM.sam.gz; cd ..; done


echo "4up;"

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	sed 's/^/test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.clusters/g' test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.clusters > test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.clusters.rename
	 
	cd cluster_unicon_60K_100
	pwd;
 
	time parallel -j 28 ' grep "^>" {} | cut -d ";" -f1 | sed "s/^>//g" > {}.ID ' < ../test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.clusters.rename;
 
	time parallel -j 28 " awk 'FNR==NR{a[\$0];next}(\$2 in a)' {}.ID ../test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.tabbed > {}.seq; cut -f1 {}.seq > {}.seq.ID; awk 'FNR==NR{a[\$0];next}(\$1 in a)' {}.seq.ID ../${f%.fastq.gz}.60K.4up.con.chrM.sam > {}.sam; samtools view -bT /Users/sw3203/Documents/Fun/Ref/Hg38/grch38_1kgmaj.fa {}.sam | samtools sort - {}; samtools index {}.bam" < ../test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.clusters.rename;
 
	cd ../..
	pwd

done



for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_unicon_60K_100;
	pwd
	
	time parallel -j 28 ' bcftools mpileup -Ov -d 10000 -L 10000 -a "FORMAT/AD,FORMAT/DP" -f /Users/sw3203/Documents/Fun/Ref/Hg38/grch38_1kgmaj.fa {}.bam > {}.vcf ' < ../test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.clusters.rename;
		

	time parallel -j 28 'bcftools view -i "FORMAT/AD[0:1]/FORMAT/DP>0.5 && INFO/DP>0" {}.vcf > {}.50.4.vcf; bgzip -c {}.50.4.vcf > {}.50.4.vcf.gz; tabix {}.50.4.vcf.gz; grep -v "#" {}.50.4.vcf > {}.50.vcf ' < ../test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.clusters.rename;
	
	 
	 cd ../..
done



for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_unicon_60K_100;
	pwd
	
	time parallel -j 28 'bcftools view -i "FORMAT/AD[0:1]/FORMAT/DP>0.6 && INFO/DP>0" {}.vcf > {}.60.2.vcf; bgzip -c {}.60.2.vcf > {}.60.2.vcf.gz; tabix {}.60.2.vcf.gz; grep -v "#" {}.60.2.vcf > {}.60.vcf ' < ../test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.clusters.rename;
	 
	 cd ../..
done


time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_unicon_60K_100;
	pwd
	
	rm 60.vcf.count
	rm call.60.vcf
	touch 60.vcf.count
	touch call.60.vcf
	
	for f in $(cat ../test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.clusters.rename); do wc -l  $f.60.vcf >> 60.vcf.count; cat $f.60.vcf >> call.60.vcf;  done
	awk '$1>0 {print $2;}' 60.vcf.count > 60.vcf.count.1
	cd ../..
	pwd
	
done


time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_unicon_60K_100;
	pwd	
	awk '$1>0 {print $2;}' 60.vcf.count > 60.vcf.count.1
	for f in $(cat 60.vcf.count.1); do echo $f; cat $f; done > 60.vcf.count.1.input
#	awk '$2~"3274" {print $0;}' 60.vcf.count.1.input
    cp  60.vcf.count.1.input ../test.umi12.f.60K.u.c.94.4up.clusters.con.f.u.c.100.60.vcf.count.1.input
	cut -f2 call.60.vcf | sort | uniq -c | sort -k1nr -k2nr
	cd ../..
	pwd
	
done

############### VARIANT calling for Ave-2 

cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101


time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; grep ">" test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.fasta | sed 's/^>//g' > test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.title; cut -d ";" -f1 test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.title > test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.title.ID; cut -d ";" -f2 test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.title | sort | uniq -c > test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.title.count; cd ..; pwd; done


time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; awk -F "\t|;" 'FNR==NR{a[$0];next}($9 in a)' test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.title.ID test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.txt | cut -f2 | sort -k1n | uniq > test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.clusters; wc -l test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.title.ID; wc -l test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.clusters; cd ..; pwd; done


#### results
#### end of results


cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_1012
#conda activate longread_umi_Intel
################# sam for each cluster, 4up  ######################
#########not the longread_umi_Intel samtool
#conda deactivate
#parallel
 time for f in $(cat nano.input.file); do  echo $f; cd ${f%.fastq.gz}; gunzip ${f%.fastq.gz}.Ave-2.con.chrM.sam.gz; cd ..; done


echo "Ave-2;"

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	sed 's/^/test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.clusters/g' test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.clusters > test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.clusters.rename
	 
	cd cluster_unicon_60K_100
	pwd;
 
	time parallel -j 28 ' grep "^>" {} | cut -d ";" -f1 | sed "s/^>//g" > {}.ID ' < ../test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.clusters.rename;
 
	time parallel -j 28 " awk 'FNR==NR{a[\$0];next}(\$2 in a)' {}.ID ../test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.tabbed > {}.seq; cut -f1 {}.seq > {}.seq.ID; awk 'FNR==NR{a[\$0];next}(\$1 in a)' {}.seq.ID ../${f%.fastq.gz}.60K.Ave-2.con.chrM.sam > {}.sam; samtools view -bT /Users/sw3203/Documents/Fun/Ref/Hg38/grch38_1kgmaj.fa {}.sam | samtools sort - {}; samtools index {}.bam" < ../test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.clusters.rename;
 
	cd ../..
	pwd

done



for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_unicon_60K_100;
	pwd
	
	time parallel -j 28 ' bcftools mpileup -Ov -d 10000 -L 10000 -a "FORMAT/AD,FORMAT/DP" -f /Users/sw3203/Documents/Fun/Ref/Hg38/grch38_1kgmaj.fa {}.bam > {}.vcf ' < ../test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.clusters.rename;
		

	time parallel -j 28 'bcftools view -i "FORMAT/AD[0:1]/FORMAT/DP>0.5 && INFO/DP>0" {}.vcf > {}.50.4.vcf; bgzip -c {}.50.4.vcf > {}.50.4.vcf.gz; tabix {}.50.4.vcf.gz; grep -v "#" {}.50.4.vcf > {}.50.vcf ' < ../test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.clusters.rename;
	
	 
	 cd ../..
done



for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_unicon_60K_100;
	pwd
	
	time parallel -j 28 'bcftools view -i "FORMAT/AD[0:1]/FORMAT/DP>0.6 && INFO/DP>0" {}.vcf > {}.60.2.vcf; bgzip -c {}.60.2.vcf > {}.60.2.vcf.gz; tabix {}.60.2.vcf.gz; grep -v "#" {}.60.2.vcf > {}.60.vcf ' < ../test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.clusters.rename;
	 
	 cd ../..
done


time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_unicon_60K_100;
	pwd
	
	rm 60.vcf.count
	rm call.60.vcf
	touch 60.vcf.count
	touch call.60.vcf
	
	for f in $(cat ../test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.clusters.rename); do wc -l  $f.60.vcf >> 60.vcf.count; cat $f.60.vcf >> call.60.vcf;  done
	awk '$1>0 {print $2;}' 60.vcf.count > 60.vcf.count.1
	cd ../..
	pwd
	
done


time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_unicon_60K_100;
	pwd	
	awk '$1>0 {print $2;}' 60.vcf.count > 60.vcf.count.1
	for f in $(cat 60.vcf.count.1); do echo $f; cat $f; done > 60.vcf.count.1.input
#	awk '$2~"3274" {print $0;}' 60.vcf.count.1.input
    cp  60.vcf.count.1.input ../test.umi12.f.60K.u.c.94.Ave-2.clusters.con.f.u.c.100.60.vcf.count.1.input
	cut -f2 call.60.vcf | sort | uniq -c | sort -k1nr -k2nr
	cd ../..
	pwd
	
done

################## variant for Ave-1


cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101


time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; grep ">" test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.fasta | sed 's/^>//g' > test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.title; cut -d ";" -f1 test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.title > test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.title.ID; cut -d ";" -f2 test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.title | sort | uniq -c > test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.title.count; cd ..; pwd; done


time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; awk -F "\t|;" 'FNR==NR{a[$0];next}($9 in a)' test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.title.ID test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.txt | cut -f2 | sort -k1n | uniq > test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.clusters; wc -l test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.title.ID; wc -l test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.clusters; cd ..; pwd; done

time for f in $(cat nano.input.file); do  cd ${f%.fastq.gz};  wc -l test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.clusters; cd ..; done


#### results
#### end of results
cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101
#conda activate longread_umi_Intel
################# sam for each cluster, 5up  ######################
#########not the longread_umi_Intel samtool
#conda deactivate
#parallel
# time for f in $(cat nano.input.file); do  echo $f; cd ${f%.fastq.gz}; gunzip ${f%.fastq.gz}.Ave-1.con.chrM.sam.gz; cd ..; done


echo "5up;"

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	sed 's/^/test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.clusters/g' test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.clusters > test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.clusters.rename
	 
	cd cluster_unicon_60K_100
	pwd;
 
	time parallel -j 28 ' grep "^>" {} | cut -d ";" -f1 | sed "s/^>//g" > {}.ID ' < ../test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.clusters.rename;
 
	time parallel -j 28 " awk 'FNR==NR{a[\$0];next}(\$2 in a)' {}.ID ../test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.tabbed > {}.seq; cut -f1 {}.seq > {}.seq.ID; awk 'FNR==NR{a[\$0];next}(\$1 in a)' {}.seq.ID ../${f%.fastq.gz}.60K.Ave-1.con.chrM.sam > {}.sam; samtools view -bT /Users/sw3203/Documents/Fun/Ref/Hg38/grch38_1kgmaj.fa {}.sam | samtools sort - {}; samtools index {}.bam" < ../test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.clusters.rename;
 
	cd ../..
	pwd

done



for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_unicon_60K_100;
	pwd
	
	time parallel -j 28 ' bcftools mpileup -Ov -d 10000 -L 10000 -a "FORMAT/AD,FORMAT/DP" -f /Users/sw3203/Documents/Fun/Ref/Hg38/grch38_1kgmaj.fa {}.bam > {}.vcf ' < ../test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.clusters.rename;
		

	time parallel -j 28 'bcftools view -i "FORMAT/AD[0:1]/FORMAT/DP>0.5 && INFO/DP>0" {}.vcf > {}.50.4.vcf; bgzip -c {}.50.4.vcf > {}.50.4.vcf.gz; tabix {}.50.4.vcf.gz; grep -v "#" {}.50.4.vcf > {}.50.vcf ' < ../test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.clusters.rename;
	
	 
	 cd ../..
done



for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_unicon_60K_100;
	pwd
	
	time parallel -j 28 'bcftools view -i "FORMAT/AD[0:1]/FORMAT/DP>0.6 && INFO/DP>0" {}.vcf > {}.60.2.vcf; bgzip -c {}.60.2.vcf > {}.60.2.vcf.gz; tabix {}.60.2.vcf.gz; grep -v "#" {}.60.2.vcf > {}.60.vcf ' < ../test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.clusters.rename;
	 
	 cd ../..
done


time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_unicon_60K_100;
	pwd
	
	rm 60.vcf.count
	rm call.60.vcf
	touch 60.vcf.count
	touch call.60.vcf
	
	for f in $(cat ../test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.clusters.rename); do wc -l  $f.60.vcf >> 60.vcf.count; cat $f.60.vcf >> call.60.vcf;  done
	awk '$1>0 {print $2;}' 60.vcf.count > 60.vcf.count.1
	cd ../..
	pwd
	
done


time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_unicon_60K_100;
	pwd	
	awk '$1>0 {print $2;}' 60.vcf.count > 60.vcf.count.1
	for f in $(cat 60.vcf.count.1); do echo $f; cat $f; done > 60.vcf.count.1.input
#	awk '$2~"3274" {print $0;}' 60.vcf.count.1.input
    cp  60.vcf.count.1.input ../test.umi12.f.60K.u.c.94.Ave-1.clusters.con.f.u.c.100.60.vcf.count.1.input
	cut -f2 call.60.vcf | sort | uniq -c | sort -k1nr -k2nr
	cd ../..
	pwd
	
done
conda deactivate

######## cons for Ave-Cal
### consensus for Ave-cal clusters 

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd	
	echo "Generate consensus for each cluster"
	## assume the con.fix.fastq files had been generated at 4up 

	rm ../test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.fix.fastq
	touch ../test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.fix.fastq
	
	for f in $(cat ../test.umi12.f.60K.u.c.94.Ave-cal.clusters.rename); do cat $f.60K.con.fix.fastq >> ../test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.fix.fastq; done 	
	
	wc -l ../test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.fix.fastq

	cd ../..
	pwd
	
done


### keep the consensus on chrM 
time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	
	pwd	
	
	awk 'BEGIN {keep = 0} /^@/ {keep = index($0, "chrM") > 0} {if (keep) print}' test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.fix.fastq > test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.fix.chrM.fastq
	
	wc -l test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.fix.fastq
	wc -l test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.fix.chrM.fastq
	
	cd ..
	pwd
	
done

time for f in $(cat nano.input.file)
do 
	
	cd ${f%.fastq.gz}; 

	
	wc -l test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.fix.chrM.fastq
	
	cd ..

	
done


### map the consensus to hg38, cluster with USEARCH

conda deactivate
conda activate samtools_env_Intel
cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101 

time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; minimap2 -x map-ont -t 8 --frag=yes --secondary=yes -N 24 -p 0.8 -a -o ${f%.fastq.gz}.60K.Ave-cal.con.chrM.sam /Users/sw3203/Documents/Fun/Ref/Hg38/grch38_1kgmaj.fa test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.fix.chrM.fastq; samtools view -bS ${f%.fastq.gz}.60K.Ave-cal.con.chrM.sam | samtools sort - ${f%.fastq.gz}.60K.Ave-cal.con.chrM.sort; samtools index ${f%.fastq.gz}.60K.Ave-cal.con.chrM.sort.bam; samtools idxstats ${f%.fastq.gz}.60K.Ave-cal.con.chrM.sort.bam; cd ..; done 
conda deactivate

echo "Ave-cal.con.chrM.sort.bam done"

conda activate longread_umi_Intel

cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101 

time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; usearch -fastx_uniques test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.fix.chrM.fastq -fastaout test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.f.u.fasta -sizeout -minuniquesize 1 -relabel Unicon -strand both -tabbedout test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.f.u.tabbed -uc test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.f.u.unicon.txt; cd ..; pwd;  done



time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; grep ">" test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.f.u.fasta | sed 's/^>//g' > test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.f.u.unicon.title; cut -d ";" -f2 test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.f.u.unicon.title | sort | uniq -c > test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.f.u.unicon.title.count; cd ..; pwd; done

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	pwd
	
	awk -F ";|=" '{print $3;}' test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.f.u.unicon.title > test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.f.u.unicon.hist.txt
	wc -l test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.f.u.unicon.hist.txt
 
	cd ..
	pwd

done

conda deactivate
conda activate samtools_env_Intel
cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101 

time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; minimap2 -x map-ont -t 8 --frag=yes --secondary=yes -N 24 -p 0.8 -a -o ${f%.fastq.gz}.60K.Ave-cal.clusters.con.f.u.sam /Users/sw3203/Documents/Fun/Ref/Hg38/grch38_1kgmaj.fa test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.f.u.fasta; samtools view -bS ${f%.fastq.gz}.60K.Ave-cal.clusters.con.f.u.sam | samtools sort - ${f%.fastq.gz}.60K.Ave-cal.clusters.con.f.u.sort; samtools index ${f%.fastq.gz}.60K.Ave-cal.clusters.con.f.u.sort.bam; samtools idxstats ${f%.fastq.gz}.60K.Ave-cal.clusters.con.f.u.sort.bam; cd ..; done 
conda deactivate

conda activate longread_umi_Intel

cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101 


time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; mkdir cluster_unicon_60K_100; usearch -cluster_fast test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.f.u.fasta -id 1.0000 -centroids test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.f.u.c.100.fasta -uc test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.f.u.c.100.txt -sizein -sizeout -strand both -minsize 1 -sort size -clusters cluster_unicon_60K_100/test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.f.u.c.100.clusters; grep ">" test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.f.u.c.100.fasta | sed 's/^>//g' > test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.f.u.c.100.title; cut -d ";" -f2 test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.f.u.c.100.title | sort | uniq -c > test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.f.u.c.100.title.count; cd ..; pwd; done


##### map the centroids -id 1.000
conda deactivate
conda activate samtools_env_Intel
cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101 

time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; minimap2 -x map-ont -t 12 --frag=yes --secondary=yes -N 24 -p 0.8 -a -o ${f%.fastq.gz}.60K.Ave-cal.clusters.con.f.u.c.100.sam /Users/sw3203/Documents/Fun/Ref/Hg38/grch38_1kgmaj.fa test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.f.u.c.100.fasta; samtools view -bS ${f%.fastq.gz}.Ave-cal.clusters.con.f.u.c.100.sam | samtools sort - ${f%.fastq.gz}.Ave-cal.clusters.con.f.u.c.100.sort; samtools index ${f%.fastq.gz}.Ave-cal.clusters.con.f.u.c.100.sort.bam; samtools idxstats ${f%.fastq.gz}.Ave-cal.clusters.con.f.u.c.100.sort.bam; cd ..; done 
conda deactivate

conda activate longread_umi_Intel

cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101 

time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; mkdir cluster_unicon_60K_1mis; usearch -cluster_fast test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.f.u.fasta -id 0.9999 -centroids test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.f.u.c.1mis.fasta -uc test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.f.u.c.1mis.txt -sizein -sizeout -strand both -minsize 1 -sort size -clusters cluster_unicon_60K_1mis/test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.f.u.c.1mis.clusters; grep ">" test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.f.u.c.1mis.fasta | sed 's/^>//g' > test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.f.u.c.1mis.title; cut -d ";" -f2 test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.f.u.c.1mis.title | sort | uniq -c > test.umi12.f.60K.u.c.94.Ave-cal.clusters.con.f.u.c.1mis.title.count; cd ..; pwd; done




time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; gzip *.sam; cd ..; done

time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; tar -zcvf ${f%.fastq.gz}.LRUMI_out.tar.gz ${f%.fastq.gz}.LRUMI_out; rm -r ${f%.fastq.gz}.LRUMI_out; cd ..; done

echo "end part 3 analysis"
date

