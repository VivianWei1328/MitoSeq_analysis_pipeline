#######################  MT-16 Mito duplex library, 60K first #############
###### this is for running on silicon chip 
#!/bin/bash
#!/usr/bin/perl

### before running, replace $PATHs to actual value
### Need conda environments for longread_umi_Intel and samtools_env_Intel
###### This pipeline use samtools 1.2 in samtools_env_Intel
###### samtools 1.21 in regular and samtools_env


cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101
pwd
source /opt/miniconda3/etc/profile.d/conda.sh

conda activate longread_umi
conda deactivate

#parallel


#time for f in $(cat input.file); do echo $f; cd ${f%.fastq.gz}; pwd; gzip *.fastq; cd ..; pwd; done

####### end of skip  ############################

###### need to use the sed under longread_umi environment ######
 
conda activate longread_umi
conda info


############ retrieve the ori sequences for UMIs. 90% similarity #######
########## MS method #########
###### the final cluster count in the title of centroid fasta file

time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; sed 's/^>\([^>]*\)>.*/>\1/' umi12f.fa > test.umi12.f.fasta; head -120000 test.umi12.f.fasta > test.umi12.f.60K.fasta; usearch -fastx_uniques test.umi12.f.60K.fasta -fastaout test.umi12.f.60K.u.fasta -sizeout -minuniquesize 1 -relabel umi -strand both -tabbedout test.umi12.f.60K.u.tabbed; mkdir cluster_90_60K; usearch -cluster_fast test.umi12.f.60K.u.fasta -id 0.90 -centroids test.umi12.f.60K.u.c.fasta -uc test.umi12.f.60K.u.c.txt -sizein -sizeout  -strand both -minsize 1 -clusters cluster_90_60K/test.umi12.f.60K.u.c.clusters; grep ">" test.umi12.f.60K.u.c.fasta | sed 's/^>//g' > test.umi12.f.60K.u.c.title; cut -d ";" -f2 test.umi12.f.60K.u.c.title | sort | uniq -c > test.umi12.f.60K.u.c.title.count; cd ..; pwd;  done

#time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; usearch -fastx_uniques test.umi12.f.60K.fasta -fastaout test.umi12.f.60K.u.fasta -sizeout -minuniquesize 1 -relabel umi -strand both -tabbedout test.umi12.f.60K.u.tabbed; mkdir cluster_90_60K; usearch -cluster_fast test.umi12.f.60K.u.fasta -id 0.90 -centroids test.umi12.f.60K.u.c.fasta -uc test.umi12.f.60K.u.c.txt -sizein -sizeout  -strand both -minsize 1 -clusters cluster_90_60K/test.umi12.f.60K.u.c.clusters; grep ">" test.umi12.f.60K.u.c.fasta | sed 's/^>//g' > test.umi12.f.60K.u.c.title; cut -d ";" -f2 test.umi12.f.60K.u.c.title | sort | uniq -c > test.umi12.f.60K.u.c.title.count; cd ..; pwd;  done





######## result #########


######### end of results #############

###############################################################################
 ############ retrieve the ori sequences for UMIs. 94% similarity 2/36 #######
 ##### assuming the deunique is done at the 90% UMI #######

time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; rm test.umi12.f.60K.u.tabbed; usearch -fastx_uniques test.umi12.f.60K.fasta -fastaout test.umi12.f.60K.u.fasta -sizeout -minuniquesize 1 -relabel umi -strand both -tabbedout test.umi12.f.60K.u.tabbed; cd ..; pwd;  done



################
########### end of results ##########
############ use 90% for Mito clusters ############
############ use 94% for Mito clusters with MinKNOW 24.11 Q>30 ############


time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; mkdir cluster_94_60K; usearch -cluster_fast test.umi12.f.60K.u.fasta -id 0.94 -centroids test.umi12.f.60K.u.c.94.fasta -uc test.umi12.f.60K.u.c.94.txt -sizein -sizeout -strand both -minsize 1 -sort size -clusters cluster_94_60K/test.umi12.f.60K.u.c.94.clusters; grep ">" test.umi12.f.60K.u.c.94.fasta | sed 's/^>//g' > test.umi12.f.60K.u.c.94.title; cut -d ";" -f2 test.umi12.f.60K.u.c.94.title | sort | uniq -c > test.umi12.f.60K.u.c.94.title.count; cd ..; pwd; done


################# result  

########### end of result  
#########################################################
############   make hist 

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	pwd
	
	awk -F ";|=" '{print $3;}' test.umi12.f.60K.u.c.94.title > test.umi12.f.60K.u.c.94.title.hist.txt
	wc -l test.umi12.f.60K.u.c.94.title.hist.txt
 
	cd ..
	pwd

done


######   4up first 

time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; mkdir cluster_94_60K_4up; usearch -cluster_fast test.umi12.f.60K.u.fasta -id 0.94 -centroids test.umi12.f.60K.u.c.94.4up.fasta -uc test.umi12.f.60K.u.c.94.4up.txt -sizein -sizeout  -strand both -minsize 4 -sort size -clusters cluster_94_60K_4up/test.umi12.f.60K.u.c.94.4up.clusters; grep ">" test.umi12.f.60K.u.c.94.4up.fasta | sed 's/^>//g' > test.umi12.f.60K.u.c.94.4up.title; cut -d ";" -f1 test.umi12.f.60K.u.c.94.4up.title > test.umi12.f.60K.u.c.94.4up.title.ID; cut -d ";" -f2 test.umi12.f.60K.u.c.94.4up.title | sort | uniq -c > test.umi12.f.60K.u.c.94.4up.title.count; cd ..; pwd; done

#### results
#### end of results

time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; awk -F "\t|;" 'FNR==NR{a[$0];next}($9 in a)' test.umi12.f.60K.u.c.94.4up.title.ID test.umi12.f.60K.u.c.94.4up.txt | cut -f2 | sort -k1n | uniq > test.umi12.f.60K.u.c.94.4up.clusters; wc -l test.umi12.f.60K.u.c.94.4up.title.ID; wc -l test.umi12.f.60K.u.c.94.4up.clusters; cd ..; pwd; done


#### results
#### end of results
cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101

################# sam for each cluster, 4up  ######################
#########not the longread_umi samtools, samtools 1.9 in samtools_env_Intel

conda activate samtools_env_Intel
#parallel

echo "4up;"

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	sed 's/^/test.umi12.f.60K.u.c.94.4up.clusters/g' test.umi12.f.60K.u.c.94.4up.clusters > test.umi12.f.60K.u.c.94.4up.clusters.rename
	 
	cd cluster_94_60K_4up;
	pwd;
 
	time parallel -j 28 ' grep "^>" {} | cut -d ";" -f1 | sed "s/^>//g" > {}.ID ' < ../test.umi12.f.60K.u.c.94.4up.clusters.rename;
 
	time parallel -j 28 " awk 'FNR==NR{a[\$0];next}(\$2 in a)' {}.ID ../test.umi12.f.60K.u.tabbed > {}.seq; cut -f1 {}.seq > {}.seq.ID; awk 'FNR==NR{a[\$0];next}(\$1 in a)' {}.seq.ID ../${f%.fastq.gz}.mini.sam > {}.sam; samtools view -bT $PATH_to_reference/Ref/Hg38/grch38_1kgmaj.fa {}.sam | samtools sort - {}; samtools index {}.bam; rm {}.sam; rm {}.seq;" < ../test.umi12.f.60K.u.c.94.4up.clusters.rename;
 
	cd ../..
	pwd

done

for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd
	
	time parallel -j 28 ' bcftools mpileup -Ov -d 10000 -L 10000 -a "FORMAT/AD,FORMAT/DP" -f $PATH_to_reference/Ref/Hg38/grch38_1kgmaj.fa {}.bam > {}.vcf ' < ../test.umi12.f.60K.u.c.94.4up.clusters.rename;
		

	time parallel -j 28 'bcftools view -i "FORMAT/AD[0:1]/FORMAT/DP>0.5 && INFO/DP>4" {}.vcf > {}.50.4.vcf; bgzip -c {}.50.4.vcf > {}.50.4.vcf.gz; tabix {}.50.4.vcf.gz; grep -v "#" {}.50.4.vcf > {}.50.vcf ' < ../test.umi12.f.60K.u.c.94.4up.clusters.rename;
	
	 
	 cd ../..
done



for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd
	
	time parallel -j 28 'bcftools view -i "FORMAT/AD[0:1]/FORMAT/DP>0.6 && INFO/DP>2" {}.vcf > {}.60.2.vcf; bgzip -c {}.60.2.vcf > {}.60.2.vcf.gz; tabix {}.60.2.vcf.gz; grep -v "#" {}.60.2.vcf > {}.60.vcf ' < ../test.umi12.f.60K.u.c.94.4up.clusters.rename;
	 
	 cd ../..
done

	

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd
	
	rm 50.vcf.count
	rm call.50.vcf
	touch 50.vcf.count
	touch call.50.vcf
	
	for f in $(cat ../test.umi12.f.60K.u.c.94.4up.clusters.rename); do wc -l  $f.50.vcf >> 50.vcf.count; cat $f.50.vcf >> call.50.vcf;  done
	awk '$1>0 {print $2;}' 50.vcf.count > 50.vcf.count.1
	cd ../..
	pwd
	
done


time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd
	
	rm 60.vcf.count
	rm call.60.vcf
	touch 60.vcf.count
	touch call.60.vcf
	
	for f in $(cat ../test.umi12.f.60K.u.c.94.4up.clusters.rename); do wc -l  $f.60.vcf >> 60.vcf.count; cat $f.60.vcf >> call.60.vcf;  done
	awk '$1>0 {print $2;}' 60.vcf.count > 60.vcf.count.1
	cd ../..
	pwd
	
done


time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd	
#	awk '$1>0 {print $2;}' 60.vcf.count > 60.vcf.count.1
	for f in $(cat 60.vcf.count.1); do echo $f; cat $f; done > 60.vcf.count.1.input
#	awk '$2~"3274" {print $0;}' 60.vcf.count.1.input
	cut -f2 call.60.vcf | sort | uniq -c | sort -k1nr -k2nr
	cd ../..
	pwd
	
done


time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	
	tar -zcvf cluster_94_60K.tar.gz cluster_94_60K
	rm -r cluster_94_60K

	cd ..
	pwd
	
done

conda deactivate 

########################################################################## 
############## cutoff at 94% avg-1 size ################

####### assuming the 4up is done ###################################
######### change the -minsize to the average coverage ###########

cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101
pwd

echo " 94% avg-1 size "
conda activate longread_umi_Intel


time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; usearch -cluster_fast test.umi12.f.60K.u.fasta -id 0.94 -centroids test.umi12.f.60K.u.c.94.Ave-1.fasta -uc test.umi12.f.60K.u.c.94.Ave-1.txt -sizein -sizeout  -strand both -minsize 5 -sort size ; grep ">" test.umi12.f.60K.u.c.94.Ave-1.fasta | sed 's/^>//g' > test.umi12.f.60K.u.c.94.Ave-1.title; cut -d ";" -f1 test.umi12.f.60K.u.c.94.Ave-1.title > test.umi12.f.60K.u.c.94.Ave-1.title.ID; cut -d ";" -f2 test.umi12.f.60K.u.c.94.Ave-1.title | sort | uniq -c > test.umi12.f.60K.u.c.94.Ave-1.title.count; awk -F "\t|;" 'FNR==NR{a[$0];next}($9 in a)' test.umi12.f.60K.u.c.94.Ave-1.title.ID test.umi12.f.60K.u.c.94.Ave-1.txt | cut -f2 | sort -k1n | uniq > test.umi12.f.60K.u.c.94.Ave-1.clusters; wc -l test.umi12.f.60K.u.c.94.Ave-1.title.ID; wc -l test.umi12.f.60K.u.c.94.Ave-1.clusters; sed 's/^/test.umi12.f.60K.u.c.94.4up.clusters/g' test.umi12.f.60K.u.c.94.Ave-1.clusters > test.umi12.f.60K.u.c.94.Ave-1.clusters.rename; cd ..; pwd; done


   
time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	pwd 
	
	sed 's/^/test.umi12.f.60K.u.c.94.4up.clusters/g' test.umi12.f.60K.u.c.94.Ave-1.clusters > test.umi12.f.60K.u.c.94.Ave-1.clusters.rename;
	
	cd ..
done
	
conda deactivate

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd
	
	rm Ave-1.50.vcf.count
	rm Ave-1.call.50.vcf
	touch Ave-1.50.vcf.count
	touch Ave-1.call.50.vcf
	
	for f in $(cat ../test.umi12.f.60K.u.c.94.Ave-1.clusters.rename); do wc -l  $f.50.vcf >> Ave-1.50.vcf.count; cat $f.50.vcf >> Ave-1.call.50.vcf; done
	awk '$1>0 {print $2;}' Ave-1.50.vcf.count > Ave-1.50.vcf.count.1
	
	cd ../..
	
	pwd
	
done

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd	
	echo "Ave-1 output"
	#awk '$1>0 {print $2;}' 50.vcf.count > 50.vcf.count.1
	for f in $(cat Ave-1.50.vcf.count.1); do echo $f; cat $f; done > Ave-1.50.vcf.count.1.input
	#grep -B1 "53428" 50.vcf.count.1.input
	#grep -B1 "534289" Ave-1.50.vcf.count.1.input
	cut -f2 Ave-1.call.50.vcf | sort | uniq -c | sort -k1nr 
	cd ../..
	pwd
	
done

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd
	
	rm Ave-1.60.vcf.count
	rm Ave-1.call.60.vcf
	touch Ave-1.60.vcf.count
	touch Ave-1.call.60.vcf
	
	for f in $(cat ../test.umi12.f.60K.u.c.94.Ave-1.clusters.rename); do wc -l  $f.60.vcf >> Ave-1.60.vcf.count; cat $f.60.vcf >> Ave-1.call.60.vcf; done
	awk '$1>0 {print $2;}' Ave-1.60.vcf.count > Ave-1.60.vcf.count.1
	
	cd ../..
	
	pwd
	
done

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd	
	echo "Ave-1 output"
	#awk '$1>0 {print $2;}' 60.vcf.count > 60.vcf.count.1
	for f in $(cat Ave-1.60.vcf.count.1); do echo $f; cat $f; done > Ave-1.60.vcf.count.1.input
	#grep -B1 "53428" 60.vcf.count.1.input
#	grep -B1 "302" Ave-1.60.vcf.count.1.input
	cut -f2 Ave-1.call.60.vcf | sort | uniq -c | sort -k1nr -k2nr
	cd ../..
	pwd
	
done

rm $Run_name-13V-run1.60K.Ave-1.count
touch $Run_name-13V-run1.60K.Ave-1.count

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd	
	echo "Ave-1 output"
	#awk '$1>0 {print $2;}' 60.vcf.count > 60.vcf.count.1
#	for f in $(cat Ave-1.60.vcf.count.1); do echo $f; cat $f; done > Ave-1.60.vcf.count.1.input
	#grep -B1 "8597" 60.vcf.count.1.input
#	head -20 Ave-1.60.vcf.count.1.input
#	head -20 Ave-1.60.vcf.count
#	head -20 Ave-1.call.60.vcf
	
	cut -f2 Ave-1.call.60.vcf | sort | uniq -c | sort -k1nr -k2nr >> ../../$Run_name-13V-run1.60K.Ave-1.count
	cd ../..
	pwd
	
done




########################################################################## 
############## cutoff at 94% avg-2 size ################

####### assuming the 4up is done ###################################
######### change the -minsize to the average coverage ###########

cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101
pwd

echo " 94% avg-2 size "
conda activate longread_umi_Intel


time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; usearch -cluster_fast test.umi12.f.60K.u.fasta -id 0.94 -centroids test.umi12.f.60K.u.c.94.Ave-2.fasta -uc test.umi12.f.60K.u.c.94.Ave-2.txt -sizein -sizeout  -strand both -minsize 6 -sort size ; grep ">" test.umi12.f.60K.u.c.94.Ave-2.fasta | sed 's/^>//g' > test.umi12.f.60K.u.c.94.Ave-2.title; cut -d ";" -f1 test.umi12.f.60K.u.c.94.Ave-2.title > test.umi12.f.60K.u.c.94.Ave-2.title.ID; cut -d ";" -f2 test.umi12.f.60K.u.c.94.Ave-2.title | sort | uniq -c > test.umi12.f.60K.u.c.94.Ave-2.title.count; awk -F "\t|;" 'FNR==NR{a[$0];next}($9 in a)' test.umi12.f.60K.u.c.94.Ave-2.title.ID test.umi12.f.60K.u.c.94.Ave-2.txt | cut -f2 | sort -k1n | uniq > test.umi12.f.60K.u.c.94.Ave-2.clusters; wc -l test.umi12.f.60K.u.c.94.Ave-2.title.ID; wc -l test.umi12.f.60K.u.c.94.Ave-2.clusters; sed 's/^/test.umi12.f.60K.u.c.94.4up.clusters/g' test.umi12.f.60K.u.c.94.Ave-2.clusters > test.umi12.f.60K.u.c.94.Ave-2.clusters.rename; cd ..; pwd; done


   
time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	pwd 
	
	sed 's/^/test.umi12.f.60K.u.c.94.4up.clusters/g' test.umi12.f.60K.u.c.94.Ave-2.clusters > test.umi12.f.60K.u.c.94.Ave-2.clusters.rename;
	
	cd ..
done
	
conda deactivate

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd
	
	rm Ave-2.50.vcf.count
	rm Ave-2.call.50.vcf
	touch Ave-2.50.vcf.count
	touch Ave-2.call.50.vcf
	
	for f in $(cat ../test.umi12.f.60K.u.c.94.Ave-2.clusters.rename); do wc -l  $f.50.vcf >> Ave-2.50.vcf.count; cat $f.50.vcf >> Ave-2.call.50.vcf; done
	awk '$1>0 {print $2;}' Ave-2.50.vcf.count > Ave-2.50.vcf.count.1
	
	cd ../..
	
	pwd
	
done

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd	
	echo "Ave-2 output"
	#awk '$1>0 {print $2;}' 50.vcf.count > 50.vcf.count.1
	for f in $(cat Ave-2.50.vcf.count.1); do echo $f; cat $f; done > Ave-2.50.vcf.count.1.input
	#grep -B1 "53428" 50.vcf.count.1.input
	#grep -B1 "534289" Ave-2.50.vcf.count.1.input
	cut -f2 Ave-2.call.50.vcf | sort | uniq -c | sort -k1nr 
	cd ../..
	pwd
	test.umi12.f.60K.u.c.94.Ave-2.clusters
done

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd
	
	rm Ave-2.60.vcf.count
	rm Ave-2.call.60.vcf
	touch Ave-2.60.vcf.count
	touch Ave-2.call.60.vcf
	
	for f in $(cat ../test.umi12.f.60K.u.c.94.Ave-2.clusters.rename); do wc -l  $f.60.vcf >> Ave-2.60.vcf.count; cat $f.60.vcf >> Ave-2.call.60.vcf; done
	awk '$1>0 {print $2;}' Ave-2.60.vcf.count > Ave-2.60.vcf.count.1
	
	cd ../..
	
	pwd
	
done

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd	
	echo "Ave-2 output"
	#awk '$1>0 {print $2;}' 60.vcf.count > 60.vcf.count.1
	for f in $(cat Ave-2.60.vcf.count.1); do echo $f; cat $f; done > Ave-2.60.vcf.count.1.input
	#grep -B1 "53428" 60.vcf.count.1.input
#	grep -B1 "302" Ave-2.60.vcf.count.1.input
	cut -f2 Ave-2.call.60.vcf | sort | uniq -c | sort -k1nr -k2nr
	cd ../..
	pwd
	
done

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd	
	echo "Ave-2 output"
	#awk '$1>0 {print $2;}' 60.vcf.count > 60.vcf.count.1
#	for f in $(cat Ave-2.60.vcf.count.1); do echo $f; cat $f; done > Ave-2.60.vcf.count.1.input
	#grep -B1 "53428" 60.vcf.count.1.input
#	grep -B2 "302" Ave-2.60.vcf.count.1.input
#	cut -f2 Ave-2.call.60.vcf | sort | uniq -c | sort -k1nr -k2nr
	cd ../..
	pwd
	
done

rm $Run_name-13V-run1.60K.Ave-2.count
touch $Run_name-13V-run1.60K.Ave-2.count

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd	
	echo "Ave-2 output"
	#awk '$1>0 {print $2;}' 60.vcf.count > 60.vcf.count.1
#	for f in $(cat Ave-2.60.vcf.count.1); do echo $f; cat $f; done > Ave-2.60.vcf.count.1.input
	#grep -B1 "53428" 60.vcf.count.1.input
#	head -20 Ave-2.60.vcf.count.1.input
#	head -20 Ave-2.60.vcf.count
#	head -20 Ave-2.call.60.vcf
	
	cut -f2 Ave-2.call.60.vcf | sort | uniq -c | sort -k1nr -k2nr >> ../../$Run_name-13V-run1.60K.Ave-2.count
	cd ../..
	pwd
	
done


################### For 90 percentile  
##########################################################################
###################### analyze 2-3 reads ############################

conda activate longread_umi_Intel
conda info
cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101
pwd

##### use the cluster files in cluster_94_60K_4up ###########


time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; awk -F ";|=" ' $3>1 && $3<4 {print $1}' test.umi12.f.60K.u.c.94.title > test.umi12.f.60K.u.c.94.2-3.title.ID; awk -F "\t|;" 'FNR==NR{a[$0];next}($9 in a)' test.umi12.f.60K.u.c.94.2-3.title.ID test.umi12.f.60K.u.c.94.txt | cut -f2 | sort -k1n | uniq > test.umi12.f.60K.u.c.94.2-3.clusters; wc -l test.umi12.f.60K.u.c.94.2-3.title.ID; wc -l test.umi12.f.60K.u.c.94.2-3.clusters; sed 's/^/test.umi12.f.60K.u.c.94.4up.clusters/g' test.umi12.f.60K.u.c.94.2-3.clusters > test.umi12.f.60K.u.c.94.2-3.clusters.rename; cd ..; pwd; done

########## results ###################


########## end of results ######################


################# sam for each cluster ##############
#########not the longread_umi_Intel samtools 
conda deactivate

conda activate samtools_env_Intel

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 	 
	cd cluster_94_60K_4up;
	pwd;
 
	time parallel -j 28 ' grep "^>" {} | cut -d ";" -f1 | sed "s/^>//g" > {}.ID ' < ../test.umi12.f.60K.u.c.94.2-3.clusters.rename;
	
	time parallel -j 28 " awk 'FNR==NR{a[\$0];next}(\$2 in a)' {}.ID ../test.umi12.f.60K.u.tabbed > {}.seq; cut -f1 {}.seq > {}.seq.ID; awk 'FNR==NR{a[\$0];next}(\$1 in a)' {}.seq.ID ../${f%.fastq.gz}.mini.sam > {}.sam; samtools view -bT $PATH_to_reference/Ref/Hg38/grch38_1kgmaj.fa {}.sam | samtools sort - {}; samtools index {}.bam; rm {}.sam; rm {}.seq;" < ../test.umi12.f.60K.u.c.94.2-3.clusters.rename;
 
 
	cd ../..
	pwd

done
conda deactivate


conda activate longread_umi_Intel

for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd
	
	time parallel -j 28 ' bcftools mpileup -Ov -d 10000 -L 10000 -a "FORMAT/AD,FORMAT/DP" -f $PATH_to_reference/Ref/Hg38/grch38_1kgmaj.fa {}.bam > {}.vcf ' < ../test.umi12.f.60K.u.c.94.2-3.clusters.rename;
		
	time parallel -j 28 'bcftools view -i "FORMAT/AD[0:1]/FORMAT/DP>0.5 && INFO/DP>1" {}.vcf > {}.50.2.vcf; bgzip -c {}.50.2.vcf > {}.50.2.vcf.gz; tabix {}.50.2.vcf.gz; grep -v "#" {}.50.2.vcf > {}.50.vcf ' < ../test.umi12.f.60K.u.c.94.2-3.clusters.rename;
	
	time parallel -j 28 'bcftools view -i "FORMAT/AD[0:1]/FORMAT/DP>0.6 && INFO/DP>1" {}.vcf > {}.60.2.vcf; bgzip -c {}.60.2.vcf > {}.60.2.vcf.gz; tabix {}.60.2.vcf.gz; grep -v "#" {}.60.2.vcf > {}.60.vcf ' < ../test.umi12.f.60K.u.c.94.2-3.clusters.rename;
	 
	 cd ../..
done


time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd
	
	rm 2-3.60.vcf.count
	rm call.2-3.60.vcf
	touch 2-3.60.vcf.count
	touch call.2-3.60.vcf
	
	for f in $(cat ../test.umi12.f.60K.u.c.94.2-3.clusters.rename); do wc -l  $f.60.vcf >> 2-3.60.vcf.count; cat $f.60.vcf >> call.2-3.60.vcf;  done
	awk '$1>0 {print $2;}' 2-3.60.vcf.count > 2-3.60.vcf.count.1
	cd ../..
	pwd
	
done


time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd	
#	awk '$1>0 {print $2;}' 2-3.60.vcf.count > 2-3.60.vcf.count.1
	for f in $(cat 2-3.60.vcf.count.1); do echo $f; cat $f; done > 2-3.60.vcf.count.1.input
	#grep -B1 "34927253" 60.vcf.count.1.input
	cut -f2 call.2-3.60.vcf | sort | uniq -c | sort -k1nr 
	cd ../..
	pwd
	
done
conda deactivate 
#################  3up 
#time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; awk -F ";|=" ' $3>1 {print $1}' test.umi12.f.60K.u.c.94.title > test.umi12.f.60K.u.c.94.3.title.ID; awk -F "\t|;" 'FNR==NR{a[$0];next}($9 in a)' test.umi12.f.60K.u.c.94.3.title.ID test.umi12.f.60K.u.c.94.txt | cut -f2 | sort -k1n | uniq > test.umi12.f.60K.u.c.94.3.clusters; wc -l test.umi12.f.60K.u.c.94.3.title.ID; wc -l test.umi12.f.60K.u.c.94.3.clusters; sed 's/^/test.umi12.f.60K.u.c.94.4up.clusters/g' test.umi12.f.60K.u.c.94.3.clusters > test.umi12.f.60K.u.c.94.3.clusters.rename; cd ..; pwd; done


cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101
pwd

echo " 90 Percentile "
conda activate longread_umi_Intel



time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; usearch -cluster_fast test.umi12.f.60K.u.fasta -id 0.94 -centroids test.umi12.f.60K.u.c.94.3up.fasta -uc test.umi12.f.60K.u.c.94.3up.txt -sizein -sizeout  -strand both -minsize 3 -sort size ; grep ">" test.umi12.f.60K.u.c.94.3up.fasta | sed 's/^>//g' > test.umi12.f.60K.u.c.94.3up.title; cut -d ";" -f1 test.umi12.f.60K.u.c.94.3up.title > test.umi12.f.60K.u.c.94.3up.title.ID; cut -d ";" -f2 test.umi12.f.60K.u.c.94.3up.title | sort | uniq -c > test.umi12.f.60K.u.c.94.3up.title.count; awk -F "\t|;" 'FNR==NR{a[$0];next}($9 in a)' test.umi12.f.60K.u.c.94.3up.title.ID test.umi12.f.60K.u.c.94.3up.txt | cut -f2 | sort -k1n | uniq > test.umi12.f.60K.u.c.94.3up.clusters; wc -l test.umi12.f.60K.u.c.94.3up.title.ID; wc -l test.umi12.f.60K.u.c.94.3up.clusters; sed 's/^/test.umi12.f.60K.u.c.94.4up.clusters/g' test.umi12.f.60K.u.c.94.3up.clusters > test.umi12.f.60K.u.c.94.3up.clusters.rename; cd ..; pwd; done


   
time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	pwd 
	
	sed 's/^/test.umi12.f.60K.u.c.94.4up.clusters/g' test.umi12.f.60K.u.c.94.3up.clusters > test.umi12.f.60K.u.c.94.3up.clusters.rename;
	
	cd ..
done
	
time for f in $(cat nano.input.file); do  cd ${f%.fastq.gz}; wc -l test.umi12.f.60K.u.c.94.3up.clusters; cd ..; done

conda deactivate

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd
	
	rm 3up.50.vcf.count
	rm 3up.call.50.vcf
	touch 3up.50.vcf.count
	touch 3up.call.50.vcf
	
	for f in $(cat ../test.umi12.f.60K.u.c.94.3up.clusters.rename); do wc -l  $f.50.vcf >> 3up.50.vcf.count; cat $f.50.vcf >> 3up.call.50.vcf; done
	awk '$1>0 {print $2;}' 3up.50.vcf.count > 3up.50.vcf.count.1
	
	cd ../..
	
	pwd
	
done

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd	
	echo "3up output"
	#awk '$1>0 {print $2;}' 50.vcf.count > 50.vcf.count.1
	for f in $(cat 3up.50.vcf.count.1); do echo $f; cat $f; done > 3up.50.vcf.count.1.input
	#grep -B1 "53428" 50.vcf.count.1.input
	#grep -B1 "534289" 3up.50.vcf.count.1.input
	cut -f2 3up.call.50.vcf | sort | uniq -c | sort -k1nr 
	cd ../..
	pwd
	
done

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd
	
	rm 3up.60.vcf.count
	rm 3up.call.60.vcf
	touch 3up.60.vcf.count
	touch 3up.call.60.vcf
	
	for f in $(cat ../test.umi12.f.60K.u.c.94.3up.clusters.rename); do wc -l  $f.60.vcf >> 3up.60.vcf.count; cat $f.60.vcf >> 3up.call.60.vcf; done
	awk '$1>0 {print $2;}' 3up.60.vcf.count > 3up.60.vcf.count.1
	
	cd ../..
	
	pwd
	
done

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd	
	echo "3up output"
	#awk '$1>0 {print $2;}' 60.vcf.count > 60.vcf.count.1
	for f in $(cat 3up.60.vcf.count.1); do echo $f; cat $f; done > 3up.60.vcf.count.1.input
	#grep -B1 "53428" 60.vcf.count.1.input
#	grep -B1 "302" 3up.60.vcf.count.1.input
	cut -f2 3up.call.60.vcf | sort | uniq -c | sort -k1nr -k2nr
	cd ../..
	pwd
	
done


rm $Run_name-13V-run1.60K.3up.count
touch $Run_name-13V-run1.60K.3up.count

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd	
	echo "3up output"
	#awk '$1>0 {print $2;}' 60.vcf.count > 60.vcf.count.1
#	for f in $(cat 3up.60.vcf.count.1); do echo $f; cat $f; done > 3up.60.vcf.count.1.input
	#grep -B1 "53428" 60.vcf.count.1.input
#	head -20 3up.60.vcf.count.1.input
#	head -20 3up.60.vcf.count
#	head -20 3up.call.60.vcf
	
	cut -f2 3up.call.60.vcf | sort | uniq -c | sort -k1nr -k2nr >> ../../$Run_name-13V-run1.60K.3up.count
	cd ../..
	pwd
	
done

####### true ave coverage by calculation 
############## cutoff at 94% avg-1 size ################

####### assuming the 4up is done ###################################
######### change the -minsize to the average coverage ###########

cd $PATH_to_basecalled_files/$Run_name-13V-run1/duplex_101
pwd

echo " 94% avg-cal size "
conda activate longread_umi_Intel

time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; ave_cov=$(grep -o 'size=[0-9]*' test.umi12.f.60K.u.c.fasta | awk -F= '{sum+=$2; n++} END {if (n>0) print int((sum/n)+0.999999); else print 0}'); echo "Average-cov: $ave_cov"; usearch -cluster_fast test.umi12.f.60K.u.fasta -id 0.94 -centroids test.umi12.f.60K.u.c.94.Ave-cal.fasta -uc test.umi12.f.60K.u.c.94.Ave-cal.txt -sizein -sizeout  -strand both -minsize "$ave_cov" -sort size ; grep ">" test.umi12.f.60K.u.c.94.Ave-cal.fasta | sed 's/^>//g' > test.umi12.f.60K.u.c.94.Ave-cal.title; cut -d ";" -f1 test.umi12.f.60K.u.c.94.Ave-cal.title > test.umi12.f.60K.u.c.94.Ave-cal.title.ID; cut -d ";" -f2 test.umi12.f.60K.u.c.94.Ave-cal.title | sort | uniq -c > test.umi12.f.60K.u.c.94.Ave-cal.title.count; awk -F "\t|;" 'FNR==NR{a[$0];next}($9 in a)' test.umi12.f.60K.u.c.94.Ave-cal.title.ID test.umi12.f.60K.u.c.94.Ave-cal.txt | cut -f2 | sort -k1n | uniq > test.umi12.f.60K.u.c.94.Ave-cal.clusters; wc -l test.umi12.f.60K.u.c.94.Ave-cal.title.ID; wc -l test.umi12.f.60K.u.c.94.Ave-cal.clusters; sed 's/^/test.umi12.f.60K.u.c.94.4up.clusters/g' test.umi12.f.60K.u.c.94.Ave-cal.clusters > test.umi12.f.60K.u.c.94.Ave-cal.clusters.rename; cd ..; pwd; done


   
time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	pwd 
	
	sed 's/^/test.umi12.f.60K.u.c.94.4up.clusters/g' test.umi12.f.60K.u.c.94.Ave-cal.clusters > test.umi12.f.60K.u.c.94.Ave-cal.clusters.rename;
	
	cd ..
done
	
conda deactivate

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd
	
	rm Ave-cal.50.vcf.count
	rm Ave-cal.call.50.vcf
	touch Ave-cal.50.vcf.count
	touch Ave-cal.call.50.vcf
	
	for f in $(cat ../test.umi12.f.60K.u.c.94.Ave-cal.clusters.rename); do wc -l  $f.50.vcf >> Ave-cal.50.vcf.count; cat $f.50.vcf >> Ave-cal.call.50.vcf; done
	awk '$1>0 {print $2;}' Ave-cal.50.vcf.count > Ave-cal.50.vcf.count.1
	
	cd ../..
	
	pwd
	
done

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd	
	echo "Ave-cal output"
	#awk '$1>0 {print $2;}' 50.vcf.count > 50.vcf.count.1
	for f in $(cat Ave-cal.50.vcf.count.1); do echo $f; cat $f; done > Ave-cal.50.vcf.count.1.input
	#grep -B1 "53428" 50.vcf.count.1.input
	#grep -B1 "73" Ave-cal.50.vcf.count.1.input
	cut -f2 Ave-cal.call.50.vcf | sort | uniq -c | sort -k1nr 
	cd ../..
	pwd
	
done

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd
	
	rm Ave-cal.60.vcf.count
	rm Ave-cal.call.60.vcf
	touch Ave-cal.60.vcf.count
	touch Ave-cal.call.60.vcf
	
	for f in $(cat ../test.umi12.f.60K.u.c.94.Ave-cal.clusters.rename); do wc -l  $f.60.vcf >> Ave-cal.60.vcf.count; cat $f.60.vcf >> Ave-cal.call.60.vcf; done
	awk '$1>0 {print $2;}' Ave-cal.60.vcf.count > Ave-cal.60.vcf.count.1
	
	cd ../..
	
	pwd
	
done

time for f in $(cat nano.input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd	
	echo "Ave-cal output"
	#awk '$1>0 {print $2;}' 60.vcf.count > 60.vcf.count.1
	for f in $(cat Ave-cal.60.vcf.count.1); do echo $f; cat $f; done > Ave-cal.60.vcf.count.1.input
	#grep -B1 "53428" 60.vcf.count.1.input
#	grep -B1 "302" Ave-cal.60.vcf.count.1.input
	cut -f2 Ave-cal.call.60.vcf | sort | uniq -c | sort -k1nr -k2nr
	
	cd ../..
	pwd
	
done

rm $Run_name-13V-run1.60K.Ave-cal.count
touch $Run_name-13V-run1.60K.Ave-cal.count

time for f in $(cat nano.input.file);
do
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_60K_4up;
	pwd	
	echo "Ave-cal output"

	
	cut -f2 Ave-cal.call.60.vcf | sort | uniq -c | sort -k1nr -k2nr >> ../../$Run_name-13V-run1.60K.Ave-cal.count;
	
	cd ../..
	pwd
	
done


time for f in $(cat nano.input.file); do echo $f;  cd ${f%.fastq.gz};  cd cluster_94_60K_4up; pwd; echo "Ave-cal output";  cut -f2 Ave-cal.call.60.vcf | sort | uniq -c | sort -k1nr -k2nr >> ../../$Run_name-13V-run1.60K.Ave-cal.count; cd ../..; done

time for f in $(cat nano.input.file); do echo $f; cd ${f%.fastq.gz}; pwd; gzip *.fastq; cd ..; done


date

echo "Part 2 finished, starting part 3:"
cd $PATH_to_script_files
pwd

time sh $Run_name-13V-run1.101.dup.part3.60K.sh
































