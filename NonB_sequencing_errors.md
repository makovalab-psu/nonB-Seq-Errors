# Non-B DNA associated errors

## General notes

This is a collection of scripts that has been used to curate data, perform analysis, and generate figures for Weissensteiner et al. 2022. The overall workflow is as follows: After mapping sequencing read data (in case of HiFi and ONT) to the human reference genome hg19, or retrieving alignment (bam) files from the Genome in a Bottle website ([giab_data_indexes/AshkenazimTrio at master Â· genome-in-a-bottle/giab_data_indexes Â· GitHub](https://github.com/genome-in-a-bottle/giab_data_indexes/tree/master/AshkenazimTrio)). Then non-B-DNA forming motifs are either retrieved from downloaded from [Advanced Biomedical Computing Center (ABCC) | non-B DB | Home](https://nonb-abcc.ncifcrf.gov/apps/site/default) or, in case of G4 motifs, generated. Following that motifs are curated, and sequence-specific features (nucleotide composition, mappability, etc.) are retrieved for each one, and a random set of controls is generated that matches motifs in size and number. Then, measures of sequencing success (single-nucleotide, insertion, and deletion mismatches, read depth and base quality), are determined for every motif and control. For every motif and control sequence, we assigned a unique ID, which consists of the non-B type, source (motif / control), chromosome, and start of the bed entry (e.g., G4Motifs_motif_1_123456)). In the downstream analysis, we used Poisson regression models to determine the effect of non-B motifs on the variability of sequencing error rates

## Read mapping

After downloading fastq files of PacBio HiFi and ONT, sequencing reads were aligned to hg19 using *minimap2*

```bash
#!/bin/bash                                                                                                                                                                   
#SBATCH --job-name=matthias_map_ont_reads                                                                                                                                     
#SBATCH --output=matthias-%j.out                                                                                                                                              
#SBATCH --error=matthias-%j.err                                                                                                                                               
#SBATCH --mem-per-cpu=16G                                                                                                                                                     
#SBATCH -n 8                                                                                                                                                                  


reference="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/ONT_error_detection/reference/hg19_formated_by_wil.fa"                                                                    
file=$1                                                                                                                                                                       
data=$2                                                                                                                                                                       
output_folder="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/bam_files/${data}"                                                                                                    
id=$(echo $file | tr "/" "\t" | cut -f 9,12 | awk '{print $1 "." $2}' )                                                                                                       
echo $id                                                                                                                                                                      
/galaxy/home/mxw1081/software/minimap2/minimap2 -ax map-ont  -t 8 $reference $file > ${output_folder}/${id}.minimap2.sam                                                      

awk '$3 != "*" || $1 ~ /^@/' ${output_folder}/${id}.minimap2.sam | \                                                                                                          
samtools view -Sb - |\                                                                                                                                                        
samtools sort - -o ${output_folder}/${id}.minimap2.sorted.bam                                                                                                                 

samtools index ${output_folder}/${id}.minimap2.sorted.bam                                                                                                                     

rm ${output_folder}/${id}.minimap2.sam 
```

## BAM file curation

Using the following scripts, we split the bam files by forward and reverse aligning reads, and by chromosome. 

```bash
#!/bin/bash                                                                                                                                                                   
#SBATCH --job-name=matthias_map_ont_reads                                                                                                                                     
#SBATCH --output=matthias-%j.out                                                                                                                                              
#SBATCH --error=matthias-%j.err                                                                                                                                               
#SBATCH --mem-per-cpu=16G                                                                                                                                                     
#SBATCH -n 8                                                                                                                                                                  

dir="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/bam_files/ONT"                                                                                                                  


samtools view -b -F 16 ${dir}/ONT_HG0002_20200914_bonito_v0.3.0.merged.bam | \                                                                                                
   samtools sort - > ${dir}/ONT_HG0002_20200914_bonito_v0.3.0_FWD.sorted.bam                                                                                                  
samtools index ${dir}/ONT_HG0002_20200914_bonito_v0.3.0_FWD.sorted.bam                                                                                                        

samtools view -b -f 16 ${dir}/ONT_HG0002_20200914_bonito_v0.3.0.merged.bam | \                                                                                                
   samtools sort - > ${dir}/ONT_HG0002_20200914_bonito_v0.3.0_REV.sorted.bam                                                                                                  
samtools index ${dir}/ONT_HG0002_20200914_bonito_v0.3.0_REV.sorted.bam                                                                                                        


array=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)                                                                                                              

for chromosome in "${array[@]}"; do                                                                                                                                           
echo $chromosome                                                                                                                                                              
samtools view -b -F 16 ${dir}/HG002.hs37d5.100x.sorted.rmdup_FWD.sorted.bam "$chromosome" >  \                                                                                
 ${dir}/fwd/${chromosome}_FWD.bam                                                                                                                                             
samtools index ${dir}/fwd/${chromosome}_FWD.bam                                                                                                                               

samtools view -b -f 16 ${dir}/HG002.hs37d5.100x.sorted.rmdup_REV.sorted.bam "$chromosome" >  \                                                                                
 ${dir}/rev/${chromosome}_REV.bam                                                                                                                                             
samtools index ${dir}/rev/${chromosome}_REV.bam                                                                                                                               

done
```

## Motif curation

### Filter motifs and generate controls

Motif bed files (downloaded from [Advanced Biomedical Computing Center (ABCC) | non-B DB | Home](https://nonb-abcc.ncifcrf.gov/apps/site/default)) are further filtered and controls are generated on the basis of this.

```bash
#!/bin/bash                                                                                                                                                                                                                                                                                              
#SBATCH --job-name=matthias_map_ont_reads                                                                                                                                                                                                                                                                
#SBATCH --output=matthias-%j.out                                                                                                                                                                                                                                                                         
#SBATCH --error=matthias-%j.err                                                                                                                                                                                                                                                                          
#SBATCH --mem-per-cpu=16G                                                                                                                                                                                                                                                                                
#SBATCH -n 4                                                                                                                                                                                                                                                                                             

input_folder="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/nonB_motifs/new_motifs/new_runs"                                                                                                                                                                                                                  
output_folder="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/nonB_motifs/homo_polymers_excluded"                                                                                                                                                                                                              

#motif_list=$1                                                                                                                                                                                                                                                                                           

homopolymer_file="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/references/hg19/homopolymers_7_min_hg19_curated.bed"                                                                                                                                                                                          
gap_file="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/references/hg19/hg19.gaps.bed"                                                                                                                                                                                                                        
genome_file="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/references/hg19/hg19_edited_genome_file.genome"                                                                                                                                                                                                    

#while IFS= read -r file; do                                                                                                                                                                                                                                                                             
#echo $file                                                                                                                                                                                                                                                                                              
#id=$(echo $file | cut -f 1 -d '.')                                                                                                                                                                                                                                                                      
#/galaxy/home/mxw1081/software/bedtools2/bin/bedtools intersect -v -a ${input_folder}/$file -b $homopolymer_file > ${output_folder}/${id}.bed                                                                                                                                                            
#done<$motif_list                                                                                                                                                                                                                                                                                        

/galaxy/home/mxw1081/software/bedtools2/bin/bedtools shuffle -chrom \                                                                                                                                                                                                                                    
   -noOverlapping \                                                                                                                                                                                                                                                                                      
   -maxTries 1000000 \                                                                                                                                                                                                                                                                                   
   -g $genome_file \                                                                                                                                                                                                                                                                                     
   -i  ${output_folder}/APhasedRepeats.bed  \                                                                                                                                                                                                                                                            
   -excl <(cat $homopolymer_file $gap_file ${output_folder}/DirectRepeats.bed ${output_folder}/InvertedRepeats.bed ${output_folder}/MirrorRepeats.bed ${output_folder}/ZDNAMotifs.bed ${output_folder}/G4Motifs.bed) | \                                                                                 
   /galaxy/home/mxw1081/software/bedtools2/bin/bedtools sort -i -   >  ${output_folder}/APhasedRepeats_controls.bed                                                                                                                                                                                      

/galaxy/home/mxw1081/software/bedtools2/bin/bedtools shuffle -chrom \                                                                                                                                                                                                                                    
   -noOverlapping \                                                                                                                                                                                                                                                                                      
   -maxTries 1000000 \                                                                                                                                                                                                                                                                                   
   -g $genome_file \                                                                                                                                                                                                                                                                                     
   -i  ${output_folder}/DirectRepeats.bed  \                                                                                                                                                                                                                                                             
   -excl <(cat $homopolymer_file $gap_file ${output_folder}/APhasedRepeats.bed ${output_folder}/InvertedRepeats.bed ${output_folder}/MirrorRepeats.bed ${output_folder}/ZDNAMotifs.bed ${output_folder}/G4Motifs.bed) | \                                                                                
   /galaxy/home/mxw1081/software/bedtools2/bin/bedtools sort -i -   >  ${output_folder}/DirectRepeats_controls.bed                                                                                                                                                                                       

/galaxy/home/mxw1081/software/bedtools2/bin/bedtools shuffle -chrom \                                                                                                                                                                                                                                    
   -noOverlapping \                                                                                                                                                                                                                                                                                      
   -maxTries 1000000 \                                                                                                                                                                                                                                                                                   
   -g $genome_file \                                                                                                                                                                                                                                                                                     
   -i  ${output_folder}/InvertedRepeats.bed  \                                                                                                                                                                                                                                                           
   -excl <(cat $homopolymer_file $gap_file ${output_folder}/APhasedRepeats.bed ${output_folder}/DirectRepeats.bed ${output_folder}/MirrorRepeats.bed ${output_folder}/ZDNAMotifs.bed ${output_folder}/G4Motifs.bed) | \                                                                                  
   /galaxy/home/mxw1081/software/bedtools2/bin/bedtools sort -i -   >  ${output_folder}/InvertedRepeats_controls.bed                                                                                                                                                                                     

/galaxy/home/mxw1081/software/bedtools2/bin/bedtools shuffle -chrom \                                                                                                                                                                                                                                    
   -noOverlapping \                                                                                                                                                                                                                                                                                      
   -maxTries 1000000 \                                                                                                                                                                                                                                                                                   
   -g $genome_file \                                                                                                                                                                                                                                                                                     
   -i  ${output_folder}/MirrorRepeats.bed  \                                                                                                                                                                                                                                                             
   -excl <(cat $homopolymer_file $gap_file ${output_folder}/APhasedRepeats.bed ${output_folder}/DirectRepeats.bed ${output_folder}/InvertedRepeats.bed ${output_folder}/ZDNAMotifs.bed ${output_folder}/G4Motifs.bed) | \                                                                                
   /galaxy/home/mxw1081/software/bedtools2/bin/bedtools sort -i -   >  ${output_folder}/MirrorRepeats_controls.bed                                                                                                                                                                                       

/galaxy/home/mxw1081/software/bedtools2/bin/bedtools shuffle -chrom \                                                                                                                                                                                                                                    
   -noOverlapping \                                                                                                                                                                                                                                                                                      
   -maxTries 1000000 \                                                                                                                                                                                                                                                                                   
   -g $genome_file \                                                                                                                                                                                                                                                                                     
   -i  ${output_folder}/G4Motifs.bed  \                                                                                                                                                                                                                                                                  
   -excl <(cat $homopolymer_file $gap_file ${output_folder}/APhasedRepeats.bed ${output_folder}/DirectRepeats.bed ${output_folder}/InvertedRepeats.bed ${output_folder}/ZDNAMotifs.bed ${output_folder}/MirrorRepeats.bed) | \                                                                           
   /galaxy/home/mxw1081/software/bedtools2/bin/bedtools sort -i -   >  ${output_folder}/G4Motifs_controls.bed                                                                                                                                                                                            

/galaxy/home/mxw1081/software/bedtools2/bin/bedtools shuffle -chrom \                                                                                                                                                                                                                                    
   -noOverlapping \                                                                                                                                                                                                                                                                                      
   -maxTries 1000000 \                                                                                                                                                                                                                                                                                   
   -g $genome_file \                                                                                                                                                                                                                                                                                     
   -i  ${output_folder}/ZDNAMotifs.bed  \                                                                                                                                                                                                                                                                
   -excl <(cat $homopolymer_file $gap_file ${output_folder}/APhasedRepeats.bed ${output_folder}/DirectRepeats.bed ${output_folder}/InvertedRepeats.bed ${output_folder}/G4Motifs.bed ${output_folder}/MirrorRepeats.bed) | \                                                                             
   /galaxy/home/mxw1081/software/bedtools2/bin/bedtools sort -i -   >  ${output_folder}/ZDNAMotifs_controls.bed 
```

### G4 motif prediction with Quadron

To computationally predict G4 motifs, we used *quadron* with default parameters. 

```r
################################################################################                                                                                              
# Requires the libraries "doMC", "foreach", "itertools", "xgboost" (0.4-4),    #                                                                                              
# "caret" and "plyr".                                                          #                                                                                              
# If not already installed in R, you can install those by typing:              #                                                                                              
# install.packages(c("doMC", "foreach", "itertools", "plyr", "caret"))         #                                                                                              
# Specific steps are needed to install the xgboost version 0.4-4, as detailed  #                                                                                              
# in the Readme file.                                                          #                                                                                              
# The default fastread==TRUE option in readfasta requires "data.table" library.#                                                                                              
################################################################################                                                                                              
#setwd("./lib")                                                                                                                                                               
#source("bitcompile.R")                                                                                                                                                       
#rm(list=ls())                                                                                                                                                                
#setwd("..")                                                                                                                                                                  
.libPaths(c(.libPaths(), "/opt/bx/R/xgboost_0.4-4"))                                                                                                                          

print("NOTE: Loading Quadron core...", quote=FALSE)                                                                                                                           
load("/galaxy/home/mxw1081/software/Quadron-master/Quadron.lib")                                                                                                              


args = commandArgs(trailingOnly=TRUE)                                                                                                                                         
Quadron(FastaFile= paste0("/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/references/hg19/single_chromosomes/", args[1], ".fa"),                                                
        OutFile  = paste0("/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/scripts/quadron/hg19_output/",args[1], ".txt"),                                                       
        nCPU     = 8,                                                                                                                                                         
        SeqPartitionBy = 1000000)                                                                                                                                             

#file.remove("Quadron.lib")                                                                                                                                                   
#rm(list=ls())                                                                                                                                                                
################################################################################
```

### Remove overlaps within motifs

With this script, we removed overlapping intervals within motif bed file. 

```bash
#!/bin/bash                                                                                                                                                                   
#SBATCH --job-name=matthias_map_ont_reads                                                                                                                                     
#SBATCH --output=matthias-%j.out                                                                                                                                              
#SBATCH --error=matthias-%j.err                                                                                                                                               
#SBATCH --mem-per-cpu=16G                                                                                                                                                     
#SBATCH -n 4                                                                                                                                                                  

input_folder="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/nonB_motifs/homo_polymers_excluded"                                                                                    
output_folder="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/nonB_motifs/overlaps_removed"  while IFS= read -r file; do                                                                                                                                                   
id=$(echo $file | cut -f 1 -d '.')                                                                                                                                            
echo $id                                                                                                                                                                      
/galaxy/home/mxw1081/software/bedtools2/bin/bedtools merge -i ${input_folder}/${id}_controls.bed -c 1 -o count > counted                                                      
awk '/\t1$/{print}' counted > filtered                                                                                                                                        
/galaxy/home/mxw1081/software/bedtools2/bin/bedtools intersect -a ${input_folder}/${id}_controls.bed -b filtered -wa | \                                                      
 /galaxy/home/mxw1081/software/bedtools2/bin/bedtools sort -i -   > ${output_folder}/${id}_controls.bed                                                                       

rm counted                                                                                                                                                                    
rm filtered                                                                                                                                                                   

entry_number=$(wc -l ${output_folder}/${id}_controls.bed | awk '{print $1}')                                                                                                  
echo $entry_number                                                                                                                                                            

/galaxy/home/mxw1081/software/bedtools2/bin/bedtools sample -n $entry_number -i ${input_folder}/${id}.bed | \                                                                 
   /galaxy/home/mxw1081/software/bedtools2/bin/bedtools sort -i -   > ${output_folder}/${id}.bed                                                                              
done<$motif_list 
```

## Per motif features

### Mappability

To calculate the average per-motif mappability scores,  we used this R script:

```r
require(GenomicRanges)
require(rtracklayer)
require(parallel)
#library(easypackages)
#libraries("ggplot2", "car", "scales", "plyr", "tidyr","msme", "dplyr", "glmmADMB", "purrr", "grid", "gridExtra", "lme4", "MatchIt", "compositions")
setwd("/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/scripts/mappability")

args = commandArgs(trailingOnly=TRUE)
input_bed<-args[1]

  windows=import(paste0(input_bed))
core_number <- detectCores()
worker_number <- core_number-1
cl <- makeCluster(worker_number,timeout=60*60*24*30) # 30 days time limit

n_windows_group <- floor(length(windows)/worker_number)
groups <- c(rep.int(1,length(windows)-worker_number*n_windows_group),
            rep(seq_len(worker_number),each=n_windows_group))
windows_groups <- split(windows,groups)

df <- data.frame(as.character(seqnames(windows)), start(windows)-1, end(windows))
colnames <- c('chrom', 'start', 'end')

colnames <- c(colnames, 'mappability')
  name<-"test"
  start=proc.time()
  signal_1Mb=Reduce(c,parLapplyLB(cl,windows_groups,
  function(windows,name){
   require(rtracklayer)
   count=import(paste0('wgEncodeCrgMapabilityAlign100mer.bigWig'),format="bw", which=windows)
   overlaps=findOverlaps(windows, count)
   overlaps_windows=queryHits(overlaps) 
   overlaps_count=subjectHits(overlaps)
   count_splitted=vector("list",length(windows))
   count_splitted[unique(overlaps_windows)]=as.list(split(count[overlaps_count],overlaps_windows))

   return(
     mapply(function(count,i){
       if(length(count)==0)
         return(0)
       y=windows[i]
       if(start(count[1])<start(y))
         start(count[1])=start(y)
       if(end(tail(count,1))>end(y))
         end(count[length(count)])=end(y)
       return(sum(count$score*width(count))/width(y))
     },count_splitted,seq_along(count_splitted))
   )
 },name=name))
  end=proc.time()
  end-start

  # SAVE  signal_1Mb somewhere
  df <- cbind(df, as.vector(unlist(signal_1Mb)))    
names(df) <- c("chr", "start", "end", "mappability")
write.table(df, file = paste0(input_bed, ".mappability"), sep = "\t", 
            row.names = FALSE, col.names = FALSE,  quote = FALSE)
stopCluster(cl)
```

In combination with this wrapper to run it on the cluster:

```bash
#!/bin/bash
#SBATCH --job-name=matthias_map_ont_reads
#SBATCH --output=matthias-%j.out 
#SBATCH --error=matthias-%j.err 
#SBATCH --mem-per-cpu=16G 
#SBATCH -n 16


motif_folder="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/nonB_motifs/RM_SD_STR_filter"

file=$1

echo ${motif_folder}/$file
awk '{print "chr" $0}'  ${motif_folder}/$file > ${file}.temp.bed 
Rscript /nfs/brubeck.bx.psu.edu/scratch6/mxw1081/scripts/mappability/mappability_on_cluster.R ${file}.temp.bed 
rm *.temp.bed
```

And to parse the output (i.e.) get a mappability value for each unique ID we use: 

```bash
#!/bin/bash

list=$1

while IFS= read -r line; do
    source=$(echo $line | awk '{print $2}') 
    type=$(echo $line | cut -f 1 -d '.'  )
    file=$(echo $line | awk '{print $1}')
    sed 's/chr//g' $file | awk -v s=$source -v t=$type '{print t "_" s "_" $1 "_" $2 "\t" $4}'
    #file=$(echo $line | cut -f 1-3 -d '.' | awk '{print $1 ".output"}')
    #awk -v s=$source -v t=$type '{print t "_" s "_" $1 "_" $2 "\t" $0 "\t" $4/$7 "\t" $5/$7 "\t" $6/$7}' ${results_folder}/$data/$orientation/$file
done<$list
```

### Base quality

To calculate the average base quality for each motif we use the following scripts:

The master script to submit many jobs in parallel. Note that this also needs the conversion tables as input (one for Illumina and one for HiFi, since those use different base quality scoring):

```bash
#!/bin/bash
#SBATCH --job-name=matthias_error_pipeline_master
#SBATCH --output=matthias-%j.out 
#SBATCH --error=matthias-%j.err 
#SBATCH --mem-per-cpu=16G 
#SBATCH -n 1

list=$1
data=$2
orientation=$3
script_folder="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/scripts/basequality_estimation"

job_limit=1000
long_read="HiFi"

if [ "$data" = "$long_read" ]; then
   echo "Data is long read. Using long read table."
   conversion_table="${script_folder}/long_read_quality_conversion_table.txt"
else
   echo "Data is short read. Using short read conversion table."
   conversion_table="${script_folder}/short_read_quality_conversion_table.txt"
fi

echo $conversion_table

while IFS= read -r file; do
echo $file
job_name=$(echo $file | cut -f 1-3 -d '.') 
#node=$(/galaxy/home/mxw1081/software/nick-bfx/slurm-wait.py --min-idle-nodes 1 --mem 16G )
#sbatch -w $node -J $job_name ${script_folder}/run_bobs_error_caller.sh $data $orientation $file
   nr_jobs=$(squeue | grep 'mxw' | wc -l)
   echo $nr_jobs
   while [ "$nr_jobs" =  "$job_limit" ]; do
      echo 'job submission limit reached.'
      nr_jobs=$(squeue | grep 'mxw' | wc -l)
      sleep 10
   done   
sbatch  --exclude nn[4-5] --nice -J $job_name ${script_folder}/get_base_q.sh $data $orientation  $conversion_table $file
sleep 0.1
done<$list
```

In the master script, the `get_base_q.sh`  script gets called, which uses samtools and a python script (see below)

```bash
#!/bin/bash
#SBATCH --job-name=matthias_mapping_quality
#SBATCH --output=matthias-%j.out 
#SBATCH --error=matthias-%j.err 
#SBATCH --mem-per-cpu=16G 
#SBATCH -n 1
#SBATCH --nice=100
## Illumina data:

data=$1
orientation=$2
conversion_table=$3

bamfolder="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/bam_files/${data}/${orientation}"
inputfolder="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/scripts/error_detection_pipeline/input"
basefolder="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/scripts/basequality_estimation"
outputfolder="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/scripts/basequality_estimation/output/${data}/${orientation}"
input_bed=$4

id=$(echo $input_bed | cut -f 1-3 -d '.')
echo $id

while IFS= read -r line; do
#echo $line

chr=$(echo $line | awk '{print $1}' | sed 's/chr//g')
from=$(echo $line | awk '{print $2}')
to=$(echo $line | awk '{print $3}')

average_base_q=$(samtools mpileup ${bamfolder}/${chr}.bam --region "${chr}:${from}-${to}" | \
 python ${basefolder}/get_base_quality.py ${conversion_table} - | \
awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }')
echo $line $average_base_q | awk '{if(NF != 4) print $1 "\t" $2 "\t" $3 "\t" "NA"; else print $0;}'  >>${outputfolder}/${id}.basequality.output

echo $line
#samtools mpileup ${bamfolder}/${chr}.bam --region "${chr}:${from}-${to}" | \ 
# python ${basefolder}/get_base_quality.py ${basefolder}/conversion_table  - #| \
#awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }')
#echo $line $average_base_q | awk '{if(NF != 4) print $1 "\t" $2 "\t" $3 "\t" "NA"; else print $0;}'  >>${basefolder}/${type}/${input}.output

done<${inputfolder}/$input_bed
```

Here is the python script `get_base_quality.py`

```python
#!/bin/python

import sys
import csv

ref = open(sys.argv[1], 'rt')

dict = {}
for line in ref:
    array = line.strip().split('\t')
    #key, value = array
    dict[array[1]] = array[0]

#file = read(sys.stdin[2], 'rt')
for line in sys.stdin:
    x = list(line.strip().split('\t')[5])
    per_line = []
    for item in x:
        per_line.append(float(dict[item]))
print(sum(per_line)/len(per_line))
```

And here are the conversion tables:

`short_read_quality_conversion_table.txt`

```txt
0    !
1    "
2    #
3    $
4    %
5    &
6    '
7    (
8    )
9    *
10    +
11    ,
12    -
13    .
14    /
15    0
16    1
17    2
18    3
19    4
20    5
21    6
22    7
23    8
24    9
25    :
26    ;
27    <
28    =
29    >
30    ?
31    @
32    A
33    B
34    C
35    D
36    E
37    F
38    G
39    H
40    I
41    J
```

`long_read_quality_conversion_table.txt`

```
0    !
1    "
2    #
3    $
4    %
5    &
6    '
7    (
8    )
9    *
10    +
11    ,
12    -
13    .
14    /
15    0
16    1
17    2
18    3
19    4
20    5
21    6
22    7
23    8
24    9
25    :
26    ;
27    <
28    =
29    >
30    ?
31    @
32    A
33    B
34    C
35    D
36    E
37    F
38    G
39    H
40    I
41    J
42    K
43    L
44    M
45    N
46    O
47    P
48    Q
49    R
50    S
51    T
52    U
53    V
54    W
55    X
56    Y
57    Z
58    [
59    \
60    ]
61    ^
62    _
63    `
64    a
65    b
66    c
67    d
68    e
69    f
70    g
71    h
72    i
73    j
74    k
75    l
76    m
77    n
78    o
79    p
80    q
81    r
82    s
83    t
84    u
85    v
86    w
87    x
88    y
89    z
90    {
91    |
92    }
93    ~
```

In the `parse_output.sh` script, each motif and control gets assigned the uniqe ID and the associated base quality value (to be able to merge results later):

```bash
#!/bin/bash

list=$1
data=$2
orientation=$3
results_folder="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/scripts/basequality_estimation/output"

while IFS= read -r line; do
    #echo $line
    source=$(echo $line | awk '{print $2}') 
    type=$(echo $line | awk '{print $1}' | cut -f 2 -d '.' | cut -f 1 -d '_')
    #echo $type $source
    file=$(echo $line | cut -f 1-3 -d '.' | awk '{print $1 ".basequality.output"}')
    awk -v s=$source -v t=$type '{print t "_" s "_" $1 "_" $2 "\t" $4}' ${results_folder}/$data/$orientation/$file
done<$list
```

### Mismatch errors

To get all mismatch errors for each motif and control regions, as well as the total read depth, we employed a python script that naively calls mismatch errors from sam files, while considerring 'blacklisted' variants that will be excluded from the reported mismatch errors. To speed up things, the input bed files are split up in smaller files with 1000 lines each. 

This is the script to split the input bed files: `prepare_input_bed.sh`

```bash
#!/bin/bash

bed_folder="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/nonB_motifs/RM_SD_STR_filter"
output_folder="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/scripts/error_detection_pipeline/input"

array=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
#chromosome=1
for chromosome in "${array[@]}"; do
echo $chromosome

  while IFS= read -r type; do
          echo $type
          split ${bed_folder}/chr$chromosome/${type}.filtered.${chromosome}.bed -l 1000 --additional-suffix=.${type}.${chromosome}.bed
          mv x* $output_folder
  done<type_list
done
```

This is the script to run the error caller per bed file:

```bash
#!/bin/bash
#SBATCH --job-name=matthias_map_ont_reads
#SBATCH --output=matthias-%j.out 
#SBATCH --error=matthias-%j.err 
#SBATCH --mem-per-cpu=16G 
#SBATCH -n 1

data=$1
orientation=$2
input_bed=$3
output_folder="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/results/${data}/${orientation}"
input_folder="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/scripts/error_detection_pipeline/input"
bam_folder="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/bam_files/${data}/${orientation}"
reference="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/references/hg19/hg19_formated_by_wil.fa"

chr=$(cat ${input_folder}/$input_bed |   head -1 | cut -f 1)
from=$(cat ${input_folder}/$input_bed | sort -n -k 2 | head -1 | cut -f 2)
to=$(cat ${input_folder}/$input_bed |  sort -n -k 3 | tail -n 1 | cut -f 3)

echo $input_bed
echo $from 
echo $to

id=$(echo $input_bed | cut -f 1-3 -d '.')

time samtools view -S ${bam_folder}/${chr}.bam \
   "${chr}:${from}-${to}" | 
 /galaxy/home/mxw1081/software/PythoBob3/simple_error_caller.py  \
  --ref=${reference} \
  --blacklist=/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/scripts/error_detection_pipeline/blacklist_vcfs/chr${chr}.blacklist.extended.vcf \
  --blacklist:report \
  --bucketlength=1K \
  --debug=usefai \
  $input_folder/${input_bed} > $output_folder/${id}.output
```

And to submit many jobs at the same time I use `master.sh` with a list of all the input bed files

```bash
#!/bin/bash

list=$1
data=$2
orientation=$3

while IFS= read -r file; do
    echo $file
    sbatch run_bobs_error_caller.sh $data $orientation $file 
done<$list
```

### Output parsing

The output from `run_bobs_error_caller.sh`` looks like this:

```
1    991452    991479    0    0    0    962
1    1322150    1322177    2    0    0    1386
1    1387760    1387784    0    0    0    845
1    1603470    1603495    2    0    0    1011
1    1651366    1651390    17    0    0    816
1    1749335    1749360    1    0    0    1113
1    1770912    1770937    3    0    0    1579
1    1895433    1895457    3    0    0    1034
1    2043660    2043686    3    0    0    1199
1    2170556    2170582    5    0    0    1399
```

To assign a unique ID to each motif and control,  `parse_output.sh` is used:

```bash
#!/bin/bash                                                                                                                                                                                                                                                                
#SBATCH --job-name=matthias_map_ont_reads                                                                                                                                                                                                                                  
#SBATCH --output=matthias-%j.out                                                                                                                                                                                                                                           
#SBATCH --error=matthias-%j.err                                                                                                                                                                                                                                            
#SBATCH --mem-per-cpu=16G                                                                                                                                                                                                                                                  
#SBATCH -n 1 

list=$1
data=$2
orientation=$3
results_folder="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/results"

while IFS= read -r line; do

    source=$(echo $line | awk '{print $2}') 
    type=$(echo $line | awk '{print $1}' | cut -f 2 -d '.' | cut -f 1 -d '_')
    #echo $type $source
    file=$(echo $line | cut -f 1-3 -d '.' | awk '{print $1 ".output"}')
    sbatch --exclude nn[4-5] --nice -J $file --error=${file}-%j.err --output=${file}-%j.out submit_job.sh  $file $source $data $orientation 
    sleep 0.01
done<$list
```

Importantly, the input list needs to have an identifier to assign the source (motif vs. control) to each entry. This looks like this:

```
xaa.APhasedRepeats.10.bed    motif
xaa.APhasedRepeats.11.bed    motif
xaa.APhasedRepeats.12.bed    motif
xaa.APhasedRepeats.13.bed    motif
xaa.APhasedRepeats.14.bed    motif
xaa.APhasedRepeats.15.bed    motif
xaa.APhasedRepeats.16.bed    motif
xaa.APhasedRepeats.17.bed    motif
xaa.APhasedRepeats.18.bed    motif
xaa.APhasedRepeats.19.bed    motif
```

This is the `submit_job.sh` script used in the loop above:

```bash
#!/bin/bash                                                                                                                                                                                                                                                                
#SBATCH --job-name=matthias_map_ont_reads                                                                                                                                                                                                                                  
#SBATCH --output=matthias-%j.out                                                                                                                                                                                                                                           
#SBATCH --error=matthias-%j.err                                                                                                                                                                                                                                            
#SBATCH --mem-per-cpu=16G                                                                                                                                                                                                                                                  
#SBATCH -n 1 

file=$1
source=$2
data=$3
orientation=$4
results_folder="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/results"

type=$(echo $file | awk '{print $1}' | cut -f 2 -d '.' | cut -f 1 -d '_')

awk -v s=$source -v t=$type '{print t "_" s "_" $1 "_" $2 "\t" t "\t" s "\t" $0 }' ${results_folder}/${data}/${orientation}/${file} > ${results_folder}/${data}/${file}.parsed.temp
python parse_output.py ${data}_basequality_parsed mappability_output_parsed nucleotide_content_parsed overlap_with_other_motifs_IDs homopolymers_IDs ${results_folder}/${data}/${file}.parsed.temp > ${results_folder}/${data}/${file}.parsed

rm ${results_folder}/${data}/${file}.parsed.temp
```

The python script that links the output of the mismatch calling with basequality, mappability, etc, is  `parse_output.py`

```python
#!/bin/python
import pandas as pd
import re
import sys
import csv

input_basequality = open(sys.argv[1], 'rt')
file_basequality = {}
for line in input_basequality:
    array = line.strip().split('\t')
    file_basequality[array[0]] = array[1]

input_mappability = open(sys.argv[2], 'rt')
file_mappability = {}
for line in input_mappability:
    array = line.strip().split('\t')
    file_mappability[array[0]] = array[1]

input_nucleotide_content = open(sys.argv[3], 'rt')
file_nucleotide_content = {}
for line in input_nucleotide_content:
    array = line.strip().split('\t')
    file_nucleotide_content[array[0]] = array[1:7]

input_overlap = open(sys.argv[4], 'rt')
input_overlap_ids = []
for line in input_overlap:
    array = line.strip().split('\t')
    input_overlap_ids.append(array[0])

input_homopolymers = open(sys.argv[5], 'rt')
input_homopolymers_ids = []
for line in input_homopolymers:
    array = line.strip().split('\t')
    input_homopolymers_ids.append(array[0])

input_mismatches = open(sys.argv[6], 'rt')
#print("id" + "\t" + "type" + "\t" + "source" + "\t" +"chr" + "\t" +  "start" + "\t" + "end" + "\t" + "snm" + "\t" +"ins" + "\t" +"DEL" + "\t" +"depth" + "\t" + "length" + "\t" + "snm_rate" + "\t" + "ins_rate" + "\t" + "del_rate" + "\t" + "basequality" + "\t" + "mappability" + "\t" + "AT_content" + "\t" + "GC_content" + "\t" + "As" + "\t" + "Cs" + "\t" + "Gs" + "\t" + "Ts" + "\t" + "overlap_with_other_motifs" + "\t" + "homopolymers") 
for line in input_mismatches:
    array = line.strip().split('\t')
    key = array[0]
    id = array[0]
    type = array[1]
    source = array[2]
    chr = array[3]
    start = float(array[4])
    end = float(array[5])
    snm = float(array[6])
    ins = float(array[7])
    DEL = float(array[8])
    depth = float(array[9])
    length = end - start
    if depth == 0:
        depth_per_bp = "NA"
        snm_rate = "NA" 
        ins_rate = "NA"   
        del_rate = "NA"   
    else:
   snm_rate = snm / depth 
        ins_rate = ins / depth
        del_rate = DEL / depth
        depth_per_bp = depth / length
    basequality = file_basequality[key]
    mappability = file_mappability[key]
    AT_content = file_nucleotide_content[key][0]
    GC_content = file_nucleotide_content[key][1]
    As = file_nucleotide_content[key][2]
    Cs = file_nucleotide_content[key][3]
    Gs = file_nucleotide_content[key][4]
    Ts = file_nucleotide_content[key][5]
    if key in input_overlap_ids:
            overlap_with_other_motifs = "1"
    else:
            overlap_with_other_motifs = "0"
    if key in input_homopolymers_ids:
            homopolymers = "1"
    else:
            homopolymers = "0"
    print(str(id) + "\t" + str(type) + "\t" + str(source) + "\t" +str(chr) + "\t" +  str(start) + "\t" + str(end) + "\t" + str(snm) + "\t" +str(ins) + "\t" +str(DEL) + "\t" +str(depth) + "\t" + str(length) + "\t" + str(snm_rate) + "\t" + str(ins_rate) + "\t" + str(del_rate) + "\t" + str(basequality) + "\t" + str(mappability) + "\t" + str(AT_content) + "\t" + str(GC_content) + "\t" + str(As) + "\t" + str(Cs) + "\t" + str(Gs) + "\t" + str(Ts) + "\t" + str(overlap_with_other_motifs) + "\t" + str(homopolymers))
```

### Homopolymers

We first predicted the occurences of homopolymers 3-7 bp using this script (retrieved from https://www.biostars.org/p/379454/#379505) and then generated an annotation bed file using `bedtools annotate -both` 

```c
/**                                                                                                                                                                                                                                                                                     
https://www.biostars.org/p/379454/#379505                                                                                                                                                                                                                                               
Code golf: detecting homopolymers of length N in the (human) genome                                                                                                                                                                                                                     
Author: Pierre Lindenbaum                                                                                                                                                                                                                                                               
compilation:                                                                                                                                                                                                                                                                            
 gcc -O3 -Wall -o biostar379454 biostar379454.c                                                                                                                                                                                                                                         

*/                                                                                                                                                                                                                                                                                      
#include <stdio.h>                                                                                                                                                                                                                                                                      
#include <fcntl.h>                                                                                                                                                                                                                                                                      
#include <sys/types.h>                                                                                                                                                                                                                                                                  
#include <sys/stat.h>                                                                                                                                                                                                                                                                   
#include <unistd.h>                                                                                                                                                                                                                                                                     
#include <sys/io.h>                                                                                                                                                                                                                                                                     
#include <sys/mman.h>                                                                                                                                                                                                                                                                   
#include <string.h>                                                                                                                                                                                                                                                                     
#include <stdlib.h>                                                                                                                                                                                                                                                                     
#include <errno.h>                                                                                                                                                                                                                                                                      
#include <ctype.h>                                                                                                                                                                                                                                                                      

#define DUMP if(len_repeat>=len) {fputs(seq_name,stdout);printf("\t%d\t%d\t%c[%d]\n",pos-len_repeat,pos,prev_c,len_repeat); }len_repeat=0;                                                                                                                                              
#define BUF_STDOUT 1000000                                                                                                                                                                                                                                                              
int main(int argc, char const *argv[]) {                                                                                                                                                                                                                                                
    char *seq;                                                                                                                                                                                                                                                                          
    char* buff=NULL;                                                                                                                                                                                                                                                                    
    size_t size,x;                                                                                                                                                                                                                                                                      
    int len=0,prev_c=-1,len_repeat=0,pos=0;                                                                                                                                                                                                                                             
    char* seq_name=NULL;                                                                                                                                                                                                                                                                
    struct stat s;                                                                                                                                                                                                                                                                      
    int fd;                                                                                                                                                                                                                                                                             
    if(argc!=3) {                                                                                                                                                                                                                                                                       
        fprintf(stderr,"Usage: %s fasta size.\n",argv[0]);                                                                                                                                                                                                                              
        return EXIT_FAILURE;                                                                                                                                                                                                                                                            
        }                                                                                                                                                                                                                                                                               

    fd = open (argv[1], O_RDONLY);                                                                                                                                                                                                                                                      
    if(fd < 0) {                                                                                                                                                                                                                                                                        
        fprintf(stderr,"Cannot open: %s %s.\n",argv[1],strerror(errno));                                                                                                                                                                                                                
        return EXIT_FAILURE;                                                                                                                                                                                                                                                            
        }                                                                                                                                                                                                                                                                               
    len = atoi(argv[2]);                                                                                                                                                                                                                                                                

    buff=(char*)malloc(BUF_STDOUT);                                                                                                                                                                                                                                                     
    if(buff==NULL) {                                                                                                                                                                                                                                                                    
        fprintf(stderr,"Out of memory\n");                                                                                                                                                                                                                                              
        return EXIT_FAILURE;                                                                                                                                                                                                                                                            
        }                                                                                                                                                                                                                                                                               
    setvbuf(stdout, buff, _IOFBF, BUF_STDOUT);                                                                                                                                                                                                                                          

    if(len < 2) {                                                                                                                                                                                                                                                                       
        fprintf(stderr,"bad length %s.\n",argv[2]);                                                                                                                                                                                                                                     
        return EXIT_FAILURE;                                                                                                                                                                                                                                                            
        }                                                                                                                                                                                                                                                                               
    /* Get the size of the file. */                                                                                                                                                                                                                                                     
    if(fstat (fd, & s)!=0) {                                                                                                                                                                                                                                                            
        fprintf(stderr,"Cannot stat: %s %s.\n",argv[1],strerror(errno));                                                                                                                                                                                                                
        return EXIT_FAILURE;                                                                                                                                                                                                                                                            
        }                                                                                                                                                                                                                                                                               
    size = s.st_size;                                                                                                                                                                                                                                                                   

    seq = (char *) mmap (0, size, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, 0);                                                                                                                                                                                                        
    x=0;                                                                                                                                                                                                                                                                                
    while(x<size) {                                                                                                                                                                                                                                                                     
        if(seq[x]=='>') {                                                                                                                                                                                                                                                               
                size_t x0=x;                                                                                                                                                                                                                                                            
                DUMP;                                                                                                                                                                                                                                                                   
                free(seq_name);                                                                                                                                                                                                                                                         
                while(seq[x]!='\n' && x < size) x++;    
```

### Nucleotide composition

To determine the number of each nucleotide in each motif and control, we used the following script:

```bash
#!/bin/bash
#SBATCH --job-name=matthias_calculate_nuc_content
#SBATCH --output=matthias-%j.out 
#SBATCH --error=matthias-%j.err 
#SBATCH --mem-per-cpu=16G 
#SBATCH -n 2


bed_folder="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/scripts/error_detection_pipeline/input"
output_folder="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/scripts/nucleotide_content/output"
input_list=$1

while IFS= read -r file; do
ID=$(echo $file |  cut -f 1-3 -d '.')
echo $ID
/galaxy/home/mxw1081/software/bedtools2/bin/bedtools nuc -fi /nfs/brubeck.bx.psu.edu/scratch6/mxw1081/references/hg19/hg19_formated_by_wil.fa -bed  ${bed_folder}/$file > ${output_folder}/${ID}.nuc_content  
done<$input_list


while IFS= read -r file; do
ID=$(echo $file |  cut -f 1-3 -d '.')
echo $ID
  tail -n +2 ${output_folder}/${ID}.nuc_content | cut -f 1-9 > ${output_folder}/temp
 rm ${output_folder}/${ID}.nuc_content
 mv ${output_folder}/temp  ${output_folder}/${ID}.nuc_content
done<$input_list
```

## STR annotation

To annotate STRs (simple tandem repeats or microsatellites) across the genome, an approach developed by Fungtamassan et al. 2015 was taken. At its core is the `microsatellite.py` script, which was used in the following wrapper (`reference.profiling.sh`)

```bash
#!/bin/bash
#SBATCH --job-name=str_annotation
#SBATCH --output=matthias-%j.out
#SBATCH --error=matthias-%j.err
#SBATCH --mem-per-cpu=16G
#SBATCH -n 1

chr_folder="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/references/hg19"
array=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)

for chr in "${array[@]}"; do 
echo $chr

INPUT=$chr
echo " "
echo " detect STR in reference genome" ## See detail in microsatellite.xml on https://github.com/Arkarachai/STR-FM
python2 microsatellite.py ${chr_folder}/${INPUT}.fa --fasta --period=1 --partialmotifs --minlength=4 --prefix=0 --suffix=0 --hamming=0 --multipleruns --flankdisplay=0  --splitbyvalidity  >${INPUT}.mono.out
python2 microsatellite.py ${chr_folder}/${INPUT}.fa --fasta --period=2 --partialmotifs --minlength=6 --prefix=0 --suffix=0 --hamming=0 --multipleruns --flankdisplay=0  --splitbyvalidity  >${INPUT}.di.out
python2 microsatellite.py ${chr_folder}/${INPUT}.fa --fasta --period=3 --partialmotifs --minlength=6 --prefix=0 --suffix=0 --hamming=0 --multipleruns --flankdisplay=0  --splitbyvalidity  >${INPUT}.tri.out
python2 microsatellite.py ${chr_folder}/${INPUT}.fa --fasta --period=4 --partialmotifs --minlength=8 --prefix=0 --suffix=0 --hamming=0 --multipleruns --flankdisplay=0  --splitbyvalidity  >${INPUT}.tetra.out

echo "formatting"
cat ${INPUT}.mono.out | awk 'BEGIN{FS="\t";OFS="\t"};{print $6,$2,$2+$1,$4,$1,length($4) }' > ${INPUT}.mono.TR
cat ${INPUT}.di.out | awk 'BEGIN{FS="\t";OFS="\t"};{print $6,$2,$2+$1,$4,$1,length($4) }' > ${INPUT}.di.TR
cat ${INPUT}.tri.out | awk 'BEGIN{FS="\t";OFS="\t"};{print $6,$2,$2+$1,$4,$1,length($4) }' > ${INPUT}.tri.TR
cat ${INPUT}.tetra.out | awk 'BEGIN{FS="\t";OFS="\t"};{print $6,$2,$2+$1,$4,$1,length($4) }' > ${INPUT}.tetra.TR
done
```

The output was then filtered by chromosome:

```bash
#!/bin/bash

# Take the output from the Fungtammasan et al 2015 script, and filter it as follows:
# mono STRs >= 8 units, di: >= 4, tri >= 3, tetra >= 3
# then also add 5 bp up and downstream to increase overlap when interesecting with motis

array=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)

for chromosome in "${array[@]}"; do 
echo $chromosome
rm ${chromosome}.STRs_filtered.bed
awk '$5 >=8' ${chromosome}.mono.TR | awk '{print $1 "\t" $2-5 "\t" $3+5}' >> ${chromosome}.STRs_filtered_5bp_flanks.bed
awk '$5 >=4' ${chromosome}.di.TR | awk '{print $1 "\t" $2-5 "\t" $3+5}' >> ${chromosome}.STRs_filtered_5bp_flanks.bed
awk '$5 >=3' ${chromosome}.tri.TR | awk '{print $1 "\t" $2-5 "\t" $3+5}' >> ${chromosome}.STRs_filtered_5bp_flanks.bed
awk '$5 >=3' ${chromosome}.tetra.TR | awk '{print $1 "\t" $2-5 "\t" $3+5}' >> ${chromosome}.STRs_filtered_5bp_flanks.bed
sort -k 2 -n ${chromosome}.STRs_filtered_5bp_flanks.bed > temp
rm ${chromosome}.STRs_filtered_5bp_flanks.bed
mv temp ${chromosome}.STRs_filtered_5bp_flanks.bed
done
```

## False-positive SNPs due to errors

At the core of this analysis is a function developed in R that takes the following parameters: an error rate, number of haploid genomes, average read depth, minimum number of reads to detect a variant, and the length of a given sequence. Based on these parameters, the number of false-positive SNPs due to sequencing errors can then be estimated. By performing these calculations per site (as opposed to averaged over the entire length of the sequence) and for a number of iterations, it is possible to capture the variation in the number of false positives instead of yielding just a point estimate.

Since doing this for a large number of sites (Mb), genomes, and iterations is computationally intensive, we've developed an approach to parallelize the process, by splitting each iteration of the simulation up into a separate job.

This is the R script that is used to perform the calculation per job:

```r
args = commandArgs(trailingOnly=TRUE)                                                                                                                                                                                                                                                                    


variantsDueToErrors_simulation <- function(rate, genomes, depth, min_reads, n_sites){                                                                                                                                                                                                                    
  # Compute the probability of observing a variant in a site due to sequencing error                                                                                                                                                                                                                     
  # rate: sequencing error rate per read base                                                                                                                                                                                                                                                            
  # genomes: number of individual haploid genomes                                                                                                                                                                                                                                                        
  # depth: average sequencing depth                                                                                                                                                                                                                                                                      
  # min_reads: minimum number of reads supporting the variant                                                                                                                                                                                                                                            
  # n_sites: number of sites corresponding to that error rate                                                                                                                                                                                                                                            

  # Sequencing error in a nucleotide: X ~ B(1,rate)                                                                                                                                                                                                                                                      
  # Assume independence among errors of different sequenced nucleotides corresponding to the same site                                                                                                                                                                                                   
  # Assume independence among sequencing errors of different haploid genomes                                                                                                                                                                                                                             
  # Number of errors on a site: Y ~ B(depth,rate)                                                                                                                                                                                                                                                        

  Y = matrix(rbinom(n = genomes*n_sites, size = depth, prob = rate), nrow = genomes, ncol = n_sites)                                                                                                                                                                                                     

  # We call a variant in an haplononBtype because of sequencing errors if Y >= min_reads                                                                                                                                                                                                                 
  # Number of haplononBtypes with a variant in a site due to sequencing errors V                                                                                                                                                                                                                         

  V = colSums( Y >= min_reads )                                                                                                                                                                                                                                                                          

  return(V)                                                                                                                                                                                                                                                                                              
}                                                                                                                                                                                                                                                                                                        

error_rate = as.numeric(args[1])                                                                                                                                                                                                                                                                         
n_sites = as.numeric((args[2]))                                                                                                                                                                                                                                                                          
depth = as.integer(args[3])                                                                                                                                                                                                                                                                              
n_genomes = as.integer(args[4])                                                                                                                                                                                                                                                                          
job=as.character(args[5])                                                                                                                                                                                                                                                                                
nonBtype=as.character(args[6])                                                                                                                                                                                                                                                                           


V = matrix(NA, nrow = 1, ncol = n_sites )                                                                                                                                                                                                                                                                
  V[1,] = variantsDueToErrors_simulation(rate = error_rate, genomes = n_genomes, depth = depth, min_reads = depth/3, n_sites = n_sites)                                                                                                                                                                  

min_variants = 1                                                                                                                                                                                                                                                                                         
exp_n_variants_singletons = rowSums( V >= min_variants )                                                                                                                                                                                                                                                 
min_variants = 0.01*n_genomes                                                                                                                                                                                                                                                                            
exp_n_variants_maf_0.01 = rowSums( V >= min_variants )                                                                                                                                                                                                                                                   
min_variants = 0.05*n_genomes                                                                                                                                                                                                                                                                            
exp_n_variants_maf_0.05 = rowSums( V >= min_variants )                                                                                                                                                                                                                                                   


cat(c(nonBtype, error_rate, n_sites, n_genomes, depth, exp_n_variants_singletons, exp_n_variants_maf_0.01, exp_n_variants_maf_0.05))                                                                                                                                                                     
```

The output of this file is simply one line that lists the nonB type, the error rate used in the calculation, the length of the sequence, the number of genomes, and the expected FP variants for singletons, MAF = 0.01, and MAF = 0.05

This is the script to submit one job:

```bash
#!/bin/bash                                                                                                                                                                                                                                                                                              
#SBATCH --job-name=matthias_map_ont_reads                                                                                                                                                                                                                                                                
#SBATCH --output=matthias-%j.out                                                                                                                                                                                                                                                                         
#SBATCH --error=matthias-%j.err                                                                                                                                                                                                                                                                          
#SBATCH -C new #SBATCH --nodes=1                                                                                                                                                                                                                                                                         
#SBATCH --ntasks=32                                                                                                                                                                                                                                                                                      
#SBATCH --mem=64G                                                                                                                                                                                                                                                                                        

#time Rscript /nfs/brubeck.bx.psu.edu/scratch6/mxw1081/scripts/false_positive_snps/probabilistic_simulation.r                                                                                                                                                                                            
rate=$1                                                                                                                                                                                                                                                                                                  
length=$2                                                                                                                                                                                                                                                                                                
depth=$3                                                                                                                                                                                                                                                                                                 
genomes=$4                                                                                                                                                                                                                                                                                               
job=$5                                                                                                                                                                                                                                                                                                   
nonBtype=$6                                                                                                                                                                                                                                                                                              

Rscript --vanilla  probabilistic_simulation_parallel.r $rate $length $depth $genomes $job $nonBtype    
```

And this is the master batch script that splits the whole task up in many single jobs:

```bash
#!/bin/bash                                                                                                                                                                                                                                                                                              
#SBATCH --job-name=matthias_map_ont_reads                                                                                                                                                                                                                                                                
#SBATCH --output=matthias-%j.out                                                                                                                                                                                                                                                                         
#SBATCH --error=matthias-%j.err                                                                                                                                                                                                                                                                          
#SBATCH -C new #SBATCH --nodes=1                                                                                                                                                                                                                                                                         

number_of_genomes=100                                                                                                                                                                                                                                                                                    
while IFS= read -r line; do                                                                                                                                                                                                                                                                              
        #echo $line                                                                                                                                                                                                                                                                                      
        nonBtype=$(echo $line | awk '{print $1}')                                                                                                                                                                                                                                                        
        error_rate=$(echo $line | awk '{print $2}')                                                                                                                                                                                                                                                      
        length=$(echo $line | awk '{print $3}')                                                                                                                                                                                                                                                          
        for depth in {3..30..3}                                                                                                                                                                                                                                                                          
           do                                                                                                                                                                                                                                                                                            
           #echo $depth                                                                                                                                                                                                                                                                                  
                for job in {1..1000}                                                                                                                                                                                                                                                                     
                do                                                                                                                                                                                                                                                                                       
                     echo  $job                                                                                                                                                                                                                                                                          
                     sbatch --exclude nn[4] --nice run_rscript.sh $error_rate $length $depth $number_of_genomes $job $nonBtype                                                                                                                                                                           
                     sleep 0.05                                                                                                                                                                                                                                                                          
                done                                                                                                                                                                                                                                                                                     
        done                                                                                                                                                                                                                                                                                             

done<input_simulations   
```

The `input_simulations`file looks like this, where the first column denotes the non-B type, the second is the estimated error rate, and the third is the total length covered by the respective type.

```
APhasedRepeats  0.002001064     5833024                                                                                                                                                                                                                                                                  
DirectRepeats   0.002299526     5688724                                                                                                                                                                                                                                                                  
G4Motifs        0.004963612     6920031                                                                                                                                                                                                                                                                  
InvertedRepeats 0.002008671     71619329                                                                                                                                                                                                                                                                 
MirrorRepeats   0.002048599     17869244                                                                                                                                                                                                                                                                 
ZDNAMotifs      0.003333478     1689492                                                                                                                                                                                                                                                                  
BDNA    0.002190000     109619844     
```

### Figures false-positive SNPs

After concatenating the output of the individual jobs, we get a data frame, which then serves as the input for the figures.

### G-runs in G4 motifs

Since middle Gs in GGGs seem to have a substantially elevated error rate, we wanted to investigate these separately.

First we used the same wrapper as in the STR identification described above to identify G-runs of length >= 3.

```bash
#!/bin/bash                                                                                                                                                                                                                                                                     
#SBATCH --job-name=str_annotation                                                                                                                                                                                                                                               
#SBATCH --output=matthias-%j.out                                                                                                                                                                                                                                                
#SBATCH --error=matthias-%j.err                                                                                                                                                                                                                                                 
#SBATCH --mem-per-cpu=16G                                                                                                                                                                                                                                                       
#SBATCH -n 1                                                                                                                                                                                                                                                                    

chr_folder="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/references/hg19"                                                                                                                                                                                                           
array=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)                                                                                                                                                                                                                

for chr in "${array[@]}"; do                                                                                                                                                                                                                                                    
echo $chr                                                                                                                                                                                                                                                                       

INPUT=$chr                                                                                                                                                                                                                                                                      
echo " "                                                                                                                                                                                                                                                                        
echo " detect STR in reference genome" ## See detail in microsatellite.xml on https://github.com/Arkarachai/STR-FM                                                                                                                                                              
python2 microsatellite.py ${chr_folder}/${INPUT}.fa --fasta --period=1 --partialmotifs --minlength=3 --prefix=0 --suffix=0 --hamming=0 --multipleruns --flankdisplay=0  --splitbyvalidity  >${INPUT}.GGGs.out                                                                   

echo "formatting"                                                                                                                                                                                                                                                               
cat ${INPUT}.GGGs.out | awk 'BEGIN{FS="\t";OFS="\t"};{print $6,$2,$2+$1,$4,$1,length($4) }' > ${INPUT}.GGG.TR                                                                                                                                                                   
done 
```

In this script, we curate the output so we only keep G-runs:

```bash
#!/bin/bash                                                                                                                                                                                                                                                                     

for file in $(ls *GGGs.out); do                                                                                                                                                                                                                                                 
        chr_id=$(echo $file | cut -f 1 -d '.')                                                                                                                                                                                                                                  
        awk '$4=="G"' $file | awk -v c=$chr_id '{print c "\t" $2 "\t" $2 + $1}'                                                                                                                                                                                                 
done  
```

Finally, we use `bedtools coverage` to find overlaps between the G4 motif annotation and the G-tracts: the first awk command filters out any G-runs that don't have any overlap with a G4 motif, in the second identifies the number of 'middle Gs' (length of G-run -2) , and the last one sums up the number of middle Gs.

```bash
/galaxy/home/mxw1081/software/bedtools2/bin/bedtools coverage -wa -b G4s_hg19.sorted.bed -a GGGs_hg19.sorted.bed  | awk '$7!=0' | awk '{print $5-2}' | awk '{sum+=$1}; END {print sum}'
```

Intersected with the G4 motifs from Illumina_curated_fully_filterd (which already has the mappability = 1 filter) this yields 1,715,082 potential middle Gs in G4 motifs.

### Empirical false-positive SNV estimation

In addition, we estimated the proportion of false-positive SNVs using the same Illumina data set that was also employed for the error detection analyses. To illustrate the effect of varying read depths on false-positives, we sub-sampled the HG002 100X data set to approximately 3, 9, 15, and 30¢·, and performed variant calling using the tool FreeBayes with default parameters and a minimum mapping quality filter set to 30.

Here is an exemplary script of how to subset and identify variants:

```bash
#!/bin/bash                                                                                                                                                                       
#SBATCH --job-name=matthias_map_ont_reads                                                                                                                                         
#SBATCH --output=matthias-%j.out                                                                                                                                                  
#SBATCH --error=matthias-%j.err                                                                                                                                                   
#SBATCH --nodes=1                                                                                                                                                                 

export CONDA_ENVS_PATH="/galaxy/home/mxw1081/conda_envs"                                                                                                                          
source /opt/anaconda/etc/profile.d/conda.sh                                                                                                                                       
conda activate false_positive_snp_calling_nonB                                                                                                                                    


bam_folder="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/bam_files/Illumina"                                                                                                          
reference_folder="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/references/hg19"                                                                                                       

output_folder="/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/scripts/false_positive_snps_nonB/vcf_files"                                                                               

#freebayes -f ${input_folder}/1.fa ${input_folder}/HG002.hs37d5.100X.sorted.rmdup.chr1.5X.bam > HG002.hs37d5.100X.sorted.rmdup.chr1.5X.bam.var.vcf                                

#freebayes --min-base-quality 30 -f ${input_folder}/1.fa ${input_folder}/HG002.hs37d5.100X.sorted.rmdup.chr1.5X.bam > HG002.hs37d5.100X.sorted.rmdup.chr1.5X.minBQ30.bam.var.vcf  

samtools view -s 0.01 -b ${bam_folder}/HG002.hs37d5.100x.sorted.rmdup.bam | \                                                                                                     
        freebayes  \                                                                                                                                                              
                  --min-mapping-quality 30 \                                                                                                                                      
                  -f ${reference_folder}/hg19_formated_by_wil.fa -c > ${output_folder}/HG002.hs37d5.100x.sorted.rmdup.bam.3X.MQ30_only_filter.vcf   
```

The output of this script (a vcf file that contains all the variants) is then taken to identify the false-positive variants:

```bash
#!/bin/bash                                                                                                                                                                                       

export CONDA_ENVS_PATH="/galaxy/home/mxw1081/conda_envs"                                                                                                                                          
source /opt/anaconda/etc/profile.d/conda.sh                                                                                                                                                       
conda activate false_positive_snp_calling_nonB                                                                                                                                                    

#TP_snps_motif="HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.breakmulti.inbed.SNPS_ONLY.OVERLAP_nonB.vcf"                                    
#TP_snps_control="HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.breakmulti.inbed.SNPS_ONLY.OVERLAP_control.vcf"                               

TP_snps_motif="HG002_GIAB_benchmark_v4.2.1_SNPs_only.OVERLAP_nonB.vcf"                                                                                                                            
TP_snps_control="HG002_GIAB_benchmark_v4.2.1_SNPs_only.OVERLAP_control.vcf"                                                                                                                       

#TP_snps="HG002.hs37d5.100x.sorted.rmdup.bam.100X.MQ30_only_filter.OVERLAP_nonB.vcf"                                                                                                              
#echo 'total number of true SNPS in non-B motifs'                                                                                                                                                 
#                                                                                                                                                                                                 
#bedtools intersect -a $TP_snps -b nonB_annotation_all_types_MERGED.bed | wc -l                                                                                                                   

depth=$1                                                                                                                                                                                          

total_n_motif=$(awk 'length($5)==1 && length($4)==1' vcf_files/HG002.hs37d5.100x.sorted.rmdup.bam.${depth}X.MQ30_only_filter.vcf | \                                                              
        cat <(grep "#" vcf_files/HG002.hs37d5.100x.sorted.rmdup.bam.${depth}X.MQ30_only_filter.vcf ) - | \                                                                                        
        bedtools intersect -a - -b nonB_annotation_all_types_MERGED.bed | wc -l)                                                                                                                  
echo -e 'total_n' '\t' $total_n_motif '\t' 'motif' '\t' ${depth} '\t' '129811185'                                                                                                                 

true_positives_motif=$(awk 'length($5)==1 && length($4)==1' vcf_files/HG002.hs37d5.100x.sorted.rmdup.bam.${depth}X.MQ30_only_filter.vcf | \                                                       
        cat <(grep "#" vcf_files/HG002.hs37d5.100x.sorted.rmdup.bam.${depth}X.MQ30_only_filter.vcf ) - | \                                                                                        
        bedtools intersect -a - -b nonB_annotation_all_types_MERGED.bed | \                                                                                                                       
        cat <(grep "#" vcf_files/HG002.hs37d5.100x.sorted.rmdup.bam.${depth}X.MQ30_only_filter.vcf ) - | \                                                                                        
        bedtools intersect -a $TP_snps_motif -b - | wc -l)                                                                                                                                        

echo -e 'tp' '\t' $true_positives_motif '\t' 'motif' '\t' ${depth} '\t' '129811185'                                                                                                               

false_positives_motif=$(awk 'length($5)==1 && length($4)==1' vcf_files/HG002.hs37d5.100x.sorted.rmdup.bam.${depth}X.MQ30_only_filter.vcf | \                                                      
        cat <(grep "#" vcf_files/HG002.hs37d5.100x.sorted.rmdup.bam.${depth}X.MQ30_only_filter.vcf ) - | \                                                                                        
        bedtools intersect -a - -b nonB_annotation_all_types_MERGED.bed | \                                                                                                                       
        cat <(grep "#" vcf_files/HG002.hs37d5.100x.sorted.rmdup.bam.${depth}X.MQ30_only_filter.vcf ) - | \                                                                                        
        bedtools intersect -v  -b $TP_snps_motif -a - | wc -l)                                                                                                                                    

echo -e 'fp' '\t' $false_positives_motif '\t' 'motif' '\t' ${depth} '\t' '129811185' 
```

## Downstream analysis

### Initial filtering and parsing data frames

Before proceeding to the downstream analysis, data frames containing error rates are filtered and linked with other sequence featuers

```r
Sys.getenv("YOUR_VAR")
library(easypackages)
libraries("ggplot2","data.table", "patchwork", "car", "cowplot", "scales","patchwork", "plyr", "tidyr","msme", "dplyr","vroom", "glmmADMB", "purrr", "grid", "gridExtra","cowplot", "lme4", "MatchIt", "compositions")
`%not_in%` <- purrr::negate(`%in%`)
setwd("/Users/matthiasweissensteiner/Dropbox/PostDoc/proj/nonB/prelim_results/new_results")

homopolymers <- vroom("homopolymers_annotated.bed", col_names = FALSE)
data.frame(homopolymers) ->  homopolymers
names(homopolymers) <- c("chr", "start", "end", "number_homopolymers", "proportion_homopolymers")
motif_list<-c("APhasedRepeats", "DirectRepeats", "G4Motifs", "InvertedRepeats", "MirrorRepeats", "ZDNAMotifs")

mappability_36mer <- vroom("mappability_36mer_output_parsed", col_names = FALSE)
data.frame(mappability_36mer) -> mappability_36mer
names(mappability_36mer) <- c("id", "mappability_36mer")

overlap_ids <- read.table("overlap_nonB_motifs_ids", header = F)
overlap_ids %>%
  select(V1) %>%
  unlist() %>% as.vector() -> overlap_ids

rm_overlap_ids <- read.table("RM_hg19_overlap_ids", header = F)
rm_overlap_ids %>%
  select(V1) %>%
  unlist() %>% as.vector() -> rm_overlap_ids

str_ids <- read.table("STR_overlap_ids", header = FALSE)
str_ids %>%
  select(V1) %>%
  unlist() %>% as.vector() -> str_overlap_ids

G4plus <- read.table("G4_plus_IDs", header = FALSE)
G4plus %>%
  select(V1) %>%
  unlist() %>% as.vector() -> G4_plus_ids

overlap_clusters <- read.table("no_clustering_IDs", header = FALSE)
overlap_clusters %>%
  select(V1) %>%
  unlist() %>% as.vector() -> no_clustering_IDs


########################################################################################################################
# Load data frames
########################################################################################################################

 Illumina
####
Illumina<- vroom("Illumina_fwd_final_output_nuc_content_transformed.table", col_names = T)
data.frame(Illumina) -> Illumina_curated
Illumina_curated$type <- as.factor(Illumina_curated$type)
Illumina_curated$source <- as.factor(Illumina_curated$source)
Illumina_curated$source <- with(Illumina_curated, relevel(source, "motif"))
Illumina_curated %>% right_join(homopolymers, by = c("chr", "start", "end")) -> Illumina_curated
Illumina_curated %>% right_join(mappability_36mer, by = c("id")) -> Illumina_curated

# Remove sites overlapping with RM hits, with mappability < 1, and length larger than 1kb
Illumina_curated %>%
  mutate(overlap_with_RM = ifelse(id %in% rm_overlap_ids, 1,0),
         overlap_with_STR = ifelse(id %in% str_overlap_ids, 1, 0)) %>%
  filter(
           mappability_36mer == 1 &
           length < 1000) -> Illumina_curated

# Subsample data frame to match numbers of motifs and controls

# First check which (motif or control) sample size is higher: 
Illumina_curated %>% 
  group_by(type, source) %>%
  tally() %>% 
  dcast(., type~source) %>% 
  mutate(control_higher = ifelse(control > motif, 1,0), 
         subsample_n = ifelse(control_higher == 1, motif, control)) -> sample_size_table

# Aphased Repeats:
Illumina_curated %>%
  filter(source == "motif" & type == "APhasedRepeats") %>% 
  sample_n(., 261798) -> temp
Illumina_curated %>%
  filter(source == "control" & type == "APhasedRepeats") %>%
  rbind.data.frame(., temp) -> AphasedRepeats

# Direct Repeats:
Illumina_curated %>%
  filter(source == "control" & type == "DirectRepeats") %>% 
  sample_n(., 360335) -> temp
Illumina_curated %>%
  filter(source == "motif" & type == "DirectRepeats") %>%
  rbind.data.frame(., temp) -> DirectRepeats

# G4 Motifs:

Illumina_curated %>%
  filter(type == "G4Motifs" & source == "motif" ) %>% 
  mutate(G4_plus = ifelse(id %in% G4_plus_ids, 1, 0)) %>% 
  filter(G4_plus == 1) %>% dim()
Illumina_curated %>%
  filter(source == "control" & type == "G4Motifs") %>% 
  sample_n(., 243017) -> temp
Illumina_curated %>%
  filter(source == "motif" & type == "G4Motifs") %>%
  mutate(G4_plus = ifelse(id %in% G4_plus_ids, 1, 0)) %>% 
  filter(G4_plus == 1) %>% 
  select(-G4_plus) %>%
  rbind.data.frame(., temp) -> G4Motifs

# Inverted Repeats:
Illumina_curated %>%
  filter(source == "motif" & type == "InvertedRepeats") %>% 
  sample_n(., 3774639) -> temp
Illumina_curated %>%
  filter(source == "control" & type == "InvertedRepeats") %>%
  rbind.data.frame(., temp) -> InvertedRepeats

# Mirror Repeats:
Illumina_curated %>%
  filter(source == "control" & type == "MirrorRepeats") %>% 
  sample_n(., 498033) -> temp
Illumina_curated %>%
  filter(source == "motif" & type == "MirrorRepeats") %>%
  rbind.data.frame(., temp) -> MirrorRepeats

# ZDNA Motifs:
Illumina_curated %>%
  filter(source == "control" & type == "ZDNAMotifs") %>% 
  sample_n(., 222534) -> temp
Illumina_curated %>%
  filter(source == "motif" & type == "ZDNAMotifs") %>%
  rbind.data.frame(., temp) -> ZDNAMotifs

rbind.data.frame(AphasedRepeats,
                 DirectRepeats,
                 G4Motifs,
                 InvertedRepeats,
                 MirrorRepeats,
                 ZDNAMotifs) -> Illumina_curated

# Add homopolymer,
# rates, and overlap with other motifs variables

Illumina_curated %>% 
  select(-smm_rate, -ins_rate, -del_rate) %>% 
  mutate(proportion_homopolymers_binary = factor(ifelse(proportion_homopolymers == 0, 0, 1)),
         smm_rate = smm / depth,
         ins_rate = ins / depth,
         del_rate = del / depth,
         overlap_with_other_nonB = ifelse(id %in% overlap_ids, 1,0),
         overlap_groups = as.factor(ifelse(source == "control", 3, ifelse(source == "motif" & overlap_with_other_nonB == 1, 2, 1)))
         ) -> Illumina_curated

Illumina_curated %>%
  write.table(., file = "Illumina_curated_fully_filtered", sep = "\t", quote = FALSE, row.names = FALSE)

Illumina_curated <- vroom("Illumina_curated_fully_filtered", col_names  = T)
data.frame(Illumina_curated) -> Illumina_curated
Illumina_curated$type <- as.factor(Illumina_curated$type)
Illumina_curated$source <- as.factor(Illumina_curated$source)
Illumina_curated$source <- with(Illumina_curated, relevel(source, "motif"))

#####

# HiFi
#####
HiFi<- vroom("HiFi_fwd_final_output_nuc_content_transformed.table", col_names = T)
data.frame(HiFi) -> HiFi_curated
HiFi_curated$type <- as.factor(HiFi_curated$type)
HiFi_curated$source <- as.factor(HiFi_curated$source)
HiFi_curated$source <- with(HiFi_curated, relevel(source, "motif"))
HiFi_curated %>% right_join(homopolymers, by = c("chr", "start", "end")) -> HiFi_curated
HiFi_curated %>% right_join(mappability_36mer, by = c("id")) -> HiFi_curated

# Remove sites overlapping with RM hits, with mappability < 1, and length larger than 1kb
HiFi_curated %>%
  mutate(overlap_with_RM = ifelse(id %in% rm_overlap_ids, 1,0),
         overlap_with_STR = ifelse(id %in% str_overlap_ids, 1, 0)) %>%
  filter(
    mappability_36mer == 1 &
      length < 1000) -> HiFi_curated

# Subsample data frame to match numbers of motifs and controls

# First check which (motif or control) sample size is higher: 
HiFi_curated %>% 
  group_by(type, source) %>%
  tally() %>% 
  dcast(., type~source) %>% 
  mutate(control_higher = ifelse(control > motif, 1,0), 
         subsample_n = ifelse(control_higher == 1, motif, control)) -> sample_size_table

# Aphased Repeats:
HiFi_curated %>%
  filter(source == "motif" & type == "APhasedRepeats") %>% 
  sample_n(., 261798) -> temp
HiFi_curated %>%
  filter(source == "control" & type == "APhasedRepeats") %>%
  rbind.data.frame(., temp) -> AphasedRepeats

# Direct Repeats:
HiFi_curated %>%
  filter(source == "control" & type == "DirectRepeats") %>% 
  sample_n(., 360335) -> temp
HiFi_curated %>%
  filter(source == "motif" & type == "DirectRepeats") %>%
  rbind.data.frame(., temp) -> DirectRepeats

# G4 Motifs:
HiFi_curated %>%
  filter(type == "G4Motifs" & source == "motif" ) %>% 
  mutate(G4_plus = ifelse(id %in% G4_plus_ids, 1, 0)) %>% 
  filter(G4_plus == 1) %>% dim()
HiFi_curated %>%
  filter(source == "control" & type == "G4Motifs") %>% 
  sample_n(., 243017) -> temp
HiFi_curated %>%
  filter(source == "motif" & type == "G4Motifs") %>%
  mutate(G4_plus = ifelse(id %in% G4_plus_ids, 1, 0)) %>% 
  filter(G4_plus == 1) %>% 
  select(-G4_plus) %>%
  rbind.data.frame(., temp) -> G4Motifs

# Inverted Repeats:
HiFi_curated %>%
  filter(source == "motif" & type == "InvertedRepeats") %>% 
  sample_n(., 3774639) -> temp
HiFi_curated %>%
  filter(source == "control" & type == "InvertedRepeats") %>%
  rbind.data.frame(., temp) -> InvertedRepeats

# Mirror Repeats:
HiFi_curated %>%
  filter(source == "control" & type == "MirrorRepeats") %>% 
  sample_n(., 498033) -> temp
HiFi_curated %>%
  filter(source == "motif" & type == "MirrorRepeats") %>%
  rbind.data.frame(., temp) -> MirrorRepeats

# ZDNA Motifs:
HiFi_curated %>%
  filter(source == "control" & type == "ZDNAMotifs") %>% 
  sample_n(., 222534) -> temp
HiFi_curated %>%
  filter(source == "motif" & type == "ZDNAMotifs") %>%
  rbind.data.frame(., temp) -> ZDNAMotifs

rbind.data.frame(AphasedRepeats,
                 DirectRepeats,
                 G4Motifs,
                 InvertedRepeats,
                 MirrorRepeats,
                 ZDNAMotifs) -> HiFi_curated

# Add homopolymer,
# rates, and overlap with other motifs variables

HiFi_curated %>% 
  mutate(proportion_homopolymers_binary = factor(ifelse(proportion_homopolymers == 0, 0, 1)),
         smm_rate = smm / depth,
         ins_rate = ins / depth,
         del_rate = del / depth,
         overlap_with_other_nonB = ifelse(id %in% overlap_ids, 1,0),
         overlap_groups = as.factor(ifelse(source == "control", 3, ifelse(source == "motif" & overlap_with_other_nonB == 1, 2, 1)))
         ) -> HiFi_curated

HiFi_curated %>%
  select(-smm_rate, -ins_rate, -del_rate) %>% 
  mutate(smm_rate = smm/depth,
         ins_rate = ins/depth,
         del_rate = del/depth) -> HiFi_curated

HiFi_curated %>%
  write.table(., file = "HiFi_curated_fully_filtered", sep = "\t", quote = FALSE, row.names = FALSE)

HiFi_curated <- vroom("HiFi_curated_fully_filtered", col_names = T)
data.frame(HiFi_curated) -> HiFi_curated
HiFi_curated$type <- as.factor(HiFi_curated$type)
HiFi_curated$source <- as.factor(HiFi_curated$source)
HiFi_curated$source <- with(HiFi_curated, relevel(source, "motif"))

#####

# ONT
#####
ONT<- vroom("ONT_fwd_final_output_nuc_content_transformed.table", col_names = T)
data.frame(ONT) -> ONT_curated
ONT_curated$type <- as.factor(ONT_curated$type)
ONT_curated$source <- as.factor(ONT_curated$source)
ONT_curated$source <- with(ONT_curated, relevel(source, "motif"))
ONT_curated %>% right_join(homopolymers, by = c("chr", "start", "end")) -> ONT_curated
ONT_curated %>% right_join(mappability_36mer, by = c("id")) -> ONT_curated

# Remove sites overlapping with RM hits, with mappability < 1, and length larger than 1kb
ONT_curated %>%
  mutate(overlap_with_RM = ifelse(id %in% rm_overlap_ids, 1,0),
         overlap_with_STR = ifelse(id %in% str_overlap_ids, 1, 0)) %>%
  filter(
    mappability_36mer == 1 &
      length < 1000) -> ONT_curated

# Subsample data frame to match numbers of motifs and controls

# First check which (motif or control) sample size is higher: 
ONT_curated %>% 
  group_by(type, source) %>%
  tally() %>% 
  dcast(., type~source) %>% 
  mutate(control_higher = ifelse(control > motif, 1,0), 
         subsample_n = ifelse(control_higher == 1, motif, control)) -> sample_size_table

# Aphased Repeats:
ONT_curated %>%
  filter(source == "motif" & type == "APhasedRepeats") %>% 
  sample_n(., 261798) -> temp
ONT_curated %>%
  filter(source == "control" & type == "APhasedRepeats") %>%
  rbind.data.frame(., temp) -> AphasedRepeats

# Direct Repeats:
ONT_curated %>%
  filter(source == "control" & type == "DirectRepeats") %>% 
  sample_n(., 360335) -> temp
ONT_curated %>%
  filter(source == "motif" & type == "DirectRepeats") %>%
  rbind.data.frame(., temp) -> DirectRepeats

# G4 Motifs:
ONT_curated %>%
  filter(type == "G4Motifs" & source == "motif" ) %>% 
  mutate(G4_plus = ifelse(id %in% G4_plus_ids, 1, 0)) %>% 
  filter(G4_plus == 1) %>% dim()
ONT_curated %>%
  filter(source == "control" & type == "G4Motifs") %>% 
  sample_n(., 243017) -> temp
ONT_curated %>%
  filter(source == "motif" & type == "G4Motifs") %>%
  mutate(G4_plus = ifelse(id %in% G4_plus_ids, 1, 0)) %>% 
  filter(G4_plus == 1) %>% 
  select(-G4_plus) %>%
  rbind.data.frame(., temp) -> G4Motifs

# Inverted Repeats:
ONT_curated %>%
  filter(source == "motif" & type == "InvertedRepeats") %>% 
  sample_n(., 3774639) -> temp
ONT_curated %>%
  filter(source == "control" & type == "InvertedRepeats") %>%
  rbind.data.frame(., temp) -> InvertedRepeats

# Mirror Repeats:
ONT_curated %>%
  filter(source == "control" & type == "MirrorRepeats") %>% 
  sample_n(., 498033) -> temp
ONT_curated %>%
  filter(source == "motif" & type == "MirrorRepeats") %>%
  rbind.data.frame(., temp) -> MirrorRepeats

# ZDNA Motifs:
ONT_curated %>%
  filter(source == "control" & type == "ZDNAMotifs") %>% 
  sample_n(., 222534) -> temp
ONT_curated %>%
  filter(source == "motif" & type == "ZDNAMotifs") %>%
  rbind.data.frame(., temp) -> ZDNAMotifs

rbind.data.frame(AphasedRepeats,
                 DirectRepeats,
                 G4Motifs,
                 InvertedRepeats,
                 MirrorRepeats,
                 ZDNAMotifs) -> ONT_curated

# Add homopolymer,
# rates, and overlap with other motifs variables

ONT_curated %>% 
  mutate(proportion_homopolymers_binary = factor(ifelse(proportion_homopolymers == 0, 0, 1)),
         smm_rate = smm / depth,
         ins_rate = ins / depth,
         del_rate = del / depth,
         overlap_with_other_nonB = ifelse(id %in% overlap_ids, 1,0),
         overlap_groups = as.factor(ifelse(source == "control", 3, ifelse(source == "motif" & overlap_with_other_nonB == 1, 2, 1)))
         ) -> ONT_curated

ONT_curated %>%
  select(-smm_rate, -ins_rate, -del_rate) %>% 
  mutate(smm_rate = smm/depth,
         ins_rate = ins/depth,
         del_rate = del/depth) -> ONT_curated

ONT_curated %>%
  write.table(., file = "ONT_curated_fully_filtered", sep = "\t", quote = FALSE, row.names = FALSE)

ONT<- vroom("ONT_curated_fully_filtered", col_names = T)
data.frame(ONT) -> ONT_curated
ONT_curated$type <- as.factor(ONT_curated$type)
ONT_curated$source <- as.factor(ONT_curated$source)
ONT_curated$source <- with(ONT_curated, relevel(source, "motif"))

#####
```

### Poisson regression modeling

The input for this script is the same filtered and curated table that is also used for producing figures. It is very similar for each technology and filtering levels, therefore only the one for Illumina and moderate filtering is shown here.

```r
library(plyr)
library(tidyr)
library(dplyr)
library(vroom)
library(reshape2)
library(MASS)

`%not_in%` <- purrr::negate(`%in%`)
setwd("/nfs/brubeck.bx.psu.edu/scratch6/mxw1081/scripts/negative_binomial_regression/2022_analysis/36mer_mappability_filter")

motif_list<-c("APhasedRepeats", "DirectRepeats", "G4Motifs", "InvertedRepeats", "MirrorRepeats", "ZDNAMotifs")

########################################################################################################################
# Technology-specific - Illumina
########################################################################################################################

Illumina<- vroom("Illumina_curated_fully_filtered", col_names = T)
data.frame(Illumina) -> Illumina_curated
Illumina_curated$type <- as.factor(Illumina_curated$type)
Illumina_curated$source <- as.factor(Illumina_curated$source)
Illumina_curated$source <- with(Illumina_curated, relevel(source, "motif"))


################################################### 
## Illumina Single-nucleotide mismatches - Generalized linear model with poisson model
#####


  deviance_expl_full_model <- c()
  deviance_expl_no_source_model <- c()
  deviance_expl_no_proportion_homopolymers_model <- c()
  deviance_expl_no_seq_comp_model<- c()
  deviance_expl_no_basequality_model<- c() 
  deviance_expl_no_overlap_model<- c() 
  deviance_expl_no_length_model<- c() 

  RCDE_no_source <- c()
  RCDE_no_proportion_homopolymers <- c() 
  RCDE_no_seq_comp <- c()
  RCDE_no_basequality <- c()
  RCDE_no_overlap <- c()
  RCDE_no_length <- c()

  BICs<-list()
  theta <- c()
  coefficients <- list()
  coefficients_pvalues <- list()
  all_poisson_models<-list()
  for (i in 1:length(motif_list)){
    print(motif_list[i])

    print("full")
    Illumina_curated %>% 
      filter(type == paste0(motif_list[i])) %>%
      na.omit() %>%
      glm(smm ~ factor(source) + 
            offset(log(depth)) +
            factor(overlap_with_other_nonB) + 
        factor(overlap_with_other_nonB)*log(length) + 
        log(length) + 
        factor(source)*log(length) + 
            V1 + V2 + V3 + basequality + 
            factor(proportion_homopolymers_binary),  
          data = .,
          family = poisson(link = "log")) -> full_poisson_model

      data.frame(cooks.distance(full_poisson_model)) %>%
      mutate(row_number = 1:dplyr::n()) %>%
      filter(cooks.distance.full_poisson_model. > 1) %>%
      dplyr::select(row_number) %>%
      unlist() %>%
      as.vector() -> influential_cases

      print("full")
    Illumina_curated %>% 
      filter(type == paste0(motif_list[i])) %>%
      na.omit() %>%
      mutate(row_number = 1:dplyr::n()) %>%
      filter(row_number %not_in% influential_cases) %>%
      glm(smm ~ factor(source) + 
            offset(log(depth)) +
            factor(overlap_with_other_nonB) + 
        factor(overlap_with_other_nonB)*log(length) + 
        log(length) + 
        factor(source)*log(length) + 
            V1 + V2 + V3 + basequality + 
            factor(proportion_homopolymers_binary),  
          data = .,
          family = poisson(link = "log")) -> full_poisson_model

    coefficients[[i]] <- data.frame(summary(full_poisson_model)$coefficients) %>% dplyr::select(Estimate) %>% unlist() %>% as.vector()
    coefficients_pvalues[[i]] <- data.frame(summary(full_poisson_model)$coefficients[,4])

    print("no_factor(source)")
    Illumina_curated %>% 
      filter(type == paste0(motif_list[i])) %>% 
      na.omit() %>%
      mutate(row_number = 1:dplyr::n()) %>%
      filter(row_number %not_in% influential_cases) %>%
      glm(smm ~ 
            offset(log(depth)) +
        log(length) + 
            V1 + V2 + V3 + basequality + 
            factor(proportion_homopolymers_binary),  
          data = .,
          family = poisson(link = "log")) -> no_source_poisson_model


    print("no_proportion_homopolymers") 
    Illumina_curated %>% 
      filter(type == paste0(motif_list[i])) %>% 
      na.omit() %>%
      mutate(row_number = 1:dplyr::n()) %>%
      filter(row_number %not_in% influential_cases) %>%
      glm(smm ~ factor(source) + 
            offset(log(depth)) +
            factor(overlap_with_other_nonB) + 
        factor(overlap_with_other_nonB)*log(length) + 
        log(length) + 
        factor(source)*log(length) + 
            V1 + V2 + V3 + basequality , 
          data = ., 
          family = poisson(link = "log"))-> no_proportion_homopolymers_poisson_model

    print("no_seq_comp")
    Illumina_curated %>% 
      filter(type == paste0(motif_list[i])) %>% 
      na.omit() %>%
      mutate(row_number = 1:dplyr::n()) %>%
      filter(row_number %not_in% influential_cases) %>%
      glm(smm ~ factor(source) + 
            offset(log(depth)) +
            factor(overlap_with_other_nonB) + 
        factor(overlap_with_other_nonB)*log(length) + 
        log(length) + 
        factor(source)*log(length) + 
        basequality + 
            factor(proportion_homopolymers_binary),  
          data = ., 
          family = poisson(link = "log")) -> no_seq_comp_poisson_model

    print("no_basequality")
    Illumina_curated %>% 
      filter(type == paste0(motif_list[i])) %>% 
      na.omit() %>%
      mutate(row_number = 1:dplyr::n()) %>%
      filter(row_number %not_in% influential_cases) %>%
      glm(smm ~ factor(source) + 
            offset(log(depth)) +
            factor(overlap_with_other_nonB) + 
        factor(overlap_with_other_nonB)*log(length) + 
        log(length) + 
        factor(source)*log(length) + 
            V1 + V2 + V3 + 
            factor(proportion_homopolymers_binary),  
          data = ., 
          family = poisson(link = "log"))-> no_basequality_poisson_model

    Illumina_curated %>% 
      filter(type == paste0(motif_list[i])) %>% 
      na.omit() %>%
      mutate(row_number = 1:dplyr::n()) %>%
      filter(row_number %not_in% influential_cases) %>%
      glm(smm ~ factor(source) + 
            offset(log(depth)) +
        log(length) + 
        factor(source)*log(length) + 
            V1 + V2 + V3 + basequality + 
            factor(proportion_homopolymers_binary),  
          data = .,
          family = poisson(link = "log")) -> no_overlap_poisson_model

    Illumina_curated %>% 
      filter(type == paste0(motif_list[i])) %>% 
      na.omit() %>%
      mutate(row_number = 1:dplyr::n()) %>%
      filter(row_number %not_in% influential_cases) %>%
      glm(smm ~ factor(source) + 
            offset(log(depth)) +
            factor(overlap_with_other_nonB) + 
            V1 + V2 + V3 + basequality + 
            factor(proportion_homopolymers_binary),  
          data = .,
          family = poisson(link = "log")) -> no_length_poisson_model


    poisson_models_per_motif<-list(full_poisson_model, 
                                   no_source_poisson_model,
                                   no_seq_comp_poisson_model,
                                   no_basequality_poisson_model,
                                   no_proportion_homopolymers_poisson_model,
                       no_overlap_poisson_model,
                                   no_length_poisson_model) 

    all_poisson_models[[i]]<-poisson_models_per_motif 

    null_deviance_full_model <- full_poisson_model$null.deviance
    ll.null.full_model <- full_poisson_model$null.deviance
    ll.residual.full_model <- full_poisson_model$deviance

    ll.residual.no_source_model <- no_source_poisson_model$deviance
    ll.null.no_source_model <- no_source_poisson_model$null.deviance

    ll.residual.no_proportion_homopolymers_model <- no_proportion_homopolymers_poisson_model$deviance
    ll.null.no_proportion_homopolymers_model <- no_proportion_homopolymers_poisson_model$null.deviance

    ll.residual.no_seq_comp_model <- no_seq_comp_poisson_model$deviance
    ll.null.no_seq_comp_model <- no_seq_comp_poisson_model$null.deviance

    ll.residual.no_basequality_model <- no_basequality_poisson_model$deviance
    ll.null.no_basequality_model <- no_basequality_poisson_model$null.deviance

    ll.residual.no_overlap_model <- no_overlap_poisson_model$deviance
    ll.null.no_overlap_model <- no_overlap_poisson_model$null.deviance

    ll.residual.no_length_model <- no_length_poisson_model$deviance
    ll.null.no_length_model <- no_length_poisson_model$null.deviance


    # Deviance explained by the full model:
    deviance_expl_full_model[i] <- 1-ll.residual.full_model/ll.null.full_model
    deviance_expl_no_source_model[i]<- 1-ll.residual.no_source_model/ll.null.no_source_model
    deviance_expl_no_proportion_homopolymers_model[i]<- 1-ll.residual.no_proportion_homopolymers_model/ll.null.no_proportion_homopolymers_model
    deviance_expl_no_seq_comp_model[i]<- 1-ll.residual.no_seq_comp_model/ll.null.no_seq_comp_model
    deviance_expl_no_basequality_model[i]<- 1-ll.residual.no_basequality_model/ll.null.no_basequality_model
    deviance_expl_no_overlap_model[i]<- 1-ll.residual.no_overlap_model/ll.null.no_overlap_model
    deviance_expl_no_length_model[i]<- 1-ll.residual.no_length_model/ll.null.no_length_model

    RCDE_no_source[i] <- ((null_deviance_full_model - ll.residual.full_model) - (null_deviance_full_model - ll.residual.no_source_model))/(null_deviance_full_model - ll.residual.full_model )
    RCDE_no_proportion_homopolymers[i] <- ((null_deviance_full_model - ll.residual.full_model) - (null_deviance_full_model - ll.residual.no_proportion_homopolymers_model))/(null_deviance_full_model - ll.residual.full_model )
    RCDE_no_seq_comp[i] <- ((null_deviance_full_model - ll.residual.full_model) - (null_deviance_full_model - ll.residual.no_seq_comp_model))/(null_deviance_full_model - ll.residual.full_model )
    RCDE_no_basequality[i] <- ((null_deviance_full_model - ll.residual.full_model) - (null_deviance_full_model - ll.residual.no_basequality_model))/(null_deviance_full_model - ll.residual.full_model )
    RCDE_no_overlap[i] <- ((null_deviance_full_model - ll.residual.full_model) - (null_deviance_full_model - ll.residual.no_overlap_model))/(null_deviance_full_model - ll.residual.full_model )
    RCDE_no_length[i] <- ((null_deviance_full_model - ll.residual.full_model) - (null_deviance_full_model - ll.residual.no_length_model))/(null_deviance_full_model - ll.residual.full_model )

    BICs[[i]] <-BIC(full_poisson_model, no_source_poisson_model, no_seq_comp_poisson_model, no_basequality_poisson_model, no_proportion_homopolymers_poisson_model, no_overlap_poisson_model, no_length_poisson_model)

  } 

  data.frame(motif_list, 
             deviance_expl_full_model ,
             deviance_expl_no_source_model ,
             RCDE_no_source ,
             deviance_expl_no_seq_comp_model, 
             RCDE_no_seq_comp,
             deviance_expl_no_basequality_model, 
             RCDE_no_basequality,
             deviance_expl_no_overlap_model ,
             RCDE_no_overlap, 
             deviance_expl_no_length_model ,
             RCDE_no_length ,
             deviance_expl_no_proportion_homopolymers_model ,
             RCDE_no_proportion_homopolymers 
  ) %>% 
    write.table(., file = "Illumina_smm_poisson_regression_model_results.txt", quote = FALSE, row.names = FALSE, sep = "\t")


  data.frame(t(data.frame(matrix(unlist(coefficients), nrow=length(coefficients), byrow=TRUE)) )) -> coefficients_output
  names(coefficients_output) <- motif_list
  coefficients_names<-row.names(data.frame(summary(full_poisson_model)$coefficients))
  coefficients_output %>% cbind(coefficients_names, .) %>%
    write.table(., file = "Illumina_smm_poisson_regression_model_coefficients.txt", quote = FALSE, row.names = FALSE, sep = "\t")

  data.frame(t(data.frame(matrix(unlist(coefficients_pvalues), nrow=length(coefficients_pvalues), byrow=TRUE)) )) -> coefficients_pvalues_output
  names(coefficients_pvalues_output) <- motif_list
  coefficients_names<-row.names(data.frame(summary(full_poisson_model)$coefficients))
  coefficients_pvalues_output %>% cbind(coefficients_names, .) %>%
    write.table(., file = "Illumina_smm_poisson_regression_model_coefficients_pvalues.txt", quote = FALSE, row.names = FALSE, sep = "\t")



################################################### 
#### Illumina Insertion mismatches - Generalized linear model with poisson model

  deviance_expl_full_model <- c()
  deviance_expl_no_source_model <- c()
  deviance_expl_no_proportion_homopolymers_model <- c()
  deviance_expl_no_seq_comp_model<- c()
  deviance_expl_no_basequality_model<- c() 
  deviance_expl_no_overlap_model<- c() 
  deviance_expl_no_length_model<- c() 

  RCDE_no_source <- c()
  RCDE_no_proportion_homopolymers <- c() 
  RCDE_no_seq_comp <- c()
  RCDE_no_basequality <- c()
  RCDE_no_overlap <- c()
  RCDE_no_length <- c()

  BICs<-list()
  theta <- c()
  coefficients <- list()
  coefficients_pvalues <- list()
  all_poisson_models<-list()
  for (i in 1:length(motif_list)){
    print(motif_list[i])

    print("full")
    Illumina_curated %>% 
      filter(type == paste0(motif_list[i])) %>%
      na.omit() %>%
      glm(ins ~ factor(source) + 
            offset(log(depth)) +
            factor(overlap_with_other_nonB) + 
        factor(overlap_with_other_nonB)*log(length) + 
        log(length) + 
        factor(source)*log(length) + 
            V1 + V2 + V3 + basequality + 
            factor(proportion_homopolymers_binary),  
          data = .,
          family = poisson(link = "log")) -> full_poisson_model

      data.frame(cooks.distance(full_poisson_model)) %>%
      mutate(row_number = 1:dplyr::n()) %>%
      filter(cooks.distance.full_poisson_model. > 1) %>%
      dplyr::select(row_number) %>%
      unlist() %>%
      as.vector() -> influential_cases

      print("full")
    Illumina_curated %>% 
      filter(type == paste0(motif_list[i])) %>%
      na.omit() %>%
      mutate(row_number = 1:dplyr::n()) %>%
      filter(row_number %not_in% influential_cases) %>%
      glm(ins ~ factor(source) + 
            offset(log(depth)) +
            factor(overlap_with_other_nonB) + 
        factor(overlap_with_other_nonB)*log(length) + 
        log(length) + 
        factor(source)*log(length) + 
            V1 + V2 + V3 + basequality + 
            factor(proportion_homopolymers_binary),  
          data = .,
          family = poisson(link = "log")) -> full_poisson_model

    coefficients[[i]] <- data.frame(summary(full_poisson_model)$coefficients) %>% dplyr::select(Estimate) %>% unlist() %>% as.vector()
    coefficients_pvalues[[i]] <- data.frame(summary(full_poisson_model)$coefficients[,4])

    print("no_factor(source)")
    Illumina_curated %>% 
      filter(type == paste0(motif_list[i])) %>% 
      na.omit() %>%
      mutate(row_number = 1:dplyr::n()) %>%
      filter(row_number %not_in% influential_cases) %>%
      glm(ins ~ 
            offset(log(depth)) +
        log(length) + 
            V1 + V2 + V3 + basequality + 
            factor(proportion_homopolymers_binary),  
          data = .,
          family = poisson(link = "log")) -> no_source_poisson_model


    print("no_proportion_homopolymers") 
    Illumina_curated %>% 
      filter(type == paste0(motif_list[i])) %>% 
      na.omit() %>%
      mutate(row_number = 1:dplyr::n()) %>%
      filter(row_number %not_in% influential_cases) %>%
      glm(ins ~ factor(source) + 
            offset(log(depth)) +
            factor(overlap_with_other_nonB) + 
        factor(overlap_with_other_nonB)*log(length) + 
        log(length) + 
        factor(source)*log(length) + 
            V1 + V2 + V3 + basequality , 
          data = ., 
          family = poisson(link = "log"))-> no_proportion_homopolymers_poisson_model

    print("no_seq_comp")
    Illumina_curated %>% 
      filter(type == paste0(motif_list[i])) %>% 
      na.omit() %>%
      mutate(row_number = 1:dplyr::n()) %>%
      filter(row_number %not_in% influential_cases) %>%
      glm(ins ~ factor(source) + 
            offset(log(depth)) +
            factor(overlap_with_other_nonB) + 
        factor(overlap_with_other_nonB)*log(length) + 
        log(length) + 
        factor(source)*log(length) + 
        basequality + 
            factor(proportion_homopolymers_binary),  
          data = ., 
          family = poisson(link = "log")) -> no_seq_comp_poisson_model

    print("no_basequality")
    Illumina_curated %>% 
      filter(type == paste0(motif_list[i])) %>% 
      na.omit() %>%
      mutate(row_number = 1:dplyr::n()) %>%
      filter(row_number %not_in% influential_cases) %>%
      glm(ins ~ factor(source) + 
            offset(log(depth)) +
            factor(overlap_with_other_nonB) + 
        factor(overlap_with_other_nonB)*log(length) + 
        log(length) + 
        factor(source)*log(length) + 
            V1 + V2 + V3 + 
            factor(proportion_homopolymers_binary),  
          data = ., 
          family = poisson(link = "log"))-> no_basequality_poisson_model

    Illumina_curated %>% 
      filter(type == paste0(motif_list[i])) %>% 
      na.omit() %>%
      mutate(row_number = 1:dplyr::n()) %>%
      filter(row_number %not_in% influential_cases) %>%
      glm(ins ~ factor(source) + 
            offset(log(depth)) +
        log(length) + 
        factor(source)*log(length) + 
            V1 + V2 + V3 + basequality + 
            factor(proportion_homopolymers_binary),  
          data = .,
          family = poisson(link = "log")) -> no_overlap_poisson_model

    Illumina_curated %>% 
      filter(type == paste0(motif_list[i])) %>% 
      na.omit() %>%
      mutate(row_number = 1:dplyr::n()) %>%
      filter(row_number %not_in% influential_cases) %>%
      glm(ins ~ factor(source) + 
            offset(log(depth)) +
            factor(overlap_with_other_nonB) + 
            V1 + V2 + V3 + basequality + 
            factor(proportion_homopolymers_binary),  
          data = .,
          family = poisson(link = "log")) -> no_length_poisson_model


    poisson_models_per_motif<-list(full_poisson_model, 
                                   no_source_poisson_model,
                                   no_seq_comp_poisson_model,
                                   no_basequality_poisson_model,
                                   no_proportion_homopolymers_poisson_model,
                       no_overlap_poisson_model,
                                   no_length_poisson_model) 

    all_poisson_models[[i]]<-poisson_models_per_motif 

    null_deviance_full_model <- full_poisson_model$null.deviance
    ll.null.full_model <- full_poisson_model$null.deviance
    ll.residual.full_model <- full_poisson_model$deviance

    ll.residual.no_source_model <- no_source_poisson_model$deviance
    ll.null.no_source_model <- no_source_poisson_model$null.deviance

    ll.residual.no_proportion_homopolymers_model <- no_proportion_homopolymers_poisson_model$deviance
    ll.null.no_proportion_homopolymers_model <- no_proportion_homopolymers_poisson_model$null.deviance

    ll.residual.no_seq_comp_model <- no_seq_comp_poisson_model$deviance
    ll.null.no_seq_comp_model <- no_seq_comp_poisson_model$null.deviance

    ll.residual.no_basequality_model <- no_basequality_poisson_model$deviance
    ll.null.no_basequality_model <- no_basequality_poisson_model$null.deviance

    ll.residual.no_overlap_model <- no_overlap_poisson_model$deviance
    ll.null.no_overlap_model <- no_overlap_poisson_model$null.deviance

    ll.residual.no_length_model <- no_length_poisson_model$deviance
    ll.null.no_length_model <- no_length_poisson_model$null.deviance


    # Deviance explained by the full model:
    deviance_expl_full_model[i] <- 1-ll.residual.full_model/ll.null.full_model
    deviance_expl_no_source_model[i]<- 1-ll.residual.no_source_model/ll.null.no_source_model
    deviance_expl_no_proportion_homopolymers_model[i]<- 1-ll.residual.no_proportion_homopolymers_model/ll.null.no_proportion_homopolymers_model
    deviance_expl_no_seq_comp_model[i]<- 1-ll.residual.no_seq_comp_model/ll.null.no_seq_comp_model
    deviance_expl_no_basequality_model[i]<- 1-ll.residual.no_basequality_model/ll.null.no_basequality_model
    deviance_expl_no_overlap_model[i]<- 1-ll.residual.no_overlap_model/ll.null.no_overlap_model
    deviance_expl_no_length_model[i]<- 1-ll.residual.no_length_model/ll.null.no_length_model

    RCDE_no_source[i] <- ((null_deviance_full_model - ll.residual.full_model) - (null_deviance_full_model - ll.residual.no_source_model))/(null_deviance_full_model - ll.residual.full_model )
    RCDE_no_proportion_homopolymers[i] <- ((null_deviance_full_model - ll.residual.full_model) - (null_deviance_full_model - ll.residual.no_proportion_homopolymers_model))/(null_deviance_full_model - ll.residual.full_model )
    RCDE_no_seq_comp[i] <- ((null_deviance_full_model - ll.residual.full_model) - (null_deviance_full_model - ll.residual.no_seq_comp_model))/(null_deviance_full_model - ll.residual.full_model )
    RCDE_no_basequality[i] <- ((null_deviance_full_model - ll.residual.full_model) - (null_deviance_full_model - ll.residual.no_basequality_model))/(null_deviance_full_model - ll.residual.full_model )
    RCDE_no_overlap[i] <- ((null_deviance_full_model - ll.residual.full_model) - (null_deviance_full_model - ll.residual.no_overlap_model))/(null_deviance_full_model - ll.residual.full_model )
    RCDE_no_length[i] <- ((null_deviance_full_model - ll.residual.full_model) - (null_deviance_full_model - ll.residual.no_length_model))/(null_deviance_full_model - ll.residual.full_model )

    BICs[[i]] <-BIC(full_poisson_model, no_source_poisson_model, no_seq_comp_poisson_model, no_basequality_poisson_model, no_proportion_homopolymers_poisson_model, no_overlap_poisson_model, no_length_poisson_model)

  } 

  data.frame(motif_list, 
             deviance_expl_full_model ,
             deviance_expl_no_source_model ,
             RCDE_no_source ,
             deviance_expl_no_seq_comp_model, 
             RCDE_no_seq_comp,
             deviance_expl_no_basequality_model, 
             RCDE_no_basequality,
             deviance_expl_no_overlap_model ,
             RCDE_no_overlap, 
             deviance_expl_no_length_model ,
             RCDE_no_length ,
             deviance_expl_no_proportion_homopolymers_model ,
             RCDE_no_proportion_homopolymers 
  ) %>% 
    write.table(., file = "Illumina_ins_poisson_regression_model_results.txt", quote = FALSE, row.names = FALSE, sep = "\t")


  data.frame(t(data.frame(matrix(unlist(coefficients), nrow=length(coefficients), byrow=TRUE)) )) -> coefficients_output
  names(coefficients_output) <- motif_list
  coefficients_names<-row.names(data.frame(summary(full_poisson_model)$coefficients))
  coefficients_output %>% cbind(coefficients_names, .) %>%
    write.table(., file = "Illumina_ins_poisson_regression_model_coefficients.txt", quote = FALSE, row.names = FALSE, sep = "\t")

  data.frame(t(data.frame(matrix(unlist(coefficients_pvalues), nrow=length(coefficients_pvalues), byrow=TRUE)) )) -> coefficients_pvalues_output
  names(coefficients_pvalues_output) <- motif_list
  coefficients_names<-row.names(data.frame(summary(full_poisson_model)$coefficients))
  coefficients_pvalues_output %>% cbind(coefficients_names, .) %>%
    write.table(., file = "Illumina_ins_poisson_regression_model_coefficients_pvalues.txt", quote = FALSE, row.names = FALSE, sep = "\t")



################################################### 
#### Illumina Deletion mismatches - Generalized linear model with poisson model
#######

  deviance_expl_full_model <- c()
  deviance_expl_no_source_model <- c()
  deviance_expl_no_proportion_homopolymers_model <- c()
  deviance_expl_no_seq_comp_model<- c()
  deviance_expl_no_basequality_model<- c() 
  deviance_expl_no_overlap_model<- c() 
  deviance_expl_no_length_model<- c() 

  RCDE_no_source <- c()
  RCDE_no_proportion_homopolymers <- c() 
  RCDE_no_seq_comp <- c()
  RCDE_no_basequality <- c()
  RCDE_no_overlap <- c()
  RCDE_no_length <- c()

  BICs<-list()
  theta <- c()
  coefficients <- list()
  coefficients_pvalues <- list()
  all_poisson_models<-list()
  for (i in 1:length(motif_list)){
    print(motif_list[i])

    print("full")
    Illumina_curated %>% 
      filter(type == paste0(motif_list[i])) %>%
      na.omit() %>%
      glm(del ~ factor(source) + 
            offset(log(depth)) +
            factor(overlap_with_other_nonB) + 
        factor(overlap_with_other_nonB)*log(length) + 
        log(length) + 
        factor(source)*log(length) + 
            V1 + V2 + V3 + basequality + 
            factor(proportion_homopolymers_binary),  
          data = .,
          family = poisson(link = "log")) -> full_poisson_model

      data.frame(cooks.distance(full_poisson_model)) %>%
      mutate(row_number = 1:dplyr::n()) %>%
      filter(cooks.distance.full_poisson_model. > 1) %>%
      dplyr::select(row_number) %>%
      unlist() %>%
      as.vector() -> influential_cases

      print("full")
    Illumina_curated %>% 
      filter(type == paste0(motif_list[i])) %>%
      na.omit() %>%
      mutate(row_number = 1:dplyr::n()) %>%
      filter(row_number %not_in% influential_cases) %>%
      glm(del ~ factor(source) + 
            offset(log(depth)) +
            factor(overlap_with_other_nonB) + 
        factor(overlap_with_other_nonB)*log(length) + 
        log(length) + 
        factor(source)*log(length) + 
            V1 + V2 + V3 + basequality + 
            factor(proportion_homopolymers_binary),  
          data = .,
          family = poisson(link = "log")) -> full_poisson_model

    coefficients[[i]] <- data.frame(summary(full_poisson_model)$coefficients) %>% dplyr::select(Estimate) %>% unlist() %>% as.vector()
    coefficients_pvalues[[i]] <- data.frame(summary(full_poisson_model)$coefficients[,4])

    print("no_factor(source)")
    Illumina_curated %>% 
      filter(type == paste0(motif_list[i])) %>% 
      na.omit() %>%
      mutate(row_number = 1:dplyr::n()) %>%
      filter(row_number %not_in% influential_cases) %>%
      glm(del ~ 
            offset(log(depth)) +
        log(length) + 
            V1 + V2 + V3 + basequality + 
            factor(proportion_homopolymers_binary),  
          data = .,
          family = poisson(link = "log")) -> no_source_poisson_model


    print("no_proportion_homopolymers") 
    Illumina_curated %>% 
      filter(type == paste0(motif_list[i])) %>% 
      na.omit() %>%
      mutate(row_number = 1:dplyr::n()) %>%
      filter(row_number %not_in% influential_cases) %>%
      glm(del ~ factor(source) + 
            offset(log(depth)) +
            factor(overlap_with_other_nonB) + 
        factor(overlap_with_other_nonB)*log(length) + 
        log(length) + 
        factor(source)*log(length) + 
            V1 + V2 + V3 + basequality , 
          data = ., 
          family = poisson(link = "log"))-> no_proportion_homopolymers_poisson_model

    print("no_seq_comp")
    Illumina_curated %>% 
      filter(type == paste0(motif_list[i])) %>% 
      na.omit() %>%
      mutate(row_number = 1:dplyr::n()) %>%
      filter(row_number %not_in% influential_cases) %>%
      glm(del ~ factor(source) + 
            offset(log(depth)) +
            factor(overlap_with_other_nonB) + 
        factor(overlap_with_other_nonB)*log(length) + 
        log(length) + 
        factor(source)*log(length) + 
        basequality + 
            factor(proportion_homopolymers_binary),  
          data = ., 
          family = poisson(link = "log")) -> no_seq_comp_poisson_model

    print("no_basequality")
    Illumina_curated %>% 
      filter(type == paste0(motif_list[i])) %>% 
      na.omit() %>%
      mutate(row_number = 1:dplyr::n()) %>%
      filter(row_number %not_in% influential_cases) %>%
      glm(del ~ factor(source) + 
            offset(log(depth)) +
            factor(overlap_with_other_nonB) + 
        factor(overlap_with_other_nonB)*log(length) + 
        log(length) + 
        factor(source)*log(length) + 
            V1 + V2 + V3 + 
            factor(proportion_homopolymers_binary),  
          data = ., 
          family = poisson(link = "log"))-> no_basequality_poisson_model

    Illumina_curated %>% 
      filter(type == paste0(motif_list[i])) %>% 
      na.omit() %>%
      mutate(row_number = 1:dplyr::n()) %>%
      filter(row_number %not_in% influential_cases) %>%
      glm(del ~ factor(source) + 
            offset(log(depth)) +
        log(length) + 
        factor(source)*log(length) + 
            V1 + V2 + V3 + basequality + 
            factor(proportion_homopolymers_binary),  
          data = .,
          family = poisson(link = "log")) -> no_overlap_poisson_model

    Illumina_curated %>% 
      filter(type == paste0(motif_list[i])) %>% 
      na.omit() %>%
      mutate(row_number = 1:dplyr::n()) %>%
      filter(row_number %not_in% influential_cases) %>%
      glm(del ~ factor(source) + 
            offset(log(depth)) +
            factor(overlap_with_other_nonB) + 
            V1 + V2 + V3 + basequality + 
            factor(proportion_homopolymers_binary),  
          data = .,
          family = poisson(link = "log")) -> no_length_poisson_model


    poisson_models_per_motif<-list(full_poisson_model, 
                                   no_source_poisson_model,
                                   no_seq_comp_poisson_model,
                                   no_basequality_poisson_model,
                                   no_proportion_homopolymers_poisson_model,
                       no_overlap_poisson_model,
                                   no_length_poisson_model) 

    all_poisson_models[[i]]<-poisson_models_per_motif 

    null_deviance_full_model <- full_poisson_model$null.deviance
    ll.null.full_model <- full_poisson_model$null.deviance
    ll.residual.full_model <- full_poisson_model$deviance

    ll.residual.no_source_model <- no_source_poisson_model$deviance
    ll.null.no_source_model <- no_source_poisson_model$null.deviance

    ll.residual.no_proportion_homopolymers_model <- no_proportion_homopolymers_poisson_model$deviance
    ll.null.no_proportion_homopolymers_model <- no_proportion_homopolymers_poisson_model$null.deviance

    ll.residual.no_seq_comp_model <- no_seq_comp_poisson_model$deviance
    ll.null.no_seq_comp_model <- no_seq_comp_poisson_model$null.deviance

    ll.residual.no_basequality_model <- no_basequality_poisson_model$deviance
    ll.null.no_basequality_model <- no_basequality_poisson_model$null.deviance

    ll.residual.no_overlap_model <- no_overlap_poisson_model$deviance
    ll.null.no_overlap_model <- no_overlap_poisson_model$null.deviance

    ll.residual.no_length_model <- no_length_poisson_model$deviance
    ll.null.no_length_model <- no_length_poisson_model$null.deviance


    # Deviance explained by the full model:
    deviance_expl_full_model[i] <- 1-ll.residual.full_model/ll.null.full_model
    deviance_expl_no_source_model[i]<- 1-ll.residual.no_source_model/ll.null.no_source_model
    deviance_expl_no_proportion_homopolymers_model[i]<- 1-ll.residual.no_proportion_homopolymers_model/ll.null.no_proportion_homopolymers_model
    deviance_expl_no_seq_comp_model[i]<- 1-ll.residual.no_seq_comp_model/ll.null.no_seq_comp_model
    deviance_expl_no_basequality_model[i]<- 1-ll.residual.no_basequality_model/ll.null.no_basequality_model
    deviance_expl_no_overlap_model[i]<- 1-ll.residual.no_overlap_model/ll.null.no_overlap_model
    deviance_expl_no_length_model[i]<- 1-ll.residual.no_length_model/ll.null.no_length_model

    RCDE_no_source[i] <- ((null_deviance_full_model - ll.residual.full_model) - (null_deviance_full_model - ll.residual.no_source_model))/(null_deviance_full_model - ll.residual.full_model )
    RCDE_no_proportion_homopolymers[i] <- ((null_deviance_full_model - ll.residual.full_model) - (null_deviance_full_model - ll.residual.no_proportion_homopolymers_model))/(null_deviance_full_model - ll.residual.full_model )
    RCDE_no_seq_comp[i] <- ((null_deviance_full_model - ll.residual.full_model) - (null_deviance_full_model - ll.residual.no_seq_comp_model))/(null_deviance_full_model - ll.residual.full_model )
    RCDE_no_basequality[i] <- ((null_deviance_full_model - ll.residual.full_model) - (null_deviance_full_model - ll.residual.no_basequality_model))/(null_deviance_full_model - ll.residual.full_model )
    RCDE_no_overlap[i] <- ((null_deviance_full_model - ll.residual.full_model) - (null_deviance_full_model - ll.residual.no_overlap_model))/(null_deviance_full_model - ll.residual.full_model )
    RCDE_no_length[i] <- ((null_deviance_full_model - ll.residual.full_model) - (null_deviance_full_model - ll.residual.no_length_model))/(null_deviance_full_model - ll.residual.full_model )

    BICs[[i]] <-BIC(full_poisson_model, no_source_poisson_model, no_seq_comp_poisson_model, no_basequality_poisson_model, no_proportion_homopolymers_poisson_model, no_overlap_poisson_model, no_length_poisson_model)

  } 

  data.frame(motif_list, 
             deviance_expl_full_model ,
             deviance_expl_no_source_model ,
             RCDE_no_source ,
             deviance_expl_no_seq_comp_model, 
             RCDE_no_seq_comp,
             deviance_expl_no_basequality_model, 
             RCDE_no_basequality,
             deviance_expl_no_overlap_model ,
             RCDE_no_overlap, 
             deviance_expl_no_length_model ,
             RCDE_no_length ,
             deviance_expl_no_proportion_homopolymers_model ,
             RCDE_no_proportion_homopolymers 
  ) %>% 
    write.table(., file = "Illumina_del_poisson_regression_model_results.txt", quote = FALSE, row.names = FALSE, sep = "\t")


  data.frame(t(data.frame(matrix(unlist(coefficients), nrow=length(coefficients), byrow=TRUE)) )) -> coefficients_output
  names(coefficients_output) <- motif_list
  coefficients_names<-row.names(data.frame(summary(full_poisson_model)$coefficients))
  coefficients_output %>% cbind(coefficients_names, .) %>%
    write.table(., file = "Illumina_del_poisson_regression_model_coefficients.txt", quote = FALSE, row.names = FALSE, sep = "\t")

  data.frame(t(data.frame(matrix(unlist(coefficients_pvalues), nrow=length(coefficients_pvalues), byrow=TRUE)) )) -> coefficients_pvalues_output
  names(coefficients_pvalues_output) <- motif_list
  coefficients_names<-row.names(data.frame(summary(full_poisson_model)$coefficients))
  coefficients_pvalues_output %>% cbind(coefficients_names, .) %>%
    write.table(., file = "Illumina_del_poisson_regression_model_coefficients_pvalues.txt", quote = FALSE, row.names = FALSE, sep = "\t")
```

### Producing figures

Input for these scripts are data tables derived from the error calling script described above.

```r


```

## 
