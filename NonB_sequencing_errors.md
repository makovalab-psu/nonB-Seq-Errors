# Non-B DNA associated errors

## General notes

This is a collection of scripts that has been used to curate data, perform analysis, and generate figures for Weissensteiner et al. 2022. The overall workflow is as follows: After mapping sequencing read data (in case of HiFi and ONT) to the human reference genome hg19, or retrieving alignment (bam) files from the Genome in a Bottle website ([giab_data_indexes/AshkenazimTrio at master · genome-in-a-bottle/giab_data_indexes · GitHub](https://github.com/genome-in-a-bottle/giab_data_indexes/tree/master/AshkenazimTrio)). Then non-B-DNA forming motifs are either retrieved from downloaded from [Advanced Biomedical Computing Center (ABCC) | non-B DB | Home](https://nonb-abcc.ncifcrf.gov/apps/site/default) or, in case of G4 motifs, generated. Following that motifs are curated, and sequence-specific features (nucleotide composition, mappability, etc.) are retrieved for each one, and a random set of controls is generated that matches motifs in size and number. Then, measures of sequencing success (single-nucleotide, insertion, and deletion mismatches, read depth and base quality), are determined for every motif and control. For every motif and control sequence, we assigned a unique ID, which consists of the non-B type, source (motif / control), chromosome, and start of the bed entry (e.g., G4Motifs_motif_1_123456)). In the downstream analysis, we used Poisson regression models to determine the effect of non-B motifs on the variability of sequencing error rates

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
Sys.getenv("YOUR_VAR")
library(easypackages)
libraries("ggplot2","data.table", "patchwork", "car", "cowplot", "scales","patchwork", "plyr", "tidyr","msme", "dplyr","vroom", "glmmADMB", "purrr", "grid", "gridExtra","cowplot", "lme4", "MatchIt", "compositions")
`%not_in%` <- purrr::negate(`%in%`)
setwd("/Users/matthiasweissensteiner/Dropbox/PostDoc/proj/Weissensteiner_et_al_2022_NonB_error_project/prelim_results/CellGenomics_submission_figures_script")
library(vroom)
library(scales)
library(ggplot2) 
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(patchwork)
`%not_in%` <- purrr::negate(`%in%`)

homopolymers <- vroom("homopolymers_annotated.bed", col_names = FALSE)
data.frame(homopolymers) ->  homopolymers
names(homopolymers) <- c("chr", "start", "end", "number_homopolymers", "proportion_homopolymers")
motif_list<-c("APhasedRepeats", "DirectRepeats", "G4Motifs", "InvertedRepeats", "MirrorRepeats", "ZDNAMotifs")

mappability_36mer <- vroom("mappability_36mer_output_parsed", col_names = FALSE)
data.frame(mappability_36mer) -> mappability_36mer
names(mappability_36mer) <- c("id", "mappability_36mer")

# Set colors for motif and control:
colors<-c("#2c7bb6", "#d7191c")
n <- c("control", "motif")
names(colors) <- n

base_size <- 9

colors_simple_section <- c("#1b9e77", "#7570b3")
n<- c("spacer_loop", "repeat_stem")
names(colors_simple_section) <- n

labeller_rcde_plot <- labeller(motif_list = c("APhasedRepeats" = "APR", "DirectRepeats"="DR",
                                              "G4Motifs" = "G4", "InvertedRepeats"="IR", 
                                              "MirrorRepeats" = "MR", "ZDNAMotifs" = "ZDNA"))

short_type_names <- c("APR", "DR", "G4", "IR", "MR", "ZDNA")


# Boxplot function:
f <- function(x) {
  r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.9))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

point <- format_format(big.mark = " ", decimal.mark = ",", scientific = FALSE)

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

# Illumina
#####

# Moderate filter
Illumina_curated <- vroom("Illumina_curated_fully_filtered", col_names  = T)
data.frame(Illumina_curated) -> Illumina_curated
Illumina_curated$type <- as.factor(Illumina_curated$type)
Illumina_curated$source <- as.factor(Illumina_curated$source)
Illumina_curated$source <- with(Illumina_curated, relevel(source, "motif"))

# Stringent filter
Illumina_ultrastringent_filter <- read.table("Illumina_ultrastringent_filter", header = T)
Illumina_ultrastringent_filter$type <- as.factor(Illumina_ultrastringent_filter$type)
Illumina_ultrastringent_filter$source <- as.factor(Illumina_ultrastringent_filter$source)
Illumina_ultrastringent_filter$source <- with(Illumina_ultrastringent_filter, relevel(source, "motif"))

# Motif subsections - moderate filter
Illumina<- vroom("Illumina.fwd.motif_subsections.results", col_names = T)
data.frame(Illumina) -> Illumina_motif_subsections
head(Illumina_motif_subsections)

names(Illumina_motif_subsections) <- 
  c("chr", "start" , "end", "smm", "ins" , "del", "depth", "type", "length", 
    "section", "id")
Illumina_motif_subsections$section <- factor(Illumina_motif_subsections$section)

Illumina_motif_subsections %>%
  mutate(smm_rate = smm / depth,
         ins_rate = ins / depth,
         del_rate = del / depth) -> Illumina_motif_subsections

Illumina_motif_subsections %>%
  mutate(simple_section = ifelse(grepl("Tract", section), "repeat_stem",
                                 ifelse(grepl("Spacer", section), "spacer_loop", 
                                        ifelse(grepl("Repeat", section), "repeat_stem",
                                               ifelse(grepl("stem", section), "repeat_stem", "spacer_loop"))))) -> Illumina_motif_subsections

# Motif subsections - stringent filter
Illumina_ultrastringent_filter_ids <- read.table("Illumina_ultrastringent_filter_ids", header = T)
Illumina_ultrastringent_filter_ids %>% 
  select(id) %>% 
  unlist() %>% as.vector() -> Illumina_ultrastringent_filter_ids

Illumina_motif_subsections %>%
  filter(id %in% Illumina_ultrastringent_filter_ids) -> Illumina_motif_subsections_ultrastringent
#####

# HiFi
#####
# Moderate filter
HiFi_curated <- vroom("HiFi_curated_fully_filtered", col_names = T)
data.frame(HiFi_curated) -> HiFi_curated
HiFi_curated$type <- as.factor(HiFi_curated$type)
HiFi_curated$source <- as.factor(HiFi_curated$source)
HiFi_curated$source <- with(HiFi_curated, relevel(source, "motif"))

# Stringent filter
HiFi_ultrastringent_filter <- read.table("HiFi_ultrastringent_filter", header = T)
HiFi_ultrastringent_filter$type <- as.factor(HiFi_ultrastringent_filter$type)
HiFi_ultrastringent_filter$source <- as.factor(HiFi_ultrastringent_filter$source)
HiFi_ultrastringent_filter$source <- with(HiFi_ultrastringent_filter, relevel(source, "motif"))

# Motif subsections - moderate filter
HiFi<- vroom("HiFi.fwd.motif_subsections.results", col_names = T)
data.frame(HiFi) -> HiFi_motif_subsections
head(HiFi_motif_subsections)

names(HiFi_motif_subsections) <- 
  c("chr", "start" , "end", "smm", "ins" , "del", "depth", "type", "length", 
    "section", "id")
HiFi_motif_subsections$section <- factor(HiFi_motif_subsections$section)

HiFi_motif_subsections %>%
  mutate(smm_rate = smm / depth,
         ins_rate = ins / depth,
         del_rate = del / depth) -> HiFi_motif_subsections

HiFi_motif_subsections %>%
  mutate(simple_section = ifelse(grepl("Tract", section), "repeat_stem",
                                 ifelse(grepl("Spacer", section), "spacer_loop", 
                                        ifelse(grepl("Repeat", section), "repeat_stem",
                                               ifelse(grepl("stem", section), "repeat_stem", "spacer_loop"))))) -> HiFi_motif_subsections

# Motif subsections - stringent filter
HiFi_ultrastringent_filter_ids <- read.table("HiFi_ultrastringent_filter_ids", header = T)
HiFi_ultrastringent_filter_ids %>% 
  select(id) %>% 
  unlist() %>% as.vector() -> HiFi_ultrastringent_filter_ids

HiFi_motif_subsections %>%
  filter(id %in% HiFi_ultrastringent_filter_ids) -> HiFi_motif_subsections_ultrastringent
#####

# ONT
#####
# Moderate filter
ONT<- vroom("ONT_curated_fully_filtered", col_names = T)
data.frame(ONT) -> ONT_curated
ONT_curated$type <- as.factor(ONT_curated$type)
ONT_curated$source <- as.factor(ONT_curated$source)
ONT_curated$source <- with(ONT_curated, relevel(source, "motif"))

# Stringent filter
ONT_ultrastringent_filter <- read.table("ONT_ultrastringent_filter", header = T)
ONT_ultrastringent_filter$type <- as.factor(ONT_ultrastringent_filter$type)
ONT_ultrastringent_filter$source <- as.factor(ONT_ultrastringent_filter$source)
ONT_ultrastringent_filter$source <- with(ONT_ultrastringent_filter, relevel(source, "motif"))

# Motif subsections - moderate filter
ONT<- vroom("ONT.fwd.motif_subsections.results", col_names = T)
data.frame(ONT) -> ONT_motif_subsections
head(ONT_motif_subsections)

names(ONT_motif_subsections) <- 
  c("chr", "start" , "end", "smm", "ins" , "del", "depth", "type", "length", 
    "section", "id")
ONT_motif_subsections$section <- factor(ONT_motif_subsections$section)

ONT_motif_subsections %>%
  mutate(smm_rate = smm / depth,
         ins_rate = ins / depth,
         del_rate = del / depth) -> ONT_motif_subsections

ONT_motif_subsections %>%
  mutate(simple_section = ifelse(grepl("Tract", section), "repeat_stem",
                                 ifelse(grepl("Spacer", section), "spacer_loop", 
                                        ifelse(grepl("Repeat", section), "repeat_stem",
                                               ifelse(grepl("stem", section), "repeat_stem", "spacer_loop"))))) -> ONT_motif_subsections

# Motif subsections - stringent filter
ONT_ultrastringent_filter_ids <- read.table("ONT_ultrastringent_filter_ids", header = T)
ONT_ultrastringent_filter_ids %>% 
  select(id) %>% 
  unlist() %>% as.vector() -> ONT_ultrastringent_filter_ids

ONT_motif_subsections %>%
  filter(id %in% ONT_ultrastringent_filter_ids) -> ONT_motif_subsections_ultrastringent

#####

# motif summary:

Illumina_curated %>% filter(source == "motif") %>% 
  group_by(type) %>% 
  tally() -> tally_table
Illumina_curated %>% 
  filter(source == "motif") %>% 
  group_by(type) %>% 
  summarize(total_length = sum(length), 
            mean_length = mean(length) , 
            median_length = median(length)) -> length_table
tally_table %>% 
  select(-type) %>% 
  cbind.data.frame(length_table, .) %>%  
  write.table(., file = "moderate_filter_motif_table.txt", quote = FALSE, row.names = FALSE, sep = "\t")


#################################################
# Raw error rates
#################################################

#Overall
#####
Illumina_curated %>%
  summarize(sum_smm = sum(smm),
            sum_ins = sum(ins),
            sum_del = sum(del),
            per_motif_smm = mean(smm_rate, na.rm=T), 
            per_motif_ins = mean(ins_rate, na.rm=T), 
            per_motif_del = mean(del_rate, na.rm=T), 
            sum_depth = sum(depth),
            raw_smm_rate = sum(smm)/sum(depth),
            raw_ins_rate = sum(ins)/sum(depth),
            raw_del_rate = sum(del)/sum(depth))  %>% 
  mutate(technology = "Illumina",
         filter_level = "moderate") -> Illumina_overall_error_rates_moderate

HiFi_curated %>%
  summarize(sum_smm = sum(smm),
            sum_ins = sum(ins),
            sum_del = sum(del),
            per_motif_smm = mean(smm_rate, na.rm = T), 
            per_motif_ins = mean(ins_rate, na.rm = T), 
            per_motif_del = mean(del_rate, na.rm = T), 
            sum_depth = sum(depth),
            raw_smm_rate = sum(smm)/sum(depth),
            raw_ins_rate = sum(ins)/sum(depth),
            raw_del_rate = sum(del)/sum(depth)) %>% 
  mutate(technology = "HiFi",
         filter_level = "moderate") -> HiFi_overall_error_rates_moderate


ONT_curated %>%
  summarize(sum_smm = sum(smm),
            sum_ins = sum(ins),
            sum_del = sum(del),
            per_motif_smm = mean(smm_rate, na.rm = T), 
            per_motif_ins = mean(ins_rate, na.rm = T), 
            per_motif_del = mean(del_rate, na.rm = T), 
            sum_depth = sum(depth),
            raw_smm_rate = sum(smm)/sum(depth),
            raw_ins_rate = sum(ins)/sum(depth),
            raw_del_rate = sum(del)/sum(depth)) %>% 
  mutate(technology = "ONT",
         filter_level = "moderate") -> ONT_overall_error_rates_moderate

rbind.data.frame(
  Illumina_overall_error_rates_moderate,
  HiFi_overall_error_rates_moderate,
  ONT_overall_error_rates_moderate
) %>%
  select(technology, filter_level, per_motif_smm, per_motif_ins, per_motif_del, raw_smm_rate, raw_ins_rate, raw_del_rate) %>% 
  melt() %>% 
  mutate(error_type = ifelse(grepl("smm", variable), "smm", ifelse(grepl("ins", variable), "ins", "del")),
         rate_type = ifelse(grepl("per", variable), "per_motif_mean", "aggregate_mean")) %>% 
  select(technology, filter_level, value, error_type, rate_type) %>% 
  dcast(., technology+filter_level+error_type~rate_type) -> overall_rates_moderate
#####

# Per motif and control
#####
Illumina_curated %>%
  group_by(type, source) %>%
  summarize(raw_smm_rate = sum(smm)/sum(depth),
            raw_ins_rate = sum(ins)/sum(depth),
            raw_del_rate = sum(del)/sum(depth),
            per_motif_smm_rate = mean(smm_rate, na.rm=T),
            per_motif_ins_rate = mean(ins_rate, na.rm = T), 
            per_motif_del_rate = mean(del_rate, na.rm = T),
            total_length = sum(length)
  ) %>%
  data.frame() %>%
  mutate(tech = "Illumina") -> raw_error_rates_Illumina_moderate

Illumina_ultrastringent_filter %>%
  group_by(type, source) %>%
  summarize(raw_smm_rate = sum(smm)/sum(depth),
            raw_ins_rate = sum(ins)/sum(depth),
            raw_del_rate = sum(del)/sum(depth),
            per_motif_smm_rate = mean(smm_rate, na.rm=T),
            per_motif_ins_rate = mean(ins_rate, na.rm = T), 
            per_motif_del_rate = mean(del_rate, na.rm = T),
            total_length = sum(length)
            ) %>%
  data.frame() %>%
  mutate(tech = "Illumina") -> raw_error_rates_Illumina_ultrastringent


  HiFi_curated %>%
  group_by(type, source) %>%
  summarize(raw_smm_rate = sum(smm)/sum(depth),
            raw_ins_rate = sum(ins)/sum(depth),
            raw_del_rate = sum(del)/sum(depth),
            per_motif_smm_rate = mean(smm_rate, na.rm=T),
            per_motif_ins_rate = mean(ins_rate, na.rm = T), 
            per_motif_del_rate = mean(del_rate, na.rm = T),
            total_length = sum(length)
  ) %>%
  data.frame() %>%
  mutate(tech = "HiFi") -> raw_error_rates_HiFi_moderate

HiFi_ultrastringent_filter %>%
  group_by(type, source) %>%
  summarize(raw_smm_rate = sum(smm)/sum(depth),
            raw_ins_rate = sum(ins)/sum(depth),
            raw_del_rate = sum(del)/sum(depth),
            per_motif_smm_rate = mean(smm_rate, na.rm=T),
            per_motif_ins_rate = mean(ins_rate, na.rm = T), 
            per_motif_del_rate = mean(del_rate, na.rm = T),
            total_length = sum(length)
  ) %>%
  data.frame() %>%
  mutate(tech = "HiFi") -> raw_error_rates_HiFi_ultrastringent

  ONT_curated %>%
  group_by(type, source) %>%
  summarize(raw_smm_rate = sum(smm)/sum(depth),
            raw_ins_rate = sum(ins)/sum(depth),
            raw_del_rate = sum(del)/sum(depth),
            per_motif_smm_rate = mean(smm_rate, na.rm=T),
            per_motif_ins_rate = mean(ins_rate, na.rm = T), 
            per_motif_del_rate = mean(del_rate, na.rm = T),
            total_length = sum(length)
  ) %>%
  data.frame() %>%
  mutate(tech = "ONT") -> raw_error_rates_ONT_moderate

ONT_ultrastringent_filter %>%
  group_by(type, source) %>%
  summarize(raw_smm_rate = sum(smm)/sum(depth),
            raw_ins_rate = sum(ins)/sum(depth),
            raw_del_rate = sum(del)/sum(depth),
            per_motif_smm_rate = mean(smm_rate, na.rm=T),
            per_motif_ins_rate = mean(ins_rate, na.rm = T), 
            per_motif_del_rate = mean(del_rate, na.rm = T),
            total_length = sum(length)
  ) %>%
  data.frame() %>%
  mutate(tech = "ONT") -> raw_error_rates_ONT_ultrastringent

rbind.data.frame(
  raw_error_rates_Illumina_moderate,
  raw_error_rates_HiFi_moderate,
  raw_error_rates_ONT_moderate) %>%
  mutate(filter_level = "moderate") -> error_rates_moderate

rbind.data.frame(
  raw_error_rates_Illumina_ultrastringent,
  raw_error_rates_HiFi_ultrastringent,
  raw_error_rates_ONT_ultrastringent) %>%
  mutate(filter_level = "stringent") -> error_rates_ultrastringent

#####

# Fold change table

# Moderate filter
#####
raw_error_rates_Illumina_moderate %>% 
  select(type, source, raw_smm_rate) %>% 
  dcast(., type~source) %>% 
  mutate(aggregate_smm_fold_change = motif/control) %>% 
  select(type, aggregate_smm_fold_change) %>% 
  melt() %>%
  mutate(tech = "Illumina") -> fold_change_table

raw_error_rates_Illumina_moderate %>% 
  select(type, source, raw_ins_rate) %>% 
  dcast(., type~source) %>% 
  mutate(aggregate_ins_fold_change = motif/control) %>% 
  select(type, aggregate_ins_fold_change) %>% 
  melt() %>% 
  mutate(tech = "Illumina") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_Illumina_moderate %>% 
  select(type, source, raw_del_rate) %>% 
  dcast(., type~source) %>% 
  mutate(aggregate_del_fold_change = motif/control) %>% 
  select(type, aggregate_del_fold_change) %>% 
  melt() %>% 
  mutate(tech = "Illumina") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_HiFi_moderate %>% 
  select(type, source, raw_smm_rate) %>% 
  dcast(., type~source) %>% 
  mutate(aggregate_smm_fold_change = motif/control) %>% 
  select(type, aggregate_smm_fold_change) %>% 
  melt() %>%
  mutate(tech = "HiFi") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_HiFi_moderate %>% 
  select(type, source, raw_ins_rate) %>% 
  dcast(., type~source) %>% 
  mutate(aggregate_ins_fold_change = motif/control) %>% 
  select(type, aggregate_ins_fold_change) %>% 
  melt() %>% 
  mutate(tech = "HiFi") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_HiFi_moderate %>% 
  select(type, source, raw_del_rate) %>% 
  dcast(., type~source) %>% 
  mutate(aggregate_del_fold_change = motif/control) %>% 
  select(type, aggregate_del_fold_change) %>% 
  melt() %>%
  mutate(tech = "HiFi") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_ONT_moderate %>% 
  select(type, source, raw_smm_rate) %>% 
  dcast(., type~source) %>% 
  mutate(aggregate_smm_fold_change = motif/control) %>% 
  select(type, aggregate_smm_fold_change) %>% 
  melt() %>%
  mutate(tech = "ONT") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_ONT_moderate %>% 
  select(type, source, raw_ins_rate) %>% 
  dcast(., type~source) %>% 
  mutate(aggregate_ins_fold_change = motif/control) %>% 
  select(type, aggregate_ins_fold_change) %>% 
  melt() %>% 
  mutate(tech = "ONT") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_ONT_moderate %>% 
  select(type, source, raw_del_rate) %>% 
  dcast(., type~source) %>% 
  mutate(aggregate_del_fold_change = motif/control) %>% 
  select(type, aggregate_del_fold_change) %>% 
  melt() %>%  
  mutate(tech = "ONT") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_Illumina_moderate %>% 
  select(type, source, per_motif_smm_rate) %>% 
  dcast(., type~source) %>% 
  mutate(per_motif_smm_fold_change = motif/control) %>% 
  select(type, per_motif_smm_fold_change) %>% 
  melt() %>%
  mutate(tech = "Illumina") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_Illumina_moderate %>% 
  select(type, source, per_motif_ins_rate) %>% 
  dcast(., type~source) %>% 
  mutate(per_motif_ins_fold_change = motif/control) %>% 
  select(type, per_motif_ins_fold_change) %>% 
  melt() %>% 
  mutate(tech = "Illumina") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_Illumina_moderate %>% 
  select(type, source, per_motif_del_rate) %>% 
  dcast(., type~source) %>% 
  mutate(per_motif_del_fold_change = motif/control) %>% 
  select(type, per_motif_del_fold_change) %>% 
  melt() %>% 
  mutate(tech = "Illumina") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_HiFi_moderate %>% 
  select(type, source, per_motif_smm_rate) %>% 
  dcast(., type~source) %>% 
  mutate(per_motif_smm_fold_change = motif/control) %>% 
  select(type, per_motif_smm_fold_change) %>% 
  melt() %>%
  mutate(tech = "HiFi") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_HiFi_moderate %>% 
  select(type, source, per_motif_ins_rate) %>% 
  dcast(., type~source) %>% 
  mutate(per_motif_ins_fold_change = motif/control) %>% 
  select(type, per_motif_ins_fold_change) %>% 
  melt() %>% 
  mutate(tech = "HiFi") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_HiFi_moderate %>% 
  select(type, source, per_motif_del_rate) %>% 
  dcast(., type~source) %>% 
  mutate(per_motif_del_fold_change = motif/control) %>% 
  select(type, per_motif_del_fold_change) %>% 
  melt() %>%
  mutate(tech = "HiFi") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_ONT_moderate %>% 
  select(type, source, per_motif_smm_rate) %>% 
  dcast(., type~source) %>% 
  mutate(per_motif_smm_fold_change = motif/control) %>% 
  select(type, per_motif_smm_fold_change) %>% 
  melt() %>%
  mutate(tech = "ONT") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_ONT_moderate %>% 
  select(type, source, per_motif_ins_rate) %>% 
  dcast(., type~source) %>% 
  mutate(per_motif_ins_fold_change = motif/control) %>% 
  select(type, per_motif_ins_fold_change) %>% 
  melt() %>% 
  mutate(tech = "ONT") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_ONT_moderate %>% 
  select(type, source, per_motif_del_rate) %>% 
  dcast(., type~source) %>% 
  mutate(per_motif_del_fold_change = motif/control) %>% 
  select(type, per_motif_del_fold_change) %>% 
  melt() %>%  
  mutate(tech = "ONT") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

fold_change_table$tech <- factor(fold_change_table$tech, levels = c("Illumina", "HiFi", "ONT"))
fold_change_table %>% 
  mutate(filter_level = "moderate") -> fold_change_table_moderate
#####

# Stringent filter
#####
raw_error_rates_Illumina_ultrastringent %>% 
  select(type, source, raw_smm_rate) %>% 
  dcast(., type~source) %>% 
  mutate(aggregate_smm_fold_change = motif/control) %>% 
  select(type, aggregate_smm_fold_change) %>% 
  melt() %>%
  mutate(tech = "Illumina") -> fold_change_table

raw_error_rates_Illumina_ultrastringent %>% 
  select(type, source, raw_ins_rate) %>% 
  dcast(., type~source) %>% 
  mutate(aggregate_ins_fold_change = motif/control) %>% 
  select(type, aggregate_ins_fold_change) %>% 
  melt() %>% 
  mutate(tech = "Illumina") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_Illumina_ultrastringent %>% 
  select(type, source, raw_del_rate) %>% 
  dcast(., type~source) %>% 
  mutate(aggregate_del_fold_change = motif/control) %>% 
  select(type, aggregate_del_fold_change) %>% 
  melt() %>% 
  mutate(tech = "Illumina") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_HiFi_ultrastringent %>% 
  select(type, source, raw_smm_rate) %>% 
  dcast(., type~source) %>% 
  mutate(aggregate_smm_fold_change = motif/control) %>% 
  select(type, aggregate_smm_fold_change) %>% 
  melt() %>%
  mutate(tech = "HiFi") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_HiFi_ultrastringent %>% 
  select(type, source, raw_ins_rate) %>% 
  dcast(., type~source) %>% 
  mutate(aggregate_ins_fold_change = motif/control) %>% 
  select(type, aggregate_ins_fold_change) %>% 
  melt() %>% 
  mutate(tech = "HiFi") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_HiFi_ultrastringent %>% 
  select(type, source, raw_del_rate) %>% 
  dcast(., type~source) %>% 
  mutate(aggregate_del_fold_change = motif/control) %>% 
  select(type, aggregate_del_fold_change) %>% 
  melt() %>%
  mutate(tech = "HiFi") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_ONT_ultrastringent %>% 
  select(type, source, raw_smm_rate) %>% 
  dcast(., type~source) %>% 
  mutate(aggregate_smm_fold_change = motif/control) %>% 
  select(type, aggregate_smm_fold_change) %>% 
  melt() %>%
  mutate(tech = "ONT") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_ONT_ultrastringent %>% 
  select(type, source, raw_ins_rate) %>% 
  dcast(., type~source) %>% 
  mutate(aggregate_ins_fold_change = motif/control) %>% 
  select(type, aggregate_ins_fold_change) %>% 
  melt() %>% 
  mutate(tech = "ONT") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_ONT_ultrastringent %>% 
  select(type, source, raw_del_rate) %>% 
  dcast(., type~source) %>% 
  mutate(aggregate_del_fold_change = motif/control) %>% 
  select(type, aggregate_del_fold_change) %>% 
  melt() %>%  
  mutate(tech = "ONT") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_Illumina_ultrastringent %>% 
  select(type, source, per_motif_smm_rate) %>% 
  dcast(., type~source) %>% 
  mutate(per_motif_smm_fold_change = motif/control) %>% 
  select(type, per_motif_smm_fold_change) %>% 
  melt() %>%
  mutate(tech = "Illumina") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_Illumina_ultrastringent %>% 
  select(type, source, per_motif_ins_rate) %>% 
  dcast(., type~source) %>% 
  mutate(per_motif_ins_fold_change = motif/control) %>% 
  select(type, per_motif_ins_fold_change) %>% 
  melt() %>% 
  mutate(tech = "Illumina") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_Illumina_ultrastringent %>% 
  select(type, source, per_motif_del_rate) %>% 
  dcast(., type~source) %>% 
  mutate(per_motif_del_fold_change = motif/control) %>% 
  select(type, per_motif_del_fold_change) %>% 
  melt() %>% 
  mutate(tech = "Illumina") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_HiFi_ultrastringent %>% 
  select(type, source, per_motif_smm_rate) %>% 
  dcast(., type~source) %>% 
  mutate(per_motif_smm_fold_change = motif/control) %>% 
  select(type, per_motif_smm_fold_change) %>% 
  melt() %>%
  mutate(tech = "HiFi") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_HiFi_ultrastringent %>% 
  select(type, source, per_motif_ins_rate) %>% 
  dcast(., type~source) %>% 
  mutate(per_motif_ins_fold_change = motif/control) %>% 
  select(type, per_motif_ins_fold_change) %>% 
  melt() %>% 
  mutate(tech = "HiFi") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_HiFi_ultrastringent %>% 
  select(type, source, per_motif_del_rate) %>% 
  dcast(., type~source) %>% 
  mutate(per_motif_del_fold_change = motif/control) %>% 
  select(type, per_motif_del_fold_change) %>% 
  melt() %>%
  mutate(tech = "HiFi") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_ONT_ultrastringent %>% 
  select(type, source, per_motif_smm_rate) %>% 
  dcast(., type~source) %>% 
  mutate(per_motif_smm_fold_change = motif/control) %>% 
  select(type, per_motif_smm_fold_change) %>% 
  melt() %>%
  mutate(tech = "ONT") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_ONT_ultrastringent %>% 
  select(type, source, per_motif_ins_rate) %>% 
  dcast(., type~source) %>% 
  mutate(per_motif_ins_fold_change = motif/control) %>% 
  select(type, per_motif_ins_fold_change) %>% 
  melt() %>% 
  mutate(tech = "ONT") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

raw_error_rates_ONT_ultrastringent %>% 
  select(type, source, per_motif_del_rate) %>% 
  dcast(., type~source) %>% 
  mutate(per_motif_del_fold_change = motif/control) %>% 
  select(type, per_motif_del_fold_change) %>% 
  melt() %>%  
  mutate(tech = "ONT") %>%
  rbind.data.frame(., fold_change_table) -> fold_change_table

fold_change_table$tech <- factor(fold_change_table$tech, levels = c("Illumina", "HiFi", "ONT"))
fold_change_table %>% 
  mutate(filter_level = "stringent") -> fold_change_table_ultrastringent
#####

# Subsection error rates
# Moderate filter 
#####
Illumina_motif_subsections %>%
  group_by(type, simple_section) %>%
  summarize(aggregate_mean_smm = sum(smm)/sum(depth),
            aggregate_mean_ins = sum(ins)/sum(depth),
            aggregate_mean_del = sum(del)/sum(depth),
            per_motif_mean_smm = mean(smm_rate, na.rm = T),
            per_motif_mean_ins = mean(ins_rate, na.rm = T),
            per_motif_mean_del = mean(del_rate, na.rm = T)) %>% 
  data.frame() -> Illumina_motif_subsections_raw_error_rates
Illumina_motif_subsections_raw_error_rates %>%
  mutate(tech = "Illumina") -> Illumina_motif_subsections_raw_error_rates

HiFi_motif_subsections %>%
  group_by(type, simple_section) %>%
  summarize(aggregate_mean_smm = sum(smm)/sum(depth),
            aggregate_mean_ins = sum(ins)/sum(depth),
            aggregate_mean_del = sum(del)/sum(depth),
            per_motif_mean_smm = mean(smm_rate, na.rm = T),
            per_motif_mean_ins = mean(ins_rate, na.rm = T),
            per_motif_mean_del = mean(del_rate, na.rm = T)) %>% 
  data.frame() -> HiFi_motif_subsections_raw_error_rates
HiFi_motif_subsections_raw_error_rates %>%
  mutate(tech = "HiFi") -> HiFi_motif_subsections_raw_error_rates


ONT_motif_subsections %>%
  group_by(type, simple_section) %>%
  summarize(aggregate_mean_smm = sum(smm)/sum(depth),
            aggregate_mean_ins = sum(ins)/sum(depth),
            aggregate_mean_del = sum(del)/sum(depth),
            per_motif_mean_smm = mean(smm_rate, na.rm = T),
            per_motif_mean_ins = mean(ins_rate, na.rm = T),
            per_motif_mean_del = mean(del_rate, na.rm = T)) %>% 
  data.frame() -> ONT_motif_subsections_raw_error_rates
ONT_motif_subsections_raw_error_rates %>%
  mutate(tech = "ONT") -> ONT_motif_subsections_raw_error_rates

rbind.data.frame(
  Illumina_motif_subsections_raw_error_rates,
  HiFi_motif_subsections_raw_error_rates,
  ONT_motif_subsections_raw_error_rates
) %>%
  melt(., measure.vars = c("aggregate_mean_smm", "aggregate_mean_ins",  "aggregate_mean_del", "per_motif_mean_smm", "per_motif_mean_ins", "per_motif_mean_del")) %>%
  mutate(filter_level = "moderate",
         error_type = ifelse(grepl("smm", variable), "smm", 
                             ifelse(grepl("ins", variable), "ins",
                                    "del")),
         mean_type = ifelse(grepl("agg", variable), "aggregate_mean", "per_motif_mean")) %>%
  select(type,error_type,tech, simple_section, mean_type, value) %>% 
  dcast(.,  tech+type+error_type+simple_section~mean_type, value.var = "value") %>%
  mutate(filter_level = "moderate") %>%
  select(type, simple_section, error_type, tech, filter_level, aggregate_mean, per_motif_mean) -> sequencing_success_table_motif_simple_subsections_combined

#sequencing_success_table_motif_simple_subsections_combined %>%
#  write.table(., file = "sequencing_success_table_motif_simple_subsections_combined.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#
#sequencing_success_table_motif_simple_subsections_combined <- read.table("sequencing_success_table_motif_simple_subsections_combined.txt", header = T)

#####

# Stringent filter
#####
Illumina_motif_subsections_ultrastringent %>%
  group_by(type, simple_section) %>%
  summarize(aggregate_mean_smm = sum(smm)/sum(depth),
            aggregate_mean_ins = sum(ins)/sum(depth),
            aggregate_mean_del = sum(del)/sum(depth),
            per_motif_mean_smm = mean(smm_rate, na.rm = T),
            per_motif_mean_ins = mean(ins_rate, na.rm = T),
            per_motif_mean_del = mean(del_rate, na.rm = T)) %>% 
  data.frame() -> Illumina_motif_subsections_ultrastringent_raw_error_rates
Illumina_motif_subsections_ultrastringent_raw_error_rates %>%
  mutate(tech = "Illumina") -> Illumina_motif_subsections_ultrastringent_raw_error_rates

HiFi_motif_subsections_ultrastringent %>%
  group_by(type, simple_section) %>%
  summarize(aggregate_mean_smm = sum(smm)/sum(depth),
            aggregate_mean_ins = sum(ins)/sum(depth),
            aggregate_mean_del = sum(del)/sum(depth),
            per_motif_mean_smm = mean(smm_rate, na.rm = T),
            per_motif_mean_ins = mean(ins_rate, na.rm = T),
            per_motif_mean_del = mean(del_rate, na.rm = T)) %>% 
  data.frame() -> HiFi_motif_subsections_ultrastringent_raw_error_rates
HiFi_motif_subsections_ultrastringent_raw_error_rates %>%
  mutate(tech = "HiFi") -> HiFi_motif_subsections_ultrastringent_raw_error_rates

ONT_motif_subsections_ultrastringent %>%
  group_by(type, simple_section) %>%
  summarize(aggregate_mean_smm = sum(smm)/sum(depth),
            aggregate_mean_ins = sum(ins)/sum(depth),
            aggregate_mean_del = sum(del)/sum(depth),
            per_motif_mean_smm = mean(smm_rate, na.rm = T),
            per_motif_mean_ins = mean(ins_rate, na.rm = T),
            per_motif_mean_del = mean(del_rate, na.rm = T)) %>% 
  data.frame() -> ONT_motif_subsections_ultrastringent_raw_error_rates
ONT_motif_subsections_ultrastringent_raw_error_rates %>%
  mutate(tech = "ONT") -> ONT_motif_subsections_ultrastringent_raw_error_rates

rbind.data.frame(
  Illumina_motif_subsections_ultrastringent_raw_error_rates,
  HiFi_motif_subsections_ultrastringent_raw_error_rates,
  ONT_motif_subsections_ultrastringent_raw_error_rates
) %>%
  melt(., measure.vars = c("aggregate_mean_smm", "aggregate_mean_ins",  "aggregate_mean_del", "per_motif_mean_smm", "per_motif_mean_ins", "per_motif_mean_del")) %>%
  mutate(
         error_type = ifelse(grepl("smm", variable), "smm", 
                             ifelse(grepl("ins", variable), "ins",
                                    "del")),
         mean_type = ifelse(grepl("agg", variable), "aggregate_mean", "per_motif_mean")) %>%
  select(type,error_type,tech, simple_section, mean_type, value) %>% 
  dcast(.,  tech+type+error_type+simple_section~mean_type, value.var = "value") %>%
  mutate(filter_level = "stringent") %>%
  select(type, simple_section, error_type, tech, filter_level, aggregate_mean, per_motif_mean) -> sequencing_success_table_motif_simple_subsections_ultrastringent_combined

#sequencing_success_table_motif_simple_subsections_ultrastringent_combined %>%
#  write.table(., file = "sequencing_success_table_motif_simple_subsections_ultrastringent_combined.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#
#sequencing_success_table_motif_simple_subsections_ultrastringent_combined <- read.table("sequencing_success_table_motif_simple_subsections_ultrastringent_combined.txt", header = T)

sequencing_success_table_motif_simple_subsections_ultrastringent_combined %>%
  dcast(., type+error_type+tech~simple_section, value.var = c("aggregate_mean")) %>%
  mutate(fold_change_aggregate_mean = repeat_stem/spacer_loop) 

# Fold-change table

# Fold change table subsections

# Moderate filter
sequencing_success_table_motif_simple_subsections_combined %>%
  dcast(., type+error_type+tech~simple_section, value.var = c("aggregate_mean")) %>%
  mutate(fold_change_rssl_aggregate_mean = repeat_stem/spacer_loop,
         fold_change_slrs_aggregate_mean = spacer_loop/repeat_stem) -> temp

sequencing_success_table_motif_simple_subsections_combined %>%
  dcast(., type+error_type+tech~simple_section, value.var = c("per_motif_mean")) %>%
  mutate(fold_change_rssl_per_motif_mean = repeat_stem/spacer_loop,
         fold_change_slrs_per_motif_mean = spacer_loop/repeat_stem) %>%
  select(fold_change_rssl_per_motif_mean, fold_change_slrs_per_motif_mean) %>%
  cbind.data.frame(temp, .) %>%
  select(type, error_type, tech, 
         fold_change_rssl_aggregate_mean,
         fold_change_slrs_aggregate_mean,
         fold_change_rssl_per_motif_mean,
         fold_change_slrs_per_motif_mean) %>%
  mutate(filter_level = "moderate") -> moderate_filter_fold_change_simple_subsections_combined

moderate_filter_fold_change_simple_subsections_combined %>%
  write.table(., file = "moderate_filter_fold_change_simple_subsections_combined.txt", sep = "\t", quote = FALSE, row.names = FALSE)
moderate_filter_fold_change_simple_subsections_combined <- read.table("moderate_filter_fold_change_simple_subsections_combined.txt", header = T)

# Stringent filter
sequencing_success_table_motif_simple_subsections_ultrastringent_combined %>%
  dcast(., type+error_type+tech~simple_section, value.var = c("aggregate_mean")) %>%
  mutate(fold_change_rssl_aggregate_mean = repeat_stem/spacer_loop,
         fold_change_slrs_aggregate_mean = spacer_loop/repeat_stem) -> temp

sequencing_success_table_motif_simple_subsections_ultrastringent_combined %>%
  dcast(., type+error_type+tech~simple_section, value.var = c("per_motif_mean")) %>%
  mutate(fold_change_rssl_per_motif_mean = repeat_stem/spacer_loop,
         fold_change_slrs_per_motif_mean = spacer_loop/repeat_stem) %>%
  select(fold_change_rssl_per_motif_mean, fold_change_slrs_per_motif_mean) %>%
  cbind.data.frame(temp, .) %>%
  select(type, error_type, tech, 
         fold_change_rssl_aggregate_mean,
         fold_change_slrs_aggregate_mean,
         fold_change_rssl_per_motif_mean,
         fold_change_slrs_per_motif_mean) %>%
  mutate(filter_level = "stringent") -> fold_change_simple_subsections_ultrastringent_combined

fold_change_simple_subsections_ultrastringent_combined %>%  
  write.table(., file = "fold_change_simple_subsections_ultrastringent_combined.txt", sep = "\t", quote = FALSE, row.names = FALSE)
fold_change_simple_subsections_ultrastringent_combined <-  read.table("fold_change_simple_subsections_ultrastringent_combined.txt", header = T)
#####

#################################################
# Read depth per bp
#################################################

# Moderate filter
#####
Illumina_curated %>%
  group_by(type, source) %>%
  summarize(raw_depth_per_bp = sum(depth)/sum(length)) %>%
  data.frame() %>%
  mutate(tech = "Illumina") -> raw_depth_per_bp_Illumina

Illumina_curated %>%
  group_by(type, source) %>%
  summarize(raw_depth_per_bp = sum(depth)/sum(length),
            per_motif_depth_per_bp = mean(depth_per_bp, na.rm = T)) %>%
  data.frame() %>%
  dcast(., type~source) %>%
  mutate(fold_change_depth_per_bp = motif/control) %>%
  select(type,fold_change_depth_per_bp) %>%
  mutate(tech = "Illumina") -> depth_per_bp_Illumina

HiFi_curated %>%
  group_by(type, source) %>%
  summarize(raw_depth_per_bp = sum(depth)/sum(length)) %>%
  data.frame() %>%
  mutate(tech = "HiFi") -> raw_depth_per_bp_HiFi

HiFi_curated %>%
  group_by(type, source) %>%
  summarize(raw_depth_per_bp = sum(depth)/sum(length),
            per_motif_depth_per_bp = mean(depth_per_bp, na.rm = T)) %>%
  data.frame() %>%
  dcast(., type~source) %>%
  mutate(fold_change_depth_per_bp = motif/control) %>%
  select(type, fold_change_depth_per_bp) %>%
  mutate(tech = "HiFi") -> depth_per_bp_HiFi

ONT_curated %>%
  group_by(type, source) %>%
  summarize(raw_depth_per_bp = sum(depth)/sum(length)) %>%
  data.frame() %>%
  mutate(tech = "ONT") -> raw_depth_per_bp_ONT

ONT_curated %>%
  group_by(type, source) %>%
  summarize(raw_depth_per_bp = sum(depth)/sum(length),
            per_motif_depth_per_bp = mean(depth_per_bp, na.rm = T)) %>%
  data.frame() %>%
  dcast(., type~source) %>%
  mutate(fold_change_depth_per_bp = motif/control) %>%
  select(type, fold_change_depth_per_bp) %>%
  mutate(tech = "ONT") -> depth_per_bp_ONT

rbind.data.frame(depth_per_bp_Illumina,
                 depth_per_bp_HiFi,
                 depth_per_bp_ONT) %>% 
  mutate(filter_level = "moderate") -> fold_change_depth_moderate
#####

# Stringent filter
#####
Illumina_ultrastringent_filter %>%
  group_by(type, source) %>%
  summarize(raw_depth_per_bp = sum(depth)/sum(length)) %>%
  data.frame() %>%
  mutate(tech = "Illumina") -> raw_depth_per_bp_Illumina_stringent  

Illumina_ultrastringent_filter %>%
  group_by(type, source) %>%
  summarize(raw_depth_per_bp = sum(depth)/sum(length)) %>%
  data.frame() %>%
  dcast(., type~source) %>%
  mutate(fold_change_depth_per_bp = motif/control) %>%
  select(type,fold_change_depth_per_bp) %>%
  mutate(tech = "Illumina") -> depth_per_bp_Illumina

HiFi_ultrastringent_filter %>%
  group_by(type, source) %>%
  summarize(raw_depth_per_bp = sum(depth)/sum(length)) %>%
  data.frame() %>%
  mutate(tech = "HiFi") -> raw_depth_per_bp_HiFi_stringent

HiFi_ultrastringent_filter %>%
  group_by(type, source) %>%
  summarize(raw_depth_per_bp = sum(depth)/sum(length)) %>%
  data.frame() %>%
  dcast(., type~source) %>%
  mutate(fold_change_depth_per_bp = motif/control) %>%
  select(type, fold_change_depth_per_bp) %>%
  mutate(tech = "HiFi") -> depth_per_bp_HiFi

ONT_ultrastringent_filter %>%
  group_by(type, source) %>%
  summarize(raw_depth_per_bp = sum(depth)/sum(length)) %>%
  data.frame() %>%
  mutate(tech = "ONT") -> raw_depth_per_bp_ONT_stringent

ONT_ultrastringent_filter %>%
  group_by(type, source) %>%
  summarize(raw_depth_per_bp = sum(depth)/sum(length)) %>%
  data.frame() %>%
  dcast(., type~source) %>%
  mutate(fold_change_depth_per_bp = motif/control) %>%
  select(type, fold_change_depth_per_bp) %>%
  mutate(tech = "ONT") -> depth_per_bp_ONT

rbind.data.frame(depth_per_bp_Illumina,
                 depth_per_bp_HiFi,
                 depth_per_bp_ONT) %>%
  mutate(filter_level = "stringent") -> fold_change_depth_ultrastringent
#####

rbind.data.frame(
  fold_change_depth_moderate,
  fold_change_depth_ultrastringent
) -> fold_change_depth_table_combined

fold_change_depth_table_combined$filter_level <-  as.factor(fold_change_depth_table_combined$filter_level)

#################################################
# Base quality
#################################################

# Aggregate mean base quality

#####

Illumina_curated %>%
  data.frame() %>%
  filter(!is.na(basequality)) %>%
  summarize(mean_bq = sum(basequality / n()))

Illumina_curated %>%
  data.frame() %>%
  filter(!is.na(basequality)) %>%
  filter(source == "control") %>% 
  summarize(mean_bq = sum(basequality / n()))

Illumina_curated %>%
  data.frame() %>%
  filter(!is.na(basequality)) %>%
  group_by(type, source) %>%
  summarize(basequality = sum(basequality)/n()) %>%
  data.frame() %>%
  dcast(., type~source) %>%
  mutate(fold_change_basequality = motif/control) %>%
  select(type,fold_change_basequality) %>%
  mutate(tech = "Illumina") -> basequality_Illumina

HiFi_curated %>%
  data.frame() %>%
  filter(!is.na(basequality)) %>%
  summarize(mean_bq = sum(basequality / n()))

HiFi_curated %>%
  data.frame() %>%
  filter(!is.na(basequality)) %>%
  filter(source == "control") %>% 
  summarize(mean_bq = sum(basequality / n()))

HiFi_curated %>%
  data.frame() %>%
  filter(!is.na(basequality)) %>%
  group_by(type, source) %>%
  summarize(basequality = sum(basequality)/n()) %>%
  data.frame() %>%
  dcast(., type~source) %>%
  mutate(fold_change_basequality = motif/control) %>%
  select(type,fold_change_basequality) %>%
  mutate(tech = "HiFi") -> basequality_HiFi

rbind.data.frame(basequality_Illumina,
                 basequality_HiFi
) -> fold_change_basequality_table_combined
#####

#################################################
# Significance tests
#################################################

# SNM, INS, DEL, depth, base quality
#####

# Moderate filter level

Illumina_curated %>% 
  group_by(type) %>%
  summarize_each(funs(t.test(.[source == "motif"], .[source=="control"])$p.value), vars = c("smm_rate", "ins_rate", "del_rate", "depth_per_bp")) %>%
  data.frame() %>%
  rename(., smm = vars1, ins = vars2, del = vars3, depth_per_bp = vars4) %>%
  mutate(tech = "Illumina", filter_level = "moderate") -> Illumina_pvalues_moderate

HiFi_curated %>% 
  group_by(type) %>%
  summarize_each(funs(t.test(.[source == "motif"], .[source=="control"])$p.value), vars = c("smm_rate", "ins_rate", "del_rate", "depth_per_bp")) %>%
  data.frame() %>%
  rename(., smm = vars1, ins = vars2, del = vars3, depth_per_bp = vars4) %>%
  mutate(tech = "HiFi", filter_level = "moderate") -> HiFi_pvalues_moderate

ONT_curated %>% 
  group_by(type) %>%
  summarize_each(funs(t.test(.[source == "motif"], .[source=="control"])$p.value), vars = c("smm_rate", "ins_rate", "del_rate", "depth_per_bp")) %>%
  data.frame() %>%
  rename(., smm = vars1, ins = vars2, del = vars3, depth_per_bp = vars4) %>%
  mutate(tech = "ONT", filter_level = "moderate") -> ONT_pvalues_moderate

# Stringent filter level

Illumina_ultrastringent_filter %>% 
  group_by(type) %>%
  summarize_each(funs(t.test(.[source == "motif"], .[source=="control"])$p.value), vars = c("smm_rate", "ins_rate", "del_rate", "depth_per_bp")) %>%
  data.frame() %>%
  rename(., smm = vars1, ins = vars2, del = vars3, depth_per_bp = vars4) %>%
  mutate(tech = "Illumina", filter_level = "stringent") -> Illumina_pvalues_ultrastringent

HiFi_ultrastringent_filter %>% 
  group_by(type) %>%
  summarize_each(funs(t.test(.[source == "motif"], .[source=="control"])$p.value), vars = c("smm_rate", "ins_rate", "del_rate", "depth_per_bp")) %>%
  data.frame() %>%
  rename(., smm = vars1, ins = vars2, del = vars3, depth_per_bp = vars4) %>%
  mutate(tech = "HiFi", filter_level = "stringent") -> HiFi_pvalues_ultrastringent

ONT_ultrastringent_filter %>% 
  group_by(type) %>%
  summarize_each(funs(t.test(.[source == "motif"], .[source=="control"])$p.value), vars = c("smm_rate", "ins_rate", "del_rate", "depth_per_bp")) %>%
  data.frame() %>%
  rename(., smm = vars1, ins = vars2, del = vars3, depth_per_bp = vars4) %>%
  mutate(tech = "ONT", filter_level = "stringent") -> ONT_pvalues_ultrastringent

rbind.data.frame(
  Illumina_pvalues_moderate,
  HiFi_pvalues_moderate,
  ONT_pvalues_moderate,
  Illumina_pvalues_ultrastringent,
  HiFi_pvalues_ultrastringent,
  ONT_pvalues_ultrastringent) ->p_value_table

p_value_table$smm <- as.character(p_value_table$smm)
p_value_table$ins <- as.character(p_value_table$ins)
p_value_table$del <- as.character(p_value_table$del)

p_value_table %>% 
  mutate(
    smm_p_adjusted = p.adjust(smm, method = "hochberg"),
    ins_p_adjusted = p.adjust(ins, method = "hochberg"),
    del_p_adjusted = p.adjust(del, method = "hochberg"),
    depth_p_adjusted = p.adjust(depth_per_bp, method = "hochberg")
  ) -> p_value_table

names<- names(p_value_table)

p_value_table %>%
  select(type, tech, filter_level, smm_p_adjusted, ins_p_adjusted, del_p_adjusted, depth_p_adjusted) %>%
  melt() -> molten_p_value_table

molten_p_value_table %>%
  mutate(error_type = ifelse(grepl("smm", variable), "smm", 
                             ifelse(grepl("del", variable), "del", "ins"))) -> molten_p_value_table

# Base quality p-values

Illumina_curated %>% 
  group_by(type) %>%
  summarize_each(funs(t.test(.[source == "motif"], .[source=="control"])$p.value), vars = c("basequality")) %>%
  data.frame() %>%
  rename(., basequality = vars) %>%
  mutate(tech = "Illumina", filter_level = "moderate") -> Illumina_pvalues_basequality

HiFi_curated %>% 
  group_by(type) %>%
  summarize_each(funs(t.test(.[source == "motif"], .[source=="control"])$p.value), vars = c("basequality")) %>%
  data.frame() %>%
  rename(., basequality = vars) %>%
  mutate(tech = "HiFi", filter_level = "moderate") -> HiFi_pvalues_basequality


rbind.data.frame( Illumina_pvalues_basequality,
                  HiFi_pvalues_basequality) %>% 
  select(-filter_level) -> pvalues_basequality_combined

pvalues_basequality_combined %>% 
  mutate( basequality_p_adjusted = p.adjust(basequality, method = "hochberg")) %>% 
  select(-basequality) -> pvalues_adjusted_basequality

fold_change_basequality_table_combined %>% 
  merge(., pvalues_adjusted_basequality, by = c("type", "tech" )) %>% 
  mutate(fontface = ifelse(basequality_p_adjusted <= 0.05, 2, 1))   -> fold_change_basequality_table_combined_with_fontface 
#####

# Motif subsections significance tests
#####

Illumina_motif_subsections %>% 
  group_by(type) %>%
  summarize_each(funs(t.test(.[simple_section == "spacer_loop"], .[simple_section=="repeat_stem"])$p.value), vars = smm_rate:ins_rate:del_rate) %>%
  data.frame() %>%
  rename(., smm = vars1, ins = vars2, del = vars3) %>%
  mutate(tech = "Illumina", filter_level = "moderate") -> Illumina_subsections_pvalues_moderate

HiFi_motif_subsections %>% 
  group_by(type) %>%
  summarize_each(funs(t.test(.[simple_section == "spacer_loop"], .[simple_section=="repeat_stem"])$p.value), vars = smm_rate:ins_rate:del_rate) %>%
  data.frame() %>%
  rename(., smm = vars1, ins = vars2, del = vars3) %>%
  mutate(tech = "HiFi", filter_level = "moderate") -> HiFi_subsections_pvalues_moderate

ONT_motif_subsections %>% 
  group_by(type) %>%
  summarize_each(funs(t.test(.[simple_section == "spacer_loop"], .[simple_section=="repeat_stem"])$p.value), vars = smm_rate:ins_rate:del_rate) %>%
  data.frame() %>%
  rename(., smm = vars1, ins = vars2, del = vars3) %>%
  mutate(tech = "ONT", filter_level = "moderate") -> ONT_subsections_pvalues_moderate

Illumina_motif_subsections_ultrastringent %>% 
  group_by(type) %>%
  summarize_each(funs(t.test(.[simple_section == "spacer_loop"], .[simple_section=="repeat_stem"])$p.value), vars = smm_rate:ins_rate:del_rate) %>%
  data.frame() %>%
  rename(., smm = vars1, ins = vars2, del = vars3) %>%
  mutate(tech = "Illumina", filter_level = "stringent") -> Illumina_subsections_pvalues_ultrastringent

HiFi_motif_subsections_ultrastringent %>% 
  group_by(type) %>%
  summarize_each(funs(t.test(.[simple_section == "spacer_loop"], .[simple_section=="repeat_stem"])$p.value), vars = smm_rate:ins_rate:del_rate) %>%
  data.frame() %>%
  rename(., smm = vars1, ins = vars2, del = vars3) %>%
  mutate(tech = "HiFi", filter_level = "stringent") -> HiFi_subsections_pvalues_ultrastringent

ONT_motif_subsections_ultrastringent %>% 
  group_by(type) %>%
  summarize_each(funs(t.test(.[simple_section == "spacer_loop"], .[simple_section=="repeat_stem"])$p.value), vars = smm_rate:ins_rate:del_rate) %>%
  data.frame() %>%
  rename(., smm = vars1, ins = vars2, del = vars3) %>%
  mutate(tech = "ONT", filter_level = "stringent") -> ONT_subsections_pvalues_ultrastringent

rbind.data.frame(
  Illumina_subsections_pvalues_moderate,
  HiFi_subsections_pvalues_moderate,
  ONT_subsections_pvalues_moderate,
  Illumina_subsections_pvalues_ultrastringent,
  HiFi_subsections_pvalues_ultrastringent,
  ONT_subsections_pvalues_ultrastringent) %>%
  mutate(
    smm_p_adjusted = p.adjust(smm, method = "hochberg"),
    ins_p_adjusted = p.adjust(ins, method = "hochberg"),
    del_p_adjusted = p.adjust(del, method = "hochberg")
  ) -> p_value_table_motif_subsections_filter_levels
#####

#################################################
# Error rates Figures
#################################################

#################################################
# Single-nucleotide mismatches

## Illumina
#####
Illumina_curated %>% 
  ggplot(aes(y = smm_rate, x = type, fill = source)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  stat_summary(fun.y=mean, geom="point", shape=20, position = position_dodge(width = 0.9), aes(group = source, color = "Per motif rate"), size = 5) + 
  geom_point(data = raw_error_rates_Illumina_moderate, aes(x = type, y = raw_smm_rate, group = source, color = "Aggregate rate"), 
             shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = 9) +
  ylab("SNM per nucleotide") + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(labels = number_format(accuracy = 0.001,decimal.mark = '.')) +
  ggtitle("Illumina") +
  scale_colour_manual(values = c("Aggregate rate" = "#feb24c"  , "Per motif rate" = "#1a1a1a"  )) +
  guides(colour = guide_legend(override.aes = list(shape = c(17, 20))),
         fill = guide_legend(override.aes = list(linetype = c(0, 0),
                                                 shape = c(NA, NA)))) + 
  theme( legend.position = "bottom",
         plot.title = element_text(hjust = 0.5),
         legend.title = element_blank(),
         axis.title.x = element_blank(),
         axis.text.x = element_blank()) -> Illumina_smm_boxplots
Illumina_smm_ymax_moderate <- layer_scales(Illumina_smm_boxplots)$y$range$range[2]

Illumina_ultrastringent_filter %>% 
  ggplot(aes(y = smm_rate, x = type, fill = source)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  stat_summary(fun.y=mean, geom="point", shape=20, position = position_dodge(width = 0.9), aes(group = source, color = "Per motif rate"), size = 5) + 
  geom_point(data = raw_error_rates_Illumina_ultrastringent, aes(x = type, y = raw_smm_rate, group = source, color = "Aggregate rate"), 
             shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = 9) +
  ylab("SNM per nucleotide") + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(labels = number_format(accuracy = 0.001,decimal.mark = '.')) + #, limits = c(0, Illumina_smm_ymax_moderate)) +
  coord_cartesian(ylim = c(0, Illumina_smm_ymax_moderate)) + 
  scale_colour_manual(values = c("Aggregate rate" = "#feb24c"  , "Per motif rate" = "#1a1a1a"  )) +
  guides(colour = guide_legend(override.aes = list(shape = c(17, 20))),
         fill = guide_legend(override.aes = list(linetype = c(0, 0),
                                                 shape = c(NA, NA)))) + 
  ggtitle("Illumina") +
  theme( legend.position = "bottom",
         plot.title = element_text(hjust = 0.5),
         legend.title = element_blank(),
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.text.x = element_blank()) -> Illumina_ultrastringent_smm_boxplots
#####

## HiFi
#####
HiFi_curated %>% 
  ggplot(aes(y = smm_rate, x = type, fill = source)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  stat_summary(fun.y=mean, geom="point", shape=20, position = position_dodge(width = 0.9), aes(group = source, color = "Per motif rate"), size = 5) + 
  geom_point(data = raw_error_rates_HiFi_moderate, aes(x = type, y = raw_smm_rate, group = source, color = "Aggregate rate"), 
             shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = 9) +
  ylab("SNM per nucleotide") + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(labels = number_format(accuracy = 0.001,decimal.mark = '.')) +
  ggtitle("HiFi") +
  scale_colour_manual(values = c("Aggregate rate" = "#feb24c"  , "Per motif rate" = "#1a1a1a"  )) +
  guides(colour = guide_legend(override.aes = list(shape = c(17, 20))),
         fill = guide_legend(override.aes = list(linetype = c(0, 0),
                                                 shape = c(NA, NA)))) + 
  theme( legend.position = "bottom",
         plot.title = element_text(hjust = 0.5),
         legend.title = element_blank(),
         axis.title.x = element_blank(),
         axis.text.x = element_blank()) -> HiFi_smm_boxplots
HiFi_smm_ymax_moderate <- layer_scales(HiFi_smm_boxplots)$y$range$range[2]

HiFi_ultrastringent_filter %>% 
  ggplot(aes(y = smm_rate, x = type, fill = source)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  stat_summary(fun.y=mean, geom="point", shape=20, position = position_dodge(width = 0.9), aes(group = source, color = "Per motif rate"), size = 5) + 
  geom_point(data = raw_error_rates_HiFi_ultrastringent, aes(x = type, y = raw_smm_rate, group = source, color = "Aggregate rate"), 
             shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = 9) +
  ylab("smm per nucleotide") + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001,decimal.mark = '.') ) + #, limits = c(0, HiFi_smm_ymax_moderate)) +
  scale_colour_manual(values = c("Aggregate rate" = "#feb24c"  , "Per motif rate" = "#1a1a1a"  )) +
  coord_cartesian(ylim = c(0, HiFi_smm_ymax_moderate)) + 
  guides(colour = guide_legend(override.aes = list(shape = c(17, 20))),
         fill = guide_legend(override.aes = list(linetype = c(0, 0),
                                                 shape = c(NA, NA)))) + 
  ggtitle("HiFi") +
  theme( legend.position = "bottom",
         plot.title = element_text(hjust = 0.5),
         legend.title = element_blank(),
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.text.x = element_blank()) -> HiFi_ultrastringent_smm_boxplots
#####

## ONT
#####
ONT_curated %>% 
  ggplot(aes(y = smm_rate, x = type, fill = source)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  stat_summary(fun.y=mean, geom="point", shape=20, position = position_dodge(width = 0.9), aes(group = source, color = "Per motif rate"), size = 5) + 
  geom_point(data = raw_error_rates_ONT_moderate, aes(x = type, y = raw_smm_rate, group = source, color = "Aggregate rate"), 
             shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = 9) +
  ylab("SNM per nucleotide") + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(labels = number_format(accuracy = 0.001,decimal.mark = '.')) +
  ggtitle("ONT") +
  scale_colour_manual(values = c("Aggregate rate" = "#feb24c"  , "Per motif rate" = "#1a1a1a"  )) +
  guides(colour = guide_legend(override.aes = list(shape = c(17, 20))),
         fill = guide_legend(override.aes = list(linetype = c(0, 0),
                                                 shape = c(NA, NA)))) + 
  theme( legend.position = "bottom",
         plot.title = element_text(hjust = 0.5),
         legend.title = element_blank(),
         axis.title.x = element_blank(),
         axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1)) -> ONT_smm_boxplots
ONT_smm_ymax_moderate <- layer_scales(ONT_smm_boxplots)$y$range$range[2]

ONT_ultrastringent_filter %>% 
  ggplot(aes(y = smm_rate, x = type, fill = source)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  stat_summary(fun.y=mean, geom="point", shape=20, position = position_dodge(width = 0.9), aes(group = source, color = "Per motif rate"), size = 5) + 
  geom_point(data = raw_error_rates_ONT_ultrastringent, aes(x = type, y = raw_smm_rate, group = source, color = "Aggregate rate"), 
             shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = 9) +
  ylab("smm per nucleotide") + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(labels = number_format(accuracy = 0.001,decimal.mark = '.')) +
  coord_cartesian(ylim =  c(0, ONT_smm_ymax_moderate)) + 
  scale_colour_manual(values = c("Aggregate rate" = "#feb24c"  , "Per motif rate" = "#1a1a1a"  )) +
  guides(colour = guide_legend(override.aes = list(shape = c(17, 20))),
         fill = guide_legend(override.aes = list(linetype = c(0, 0),
                                                 shape = c(NA, NA)))) + 
  ggtitle("ONT") +
  theme( legend.position = "bottom",
         plot.title = element_text(hjust = 0.5),
         legend.title = element_blank(),
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1)) -> ONT_ultrastringent_smm_boxplots
#####

SNM_boxplots <-  wrap_elements(grid::textGrob('Moderate filter')) + wrap_elements(grid::textGrob('Stringent filter')) + 
  Illumina_smm_boxplots + Illumina_ultrastringent_smm_boxplots +  
  HiFi_smm_boxplots + HiFi_ultrastringent_smm_boxplots +
  ONT_smm_boxplots + ONT_ultrastringent_smm_boxplots + plot_layout(nrow = 4, heights = c(1,4,4,4),
                                                                   guides = 'collect') & theme(legend.position = "bottom") 

ggsave(paste0("smm_moderate_ultrastringent_filter.", format(Sys.time(), "%Y-%m-%d"), ".png" ), height= 150, width = 185,units = "mm", SNM_boxplots)

#################################################
# Single-nucleotide mismatches - Motif subsection boxplots

## Illumina
#####
Illumina_motif_subsections %>%
  ggplot(aes(y = smm_rate, x = type, fill = simple_section)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  geom_point(data = Illumina_motif_subsections_raw_error_rates, aes(x = type, y = per_motif_mean_smm, group = simple_section, color = "Per motif rate"), 
             size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  geom_point(data = Illumina_motif_subsections_raw_error_rates, aes(x = type, y = aggregate_mean_smm, group = simple_section, color = "Aggregate rate"), 
             shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = base_size) +
  scale_fill_manual(values = colors_simple_section, labels = c("Repeat / stem", "Spacer / loop")) + 
  scale_colour_manual(values = c("Aggregate rate" = "#feb24c"  , "Per motif rate" = "#1a1a1a"  )) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001,decimal.mark = '.')) +
  guides(colour = guide_legend(override.aes = list(shape = c(17, 20))),
         fill = guide_legend(override.aes = list(linetype = c(0, 0),
                                                 shape = c(NA, NA)))) + 
  ggtitle("Illumina") +
  ylab("SNM per nucleotide") + 
  theme( legend.position = "bottom",
         plot.title = element_text(hjust = 0.5),
         axis.title.x = element_blank(),
         legend.title = element_blank(),
         axis.text.x = element_blank()) -> Illumina_smm_subsections_moderate
Illumina_smm_subsections_ymax_moderate <- layer_scales(Illumina_smm_subsections_moderate)$y$range$range[2]

Illumina_motif_subsections_ultrastringent %>%
  ggplot(aes(y = smm_rate, x = type, fill = simple_section)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  geom_point(data = Illumina_motif_subsections_ultrastringent_raw_error_rates, aes(x = type, y = per_motif_mean_smm, group = simple_section, color = "Per motif rate"), 
             size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  geom_point(data = Illumina_motif_subsections_ultrastringent_raw_error_rates, aes(x = type, y = aggregate_mean_smm, group = simple_section, color = "Aggregate rate"), 
             shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = base_size) +
  scale_fill_manual(values = colors_simple_section, labels = c("Repeat / stem", "Spacer / loop")) + 
  scale_colour_manual(values = c("Aggregate rate" = "#feb24c"  , "Per motif rate" = "#1a1a1a"  )) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001,decimal.mark = '.')) +
  coord_cartesian(ylim =c(0, Illumina_smm_subsections_ymax_moderate)) + 
  guides(colour = guide_legend(override.aes = list(shape = c(17, 20))),
         fill = guide_legend(override.aes = list(linetype = c(0, 0),
                                                 shape = c(NA, NA)))) + 
  ggtitle("Illumina") +
  ylab("SNM per nucleotide") + 
  theme( legend.position = "bottom",
         plot.title = element_text(hjust = 0.5),
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         legend.title = element_blank(),
         axis.text.x = element_blank()) -> Illumina_smm_subsections_ultrastringent
#####

## HiFi 
#####
HiFi_motif_subsections %>%
  ggplot(aes(y = smm_rate, x = type, fill = simple_section)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  geom_point(data = HiFi_motif_subsections_raw_error_rates, aes(x = type, y = per_motif_mean_smm, group = simple_section, color = "Per motif rate"), 
             size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  geom_point(data = HiFi_motif_subsections_raw_error_rates, aes(x = type, y = aggregate_mean_smm, group = simple_section, color = "Aggregate rate"), 
             shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = base_size) +
  scale_fill_manual(values = colors_simple_section, labels = c("Repeat / stem", "Spacer / loop")) + 
  scale_colour_manual(values = c("Aggregate rate" = "#feb24c"  , "Per motif rate" = "#1a1a1a"  )) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001,decimal.mark = '.') )+
  guides(colour = guide_legend(override.aes = list(shape = c(17, 20))),
         fill = guide_legend(override.aes = list(linetype = c(0, 0),
                                                 shape = c(NA, NA)))) + 
  ggtitle("HiFi") +
  ylab("SNM per nucleotide") + 
  theme( legend.position = "bottom",
         plot.title = element_text(hjust = 0.5),
         axis.title.x = element_blank(),
         legend.title = element_blank(),
         axis.text.x = element_blank()) -> HiFi_smm_subsections_moderate
HiFi_smm_subsections_ymax_moderate <- layer_scales(HiFi_smm_subsections_moderate)$y$range$range[2]

HiFi_motif_subsections_ultrastringent %>%
  ggplot(aes(y = smm_rate, x = type, fill = simple_section)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  geom_point(data = HiFi_motif_subsections_ultrastringent_raw_error_rates, aes(x = type, y = per_motif_mean_smm, group = simple_section, color = "Per motif rate"), 
             size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  geom_point(data = HiFi_motif_subsections_ultrastringent_raw_error_rates, aes(x = type, y = aggregate_mean_smm, group = simple_section, color = "Aggregate rate"), 
             shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = base_size) +
  scale_fill_manual(values = colors_simple_section, labels = c("Repeat / stem", "Spacer / loop")) + 
  scale_colour_manual(values = c("Aggregate rate" = "#feb24c"  , "Per motif rate" = "#1a1a1a"  )) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001,decimal.mark = '.')) +
  coord_cartesian(ylim =  c(0, HiFi_smm_subsections_ymax_moderate)) + 
  guides(colour = guide_legend(override.aes = list(shape = c(17, 20))),
         fill = guide_legend(override.aes = list(linetype = c(0, 0),
                                                 shape = c(NA, NA)))) + 
  ggtitle("HiFi") +
  ylab("SNM per nucleotide") + 
  theme( legend.position = "bottom",
         plot.title = element_text(hjust = 0.5),
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         legend.title = element_blank(),
         axis.text.x = element_blank()) -> HiFi_smm_subsections_ultrastringent
#####

## ONT
#####
ONT_motif_subsections %>%
  ggplot(aes(y = smm_rate, x = type, fill = simple_section)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  geom_point(data = ONT_motif_subsections_raw_error_rates, aes(x = type, y = per_motif_mean_smm, group = simple_section, color = "Per motif rate"), 
             size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  geom_point(data = ONT_motif_subsections_raw_error_rates, aes(x = type, y = aggregate_mean_smm, group = simple_section, color = "Aggregate rate"), 
             shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = base_size) +
  scale_fill_manual(values = colors_simple_section, labels = c("Repeat / stem", "Spacer / loop")) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001,decimal.mark = '.')) +
  scale_colour_manual(values = c("Aggregate rate" = "#feb24c"  , "Per motif rate" = "#1a1a1a"  )) +
  guides(colour = guide_legend(override.aes = list(shape = c(17, 20))),
         fill = guide_legend(override.aes = list(linetype = c(0, 0),
                                                 shape = c(NA, NA)))) + 
  ggtitle("ONT") +
  ylab("SNM per nucleotide") + 
  theme( legend.position = "bottom",
         plot.title = element_text(hjust = 0.5),
         axis.title.x = element_blank(),
         legend.title = element_blank(),
         axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1)) -> ONT_smm_subsections_moderate
ONT_smm_subsections_ymax_moderate <- layer_scales(ONT_smm_subsections_moderate)$y$range$range[2]

ONT_motif_subsections_ultrastringent %>%
  ggplot(aes(y = smm_rate, x = type, fill = simple_section)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  geom_point(data = ONT_motif_subsections_ultrastringent_raw_error_rates, aes(x = type, y = per_motif_mean_smm, group = simple_section, color = "Per motif rate"), 
             size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  geom_point(data = ONT_motif_subsections_ultrastringent_raw_error_rates, aes(x = type, y = aggregate_mean_smm, group = simple_section, color = "Aggregate rate"), 
             shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = base_size) +
  scale_fill_manual(values = colors_simple_section, labels = c("Repeat / stem", "Spacer / loop")) + 
  scale_colour_manual(values = c("Aggregate rate" = "#feb24c"  , "Per motif rate" = "#1a1a1a"  )) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001,decimal.mark = '.')) +
  coord_cartesian(ylim =  c(0, ONT_smm_subsections_ymax_moderate)) + 
  guides(colour = guide_legend(override.aes = list(shape = c(17, 20))),
         fill = guide_legend(override.aes = list(linetype = c(0, 0),
                                                 shape = c(NA, NA)))) + 
  ggtitle("ONT") +
  ylab("SNM per nucleotide") + 
  theme( legend.position = "bottom",
         plot.title = element_text(hjust = 0.5),
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         legend.title = element_blank(),
         axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1)) -> ONT_smm_subsections_ultrastringent
#####

figure_smm_subsections<- wrap_elements(grid::textGrob('Moderate filter')) + wrap_elements(grid::textGrob('Stringent filter')) +
  Illumina_smm_subsections_moderate + Illumina_smm_subsections_ultrastringent + 
  HiFi_smm_subsections_moderate + HiFi_smm_subsections_ultrastringent + 
  ONT_smm_subsections_moderate + ONT_smm_subsections_ultrastringent +  
  plot_layout(nrow = 4, heights = c(1,4,4,4), guides = 'collect') & theme(legend.position = "bottom") 

ggsave(paste0("Motif_subsections_simplified_filter_comparison_smm.", format(Sys.time(), "%Y-%m-%d"), ".png" ), height= 150, width = 185,units = "mm",figure_smm_subsections)

#################################################
# Deletion errors

## Illumina
#####
Illumina_curated %>% 
  ggplot(aes(y = del_rate, x = type, fill = source)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  stat_summary(fun.y=mean, geom="point", shape=20, position = position_dodge(width = 0.9), aes(group = source, color = "Per motif rate"), size = 5) + 
  geom_point(data = raw_error_rates_Illumina_moderate, aes(x = type, y = raw_del_rate, group = source, color = "Aggregate rate"), 
             shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = 9) +
  ylab("DEL per nucleotide") + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.0001,decimal.mark = '.')) +
  ggtitle("Illumina") +
  scale_colour_manual(values = c("Aggregate rate" = "#feb24c"  , "Per motif rate" = "#1a1a1a"  )) +
  guides(colour = guide_legend(override.aes = list(shape = c(17, 20))),
         fill = guide_legend(override.aes = list(linetype = c(0, 0),
                                                 shape = c(NA, NA)))) + 
  theme( legend.position = "bottom",
         plot.title = element_text(hjust = 0.5),
         legend.title = element_blank(),
         axis.title.x = element_blank(),
         axis.text.x = element_blank()) -> Illumina_del_boxplots
Illumina_del_ymax_moderate <- layer_scales(Illumina_del_boxplots)$y$range$range[2]

Illumina_ultrastringent_filter %>% 
  ggplot(aes(y = del_rate, x = type, fill = source)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  stat_summary(fun.y=mean, geom="point", shape=20, position = position_dodge(width = 0.9), aes(group = source, color = "Per motif rate"), size = 5) + 
  geom_point(data = raw_error_rates_Illumina_ultrastringent, aes(x = type, y = raw_del_rate, group = source, color = "Aggregate rate"), 
             shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = 9) +
  ylab("DEL per nucleotide") + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.0001,decimal.mark = '.')) +
  coord_cartesian(ylim =  c(0, Illumina_del_ymax_moderate)) + 
  scale_colour_manual(values = c("Aggregate rate" = "#feb24c"  , "Per motif rate" = "#1a1a1a"  )) +
  guides(colour = guide_legend(override.aes = list(shape = c(17, 20))),
         fill = guide_legend(override.aes = list(linetype = c(0, 0),
                                                 shape = c(NA, NA)))) + 
  ggtitle("Illumina") +
  theme( legend.position = "bottom",
         plot.title = element_text(hjust = 0.5),
         legend.title = element_blank(),
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.text.x = element_blank()) -> Illumina_ultrastringent_del_boxplots

#####

## HiFi 
#####
HiFi_curated %>% 
  ggplot(aes(y = del_rate, x = type, fill = source)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  stat_summary(fun.y=mean, geom="point", shape=20, position = position_dodge(width = 0.9), aes(group = source, color = "Per motif rate"), size = 5) + 
  geom_point(data = raw_error_rates_HiFi_moderate, aes(x = type, y = raw_del_rate, group = source, color = "Aggregate rate"), 
             shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = 9) +
  ylab("DEL per nucleotide") + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.0001,decimal.mark = '.')) +
  ggtitle("HiFi") +
  scale_colour_manual(values = c("Aggregate rate" = "#feb24c"  , "Per motif rate" = "#1a1a1a"  )) +
  guides(colour = guide_legend(override.aes = list(shape = c(17, 20))),
         fill = guide_legend(override.aes = list(linetype = c(0, 0),
                                                 shape = c(NA, NA)))) + 
  theme( legend.position = "bottom",
         plot.title = element_text(hjust = 0.5),
         legend.title = element_blank(),
         axis.title.x = element_blank(),
         axis.text.x = element_blank()) -> HiFi_del_boxplots
HiFi_del_ymax_moderate <- layer_scales(HiFi_del_boxplots)$y$range$range[2]

HiFi_ultrastringent_filter %>% 
  ggplot(aes(y = del_rate, x = type, fill = source)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  stat_summary(fun.y=mean, geom="point", shape=20, position = position_dodge(width = 0.9), aes(group = source, color = "Per motif rate"), size = 5) + 
  geom_point(data = raw_error_rates_HiFi_ultrastringent, aes(x = type, y = raw_del_rate, group = source, color = "Aggregate rate"), 
             shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = 9) +
  ylab("DEL per nucleotide") + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.0001,decimal.mark = '.')) +
  coord_cartesian(ylim =  c(0, HiFi_del_ymax_moderate)) + 
  scale_colour_manual(values = c("Aggregate rate" = "#feb24c"  , "Per motif rate" = "#1a1a1a"  )) +
  guides(colour = guide_legend(override.aes = list(shape = c(17, 20))),
         fill = guide_legend(override.aes = list(linetype = c(0, 0),
                                                 shape = c(NA, NA)))) + 
  ggtitle("HiFi") +
  theme( legend.position = "bottom",
         plot.title = element_text(hjust = 0.5),
         legend.title = element_blank(),
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.text.x = element_blank()) -> HiFi_ultrastringent_del_boxplots
#####

## ONT 
#####
ONT_curated %>% 
  ggplot(aes(y = del_rate, x = type, fill = source)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  stat_summary(fun.y=mean, geom="point", shape=20, position = position_dodge(width = 0.9), aes(group = source, color = "Per motif rate"), size = 5) + 
  geom_point(data = raw_error_rates_ONT_moderate, aes(x = type, y = raw_del_rate, group = source, color = "Aggregate rate"), 
             shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = 9) +
  ylab("DEL per nucleotide") + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.0001,decimal.mark = '.')) +
  ggtitle("ONT") +
  scale_colour_manual(values = c("Aggregate rate" = "#feb24c"  , "Per motif rate" = "#1a1a1a"  )) +
  guides(colour = guide_legend(override.aes = list(shape = c(17, 20))),
         fill = guide_legend(override.aes = list(linetype = c(0, 0),
                                                 shape = c(NA, NA)))) + 
  theme( legend.position = "bottom",
         plot.title = element_text(hjust = 0.5),
         legend.title = element_blank(),
         axis.title.x = element_blank(),
         axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1)) -> ONT_del_boxplots
ONT_del_ymax_moderate <- layer_scales(ONT_del_boxplots)$y$range$range[2]

ONT_ultrastringent_filter %>% 
  ggplot(aes(y = del_rate, x = type, fill = source)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  stat_summary(fun.y=mean, geom="point", shape=20, position = position_dodge(width = 0.9), aes(group = source, color = "Per motif rate"), size = 5) + 
  geom_point(data = raw_error_rates_ONT_ultrastringent, aes(x = type, y = raw_del_rate, group = source, color = "Aggregate rate"), 
             shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = 9) +
  ylab("DEL per nucleotide") + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.0001,decimal.mark = '.')) +
  coord_cartesian(ylim =  c(0, ONT_del_ymax_moderate)) + 
  scale_colour_manual(values = c("Aggregate rate" = "#feb24c"  , "Per motif rate" = "#1a1a1a"  )) +
  guides(colour = guide_legend(override.aes = list(shape = c(17, 20))),
         fill = guide_legend(override.aes = list(linetype = c(0, 0),
                                                 shape = c(NA, NA)))) + 
  ggtitle("ONT") +
  theme( legend.position = "bottom",
         plot.title = element_text(hjust = 0.5),
         legend.title = element_blank(),
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1)) -> ONT_ultrastringent_del_boxplots
#####

DEL_boxplots <-  wrap_elements(grid::textGrob('Moderate filter')) + wrap_elements(grid::textGrob('Stringent filter')) + 
  Illumina_del_boxplots + Illumina_ultrastringent_del_boxplots +  
  HiFi_del_boxplots + HiFi_ultrastringent_del_boxplots +
  ONT_del_boxplots + ONT_ultrastringent_del_boxplots + plot_layout(nrow = 4, heights = c(1,4,4,4),
                                                                   guides = 'collect') & theme(legend.position = "bottom") 

ggsave(paste0("Deletions_moderate_ultrastringent_filter.", format(Sys.time(), "%Y-%m-%d"), ".png" ), height= 150, width = 185,units = "mm", DEL_boxplots)

#################################################
# Insertion errors

## Illumina
#####
Illumina_curated %>% 
  ggplot(aes(y = ins_rate, x = type, fill = source)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  stat_summary(fun.y=mean, geom="point", shape=20, position = position_dodge(width = 0.9), aes(group = source, color = "Per motif rate"), size = 5) + 
  geom_point(data = raw_error_rates_Illumina_moderate, aes(x = type, y = raw_ins_rate, group = source, color = "Aggregate rate"), 
             shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = 9) +
  ylab("INS per nucleotide") + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.0001,decimal.mark = '.')) +
  ggtitle("Illumina") +
  scale_colour_manual(values = c("Aggregate rate" = "#feb24c"  , "Per motif rate" = "#1a1a1a"  )) +
  guides(colour = guide_legend(override.aes = list(shape = c(17, 20))),
         fill = guide_legend(override.aes = list(linetype = c(0, 0),
                                                 shape = c(NA, NA)))) + 
  theme( legend.position = "bottom",
         plot.title = element_text(hjust = 0.5),
         legend.title = element_blank(),
         axis.title.x = element_blank(),
         axis.text.x = element_blank()) -> Illumina_ins_boxplots
Illumina_ins_ymax_moderate <- layer_scales(Illumina_ins_boxplots)$y$range$range[2]

Illumina_ultrastringent_filter %>% 
  ggplot(aes(y = ins_rate, x = type, fill = source)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  stat_summary(fun.y=mean, geom="point", shape=20, position = position_dodge(width = 0.9), aes(group = source, color = "Per motif rate"), size = 5) + 
  geom_point(data = raw_error_rates_Illumina_ultrastringent, aes(x = type, y = raw_ins_rate, group = source, color = "Aggregate rate"), 
             shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = 9) +
  ylab("ins per nucleotide") + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.0001,decimal.mark = '.')) +
  coord_cartesian(ylim = c(0, Illumina_ins_ymax_moderate)) + 
  scale_colour_manual(values = c("Aggregate rate" = "#feb24c"  , "Per motif rate" = "#1a1a1a"  )) +
  guides(colour = guide_legend(override.aes = list(shape = c(17, 20))),
         fill = guide_legend(override.aes = list(linetype = c(0, 0),
                                                 shape = c(NA, NA)))) + 
  ggtitle("Illumina") +
  theme( legend.position = "bottom",
         plot.title = element_text(hjust = 0.5),
         legend.title = element_blank(),
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.text.x = element_blank()) -> Illumina_ultrastringent_ins_boxplots
#####

## HiFi 
#####
HiFi_curated %>% 
  ggplot(aes(y = ins_rate, x = type, fill = source)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  stat_summary(fun.y=mean, geom="point", shape=20, position = position_dodge(width = 0.9), aes(group = source, color = "Per motif rate"), size = 5) + 
  geom_point(data = raw_error_rates_HiFi_moderate, aes(x = type, y = raw_ins_rate, group = source, color = "Aggregate rate"), 
             shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = 9) +
  ylab("INS per nucleotide") + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.0001,decimal.mark = '.')) +
  ggtitle("HiFi") +
  scale_colour_manual(values = c("Aggregate rate" = "#feb24c"  , "Per motif rate" = "#1a1a1a"  )) +
  guides(colour = guide_legend(override.aes = list(shape = c(17, 20))),
         fill = guide_legend(override.aes = list(linetype = c(0, 0),
                                                 shape = c(NA, NA)))) + 
  theme( legend.position = "bottom",
         plot.title = element_text(hjust = 0.5),
         legend.title = element_blank(),
         axis.title.x = element_blank(),
         axis.text.x = element_blank()) -> HiFi_ins_boxplots
HiFi_ins_ymax_moderate <- layer_scales(HiFi_ins_boxplots)$y$range$range[2]

HiFi_ultrastringent_filter %>% 
  ggplot(aes(y = ins_rate, x = type, fill = source)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  stat_summary(fun.y=mean, geom="point", shape=20, position = position_dodge(width = 0.9), aes(group = source, color = "Per motif rate"), size = 5) + 
  geom_point(data = raw_error_rates_HiFi_ultrastringent, aes(x = type, y = raw_ins_rate, group = source, color = "Aggregate rate"), 
             shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = 9) +
  ylab("ins per nucleotide") + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.0001,decimal.mark = '.')) +
  coord_cartesian(ylim =  c(0, HiFi_ins_ymax_moderate)) + 
  scale_colour_manual(values = c("Aggregate rate" = "#feb24c"  , "Per motif rate" = "#1a1a1a"  )) +
  guides(colour = guide_legend(override.aes = list(shape = c(17, 20))),
         fill = guide_legend(override.aes = list(linetype = c(0, 0),
                                                 shape = c(NA, NA)))) + 
  ggtitle("HiFi") +
  theme( legend.position = "bottom",
         plot.title = element_text(hjust = 0.5),
         legend.title = element_blank(),
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.text.x = element_blank()) -> HiFi_ultrastringent_ins_boxplots

#####

## ONT 
#####
ONT_curated %>% 
  ggplot(aes(y = ins_rate, x = type, fill = source)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  stat_summary(fun.y=mean, geom="point", shape=20, position = position_dodge(width = 0.9), aes(group = source, color = "Per motif rate"), size = 5) + 
  geom_point(data = raw_error_rates_ONT_moderate, aes(x = type, y = raw_ins_rate, group = source, color = "Aggregate rate"), 
             shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = 9) +
  ylab("INS per nucleotide") + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.0001,decimal.mark = '.')) +
  ggtitle("ONT") +
  scale_colour_manual(values = c("Aggregate rate" = "#feb24c"  , "Per motif rate" = "#1a1a1a"  )) +
  guides(colour = guide_legend(override.aes = list(shape = c(17, 20))),
         fill = guide_legend(override.aes = list(linetype = c(0, 0),
                                                 shape = c(NA, NA)))) + 
  theme( legend.position = "bottom",
         plot.title = element_text(hjust = 0.5),
         legend.title = element_blank(),
         axis.title.x = element_blank(),
         axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1)) -> ONT_ins_boxplots
ONT_ins_ymax_moderate <- layer_scales(ONT_ins_boxplots)$y$range$range[2]

ONT_ultrastringent_filter %>% 
  ggplot(aes(y = ins_rate, x = type, fill = source)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  stat_summary(fun.y=mean, geom="point", shape=20, position = position_dodge(width = 0.9), aes(group = source, color = "Per motif rate"), size = 5) + 
  geom_point(data = raw_error_rates_ONT_ultrastringent, aes(x = type, y = raw_ins_rate, group = source, color = "Aggregate rate"), 
             shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = 9) +
  ylab("ins per nucleotide") + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.0001,decimal.mark = '.')) +
  coord_cartesian(ylim =  c(0, ONT_ins_ymax_moderate)) + 
  scale_colour_manual(values = c("Aggregate rate" = "#feb24c"  , "Per motif rate" = "#1a1a1a"  )) +
  guides(colour = guide_legend(override.aes = list(shape = c(17, 20))),
         fill = guide_legend(override.aes = list(linetype = c(0, 0),
                                                 shape = c(NA, NA)))) + 
  ggtitle("ONT") +
  theme( legend.position = "bottom",
         plot.title = element_text(hjust = 0.5),
         legend.title = element_blank(),
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1)) -> ONT_ultrastringent_ins_boxplots
#####

INS_boxplots <-  wrap_elements(grid::textGrob('Moderate filter')) + wrap_elements(grid::textGrob('Stringent filter')) + 
  Illumina_ins_boxplots + Illumina_ultrastringent_ins_boxplots +  
  HiFi_ins_boxplots + HiFi_ultrastringent_ins_boxplots +
  ONT_ins_boxplots + ONT_ultrastringent_ins_boxplots + plot_layout(nrow = 4, heights = c(1,4,4,4),
                                                                   guides = 'collect') & theme(legend.position = "bottom") 

ggsave(paste0("Insertions_moderate_ultrastringent_filter.", format(Sys.time(), "%Y-%m-%d"), ".png" ), height= 150, width = 185,units = "mm", INS_boxplots)

#################################################
# Heatmaps

# Fold-change table heatmap
#####

rbind.data.frame(fold_change_table_moderate, fold_change_table_ultrastringent) -> fold_change_table_combined

fold_change_table_combined$tech <- as.factor(fold_change_table_combined$tech)
fold_change_table_combined$tech <- with(fold_change_table_combined, relevel(tech, "Illumina")) 

fold_change_table_combined %>% 
  filter(variable %in% c("per_motif_del_fold_change", "per_motif_ins_fold_change","per_motif_smm_fold_change", "depth")) %>%
  mutate(error_type = ifelse(grepl("smm", variable), "smm", 
                             ifelse(grepl("del", variable), "del", "ins"))) %>%
  merge(., molten_p_value_table, by = c("type", "tech", "filter_level", "error_type")) %>% 
  mutate(fontface = ifelse(value.y <= 0.05, 2, 1)) %>% 
  select(type, tech, filter_level, error_type, variable.x, value.x, fontface) %>% 
  rename(variable = variable.x, value = value.x) %>%
  select(-error_type) -> fold_change_table_combined_with_fontface 

# Per motif mean fold change
fold_change_table_combined %>%
  filter(variable %in% c("aggregate_smm_fold_change", "aggregate_ins_fold_change", "aggregate_del_fold_change")) %>%
  mutate(fontface = 1) %>% 
  rbind.data.frame(., fold_change_table_combined_with_fontface) -> fold_change_table_combined_merged_with_fontface

fold_change_table_combined_merged_with_fontface %>%
  mutate(fill_color = log(value)) %>%
  filter(variable %in% c( "per_motif_smm_fold_change" )) %>% 
  ggplot(aes(y = variable ,x = type)) + 
  geom_tile(aes(fill = fill_color)) +
  scale_fill_gradient2( low = "#1A85FF", high =  "#D41159") + 
  geom_text(aes(y = variable, x = type,fontface = fontface, label=sprintf("%0.2f", round(value, 2))))  + 
  facet_grid(tech~filter_level) + 
  ylab("Fold change per motif SNM error rate") + 
  theme(axis.title.x = element_blank() ,
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=0.9),
        legend.position = "none") 

ggsave(paste0("smm_fold_changes_heatmap_filter_comparison.", format(Sys.time(), "%Y-%m-%d"), ".png" ), height= 100, width = 185,units = "mm")

fold_change_table_combined_merged_with_fontface %>%
  mutate(fill_color = log(value)) %>%
  filter(variable %in% c( "per_motif_ins_fold_change" )) %>% 
  ggplot(aes(y = variable ,x = type)) + 
  geom_tile(aes(fill = fill_color)) +
  scale_fill_gradient2( low = "#1A85FF", high =  "#D41159") + 
  geom_text(aes(y = variable, x = type,fontface = fontface, label= sprintf("%0.2f", round(value, 2))))  + 
  facet_grid(tech~filter_level) + 
  ylab("Fold change per motif insertion error rate") + 
  theme(axis.title.x = element_blank() ,
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=0.9),
        legend.position = "none") 

ggsave(paste0("ins_fold_changes_heatmap_filter_comparison.", format(Sys.time(), "%Y-%m-%d"), ".png" ), height= 100, width = 185,units = "mm")

fold_change_table_combined_merged_with_fontface  %>%
  mutate(fill_color = log(value)) %>%
  filter(variable %in% c( "per_motif_del_fold_change" )) %>% 
  ggplot(aes(y = variable ,x = type)) + 
  geom_tile(aes(fill = fill_color)) +
  scale_fill_gradient2( low = "#1A85FF", high =  "#D41159") + 
  geom_text(aes(y = variable, x = type,fontface = fontface, label= sprintf("%0.2f",round(value, 2))))  + 
  facet_grid(tech~filter_level) + 
  ylab("Fold change per motif deletion error rate") + 
  theme(axis.title.x = element_blank() ,
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=0.9),
        legend.position = "none") 

ggsave(paste0("del_fold_changes_heatmap_filter_comparison.", format(Sys.time(), "%Y-%m-%d"), ".png" ), height= 100, width = 185,units = "mm")
#####

# Fold change table heatmap motif subsections
#####
p_value_table_motif_subsections_filter_levels %>%
  select(type, smm_p_adjusted, tech, filter_level ) %>% 
  mutate(error_type = "smm") -> temp

rbind.data.frame(moderate_filter_fold_change_simple_subsections_combined,
                 fold_change_simple_subsections_ultrastringent_combined) %>% 
  filter(error_type == "smm") %>%  
  select(type, tech, error_type, filter_level,fold_change_slrs_per_motif_mean ) %>% 
  merge(., temp, by = c("type", "tech", "error_type", "filter_level")) %>% 
  mutate(fontface = ifelse(smm_p_adjusted >=0.05,1, 2)) %>% 
  mutate(fill_color = log(fold_change_slrs_per_motif_mean)) %>%
  ggplot(aes(y = error_type ,x = type)) + 
  geom_tile(aes(fill = fill_color)) +
  scale_fill_gradient2( low = "#1A85FF", high =  "#D41159") + 
  geom_text(aes(y = error_type, x = type,fontface = fontface,label= sprintf("%0.2f", round(fold_change_slrs_per_motif_mean, 2))))  + 
  ylab("Fold change per motif SNM rate") + 
  facet_grid(relevel(as.factor(tech), "Illumina")~filter_level) + 
  theme(axis.title.x = element_blank() ,
        axis.text.y = element_blank() , 
        axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1),
        legend.position = "none") 

ggsave(paste0("motif_subsections_fold_change_heatmap_with_pvalues.", format(Sys.time(), "%Y-%m-%d"), ".png" ), height= 100, width = 185,units = "mm")
#####

# Heatmap aggreagte error rates
#####

fold_change_table_combined %>%
  filter(variable %in% c("aggregate_smm_fold_change", "aggregate_ins_fold_change", "aggregate_del_fold_change")) %>%
  mutate(fill_color = log(value)) %>%
  ggplot(aes(y = variable ,x = type)) + 
  geom_tile(aes(fill = fill_color)) +
  scale_fill_gradient2( low = "#1A85FF", high =  "#D41159") + 
  geom_text(aes(y = variable, x = type,label=  sprintf("%0.2f", round(value, 2))))  + 
  scale_y_discrete(limits = c("aggregate_ins_fold_change",  "aggregate_del_fold_change", "aggregate_smm_fold_change"),
                   labels = c("aggregate_ins_fold_change" = "INS",  "aggregate_del_fold_change" = "DEL", "aggregate_smm_fold_change" = "SNM")) +
  facet_grid(tech~filter_level) + 
  theme(axis.title = element_blank() ,
        axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1),
        legend.position = "none")
ggsave(paste0("aggregate_mean_fold_changes_heatmap_filter_comparison.", format(Sys.time(), "%Y-%m-%d"), ".png" ), height= 150, width = 185,units = "mm")

#####

#################################################
# Read depth per bp

# Boxplots

## Illumina
#####
Illumina_curated %>% 
  ggplot(aes(y = depth_per_bp, x = type, fill = source)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  stat_summary(fun.y=mean, geom="point", shape=20, position = position_dodge(width = 0.9), aes(group = source), size = 5) + 
  geom_point(data = raw_depth_per_bp_Illumina, aes(x = type, y = raw_depth_per_bp, group = source), 
             color = "orange", shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = 9) +
  ylab("Read depth per bp") + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(labels = scales::comma) +
  ggtitle("Illumina") +
  theme( legend.position = "none",
         plot.title = element_text(hjust = 0.5),
         legend.title = element_blank(),
         axis.title.x = element_blank(),
         axis.text.x = element_blank()) -> Illumina_moderate_depth_boxplots

Illumina_ultrastringent_filter %>% 
  ggplot(aes(y = depth_per_bp, x = type, fill = source)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  stat_summary(fun.y=mean, geom="point", shape=20, position = position_dodge(width = 0.9), aes(group = source), size = 5) + 
  geom_point(data = raw_depth_per_bp_Illumina, aes(x = type, y = raw_depth_per_bp, group = source), 
             color = "orange", shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = 9) +
  ylab("Read depth per bp") + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(labels = scales::comma) +
  ggtitle("Illumina") +
  theme( legend.position = "none",
         plot.title = element_text(hjust = 0.5),
         legend.title = element_blank(),
         axis.title.x = element_blank(),
         axis.text.x = element_blank()) -> Illumina_ultrastringent_depth_boxplots

#####

## HiFi
#####
HiFi_curated %>% 
  ggplot(aes(y = depth_per_bp, x = type, fill = source)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  stat_summary(fun.y=mean, geom="point", shape=20, position = position_dodge(width = 0.9), aes(group = source), size = 5) + 
  geom_point(data = raw_depth_per_bp_HiFi, aes(x = type, y = raw_depth_per_bp, group = source), 
             color = "orange", shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = 9) +
  ylab("Read depth per bp") + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(labels = scales::comma) +
  ggtitle("HiFi") +
  theme( legend.position = "none",
         plot.title = element_text(hjust = 0.5),
         legend.title = element_blank(),
         axis.title.x = element_blank(),
         axis.text.x = element_blank()) -> HiFi_moderate_depth_boxplots

HiFi_ultrastringent_filter %>% 
  ggplot(aes(y = depth_per_bp, x = type, fill = source)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  stat_summary(fun.y=mean, geom="point", shape=20, position = position_dodge(width = 0.9), aes(group = source), size = 5) + 
  geom_point(data = raw_depth_per_bp_HiFi, aes(x = type, y = raw_depth_per_bp, group = source), 
             color = "orange", shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = 9) +
  ylab("Read depth per bp") + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(labels = scales::comma) +
  ggtitle("HiFi") +
  theme( legend.position = "none",
         plot.title = element_text(hjust = 0.5),
         legend.title = element_blank(),
         axis.title.x = element_blank(),
         axis.text.x = element_blank()) -> HiFi_ultrastringent_depth_boxplots

#####

## ONT
#####
ONT_curated %>% 
  ggplot(aes(y = depth_per_bp, x = type, fill = source)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  stat_summary(fun.y=mean, geom="point", shape=20, position = position_dodge(width = 0.9), aes(group = source), size = 5) + 
  geom_point(data = raw_depth_per_bp_ONT, aes(x = type, y = raw_depth_per_bp, group = source), 
             color = "orange", shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = 9) +
  ylab("Read depth per bp") + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(labels = scales::comma) +
  ggtitle("ONT") +
  theme( legend.position = "none",
         plot.title = element_text(hjust = 0.5),
         legend.title = element_blank(),
         axis.title.x = element_blank(),
         axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1) ) -> ONT_moderate_depth_boxplots

ONT_ultrastringent_filter %>% 
  ggplot(aes(y = depth_per_bp, x = type, fill = source)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  stat_summary(fun.y=mean, geom="point", shape=20, position = position_dodge(width = 0.9), aes(group = source), size = 5) + 
  geom_point(data = raw_depth_per_bp_ONT, aes(x = type, y = raw_depth_per_bp, group = source), 
             color = "orange", shape = 17, size = 2.5,  stat = "identity",position = position_dodge(width = 0.9)) + 
  theme_classic(base_size = 9) +
  ylab("Read depth per bp") + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(labels = scales::comma) +
  ggtitle("ONT") +
  theme( legend.position = "none",
         plot.title = element_text(hjust = 0.5),
         legend.title = element_blank(),
         axis.title.x = element_blank(),
         axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1)) -> ONT_ultrastringent_depth_boxplots

#####

depth_boxplots <- (Illumina_moderate_depth_boxplots | Illumina_ultrastringent_depth_boxplots) / 
  (HiFi_moderate_depth_boxplots | HiFi_ultrastringent_depth_boxplots) /
  ( ONT_moderate_depth_boxplots | ONT_ultrastringent_depth_boxplots)

ggsave(paste0("Figure_5_depth_filter_comparison.", format(Sys.time(), "%Y-%m-%d"), ".png" ), height= 150, width = 185,units = "mm", depth_boxplots)

# Heatmaps
#####
molten_p_value_table %>%  
  filter(variable == "depth_p_adjusted") -> p_values_adjusted_depth

fold_change_depth_table_combined %>% 
  merge(., p_values_adjusted_depth, by = c("type", "tech", "filter_level")) %>% 
  mutate(fontface = ifelse(value <= 0.05, 2, 1)) %>% 
  ggplot(aes(y = tech ,x = type)) + 
  geom_tile(aes(fill = fold_change_depth_per_bp)) +
  scale_fill_gradient2(  high = "#1A85FF", low =  "#D41159" , midpoint = 1) + 
  geom_text(aes(y = tech,fontface = fontface,  x = type,label= sprintf("%0.2f", round(fold_change_depth_per_bp, digits = 2))))  + 
  scale_y_discrete(limits = c("ONT", "HiFi", "Illumina"),
                   labels = c("ONT", "HiFi", "Illumina")) +
  facet_grid(.~filter_level) + 
  theme(axis.title = element_blank() ,
        axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1),
        legend.position = "none")

ggsave(paste0("depth_fold_changes_heatmap_filter_comparison.", format(Sys.time(), "%Y-%m-%d"), ".png" ), height= 100, width = 185,units = "mm")
#####

#################################################
# Base quality

# Boxplots

## Illumina 
#####

Illumina_curated %>%
  ggplot(aes(y = basequality, x = type, fill = source)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  stat_summary(fun.y=mean, geom="point", shape=20, position = position_dodge(width = 0.9), aes(group = source, color = "Per motif rate"), size = 5) + 
  theme_classic(base_size = 9) +
  ylab("Mean base quality") + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(labels = scales::comma) +
  ggtitle("Illumina") +
  scale_colour_manual(values = c( "Per motif rate" = "#1a1a1a"  )) +
  guides(colour = guide_legend(override.aes = list(shape = c( 20))),
         fill = guide_legend(override.aes = list(linetype = c( 0),
                                                 shape = c( NA)))) + 
  theme( legend.position = "bottom",
         plot.title = element_text(hjust = 0.5),
         legend.title = element_blank(),
         axis.title.x = element_blank(),
         axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1)) -> Illumina_basequality_boxplots

#####

## HiFi
#####
HiFi_curated %>%
  ggplot(aes(y = basequality, x = type, fill = source)) + 
  stat_summary(fun.data = f, geom = "boxplot", position = "dodge2") +
  stat_summary(fun.y=mean, geom="point", shape=20, position = position_dodge(width = 0.9), aes(group = source, color = "Per motif rate"), size = 5) + 
  theme_classic(base_size = 9) +
  ylab("Mean base quality") + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(labels = scales::comma) +
  ggtitle("HiFi") +
  scale_colour_manual(values = c( "Per motif rate" = "#1a1a1a"  )) +
  guides(colour = guide_legend(override.aes = list(shape = c( 20))),
         fill = guide_legend(override.aes = list(linetype = c( 0),
                                                 shape = c( NA)))) + 
  theme( legend.position = "bottom",
         plot.title = element_text(hjust = 0.5),
         legend.title = element_blank(),
         axis.title.x = element_blank(),
         axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1)) -> HiFi_basequality_boxplots
#####

# Heatmap

fold_change_basequality_table_combined %>%
  ggplot(aes(y = tech ,x = type)) + 
  geom_tile(aes(fill = fold_change_basequality)) +
  scale_fill_gradient2( high = "#1A85FF", low =  "#D41159" , midpoint = 1) + 
  geom_text(aes(y = tech, x = type,label= sprintf("%0.2f", round(fold_change_basequality, digits = 2))), fontface = "bold")  + 
  theme(axis.title = element_blank() ,
        axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1),
        legend.position = "none")

ggsave(paste0("basequality_fold_changes_heatmap_filter_comparison.", format(Sys.time(), "%Y-%m-%d"), ".png" ), height= 80, width = 185,units = "mm")


#################################################
# Probabilistic function to estimate False-positive SNVs
#################################################

# This is the function that is at the core of this estimation:

variantsDueToErrors <- function(rate, genomes, depth, min_reads, min_variants, n_sites){
  # Compute the probability of observing a variant in a site due to sequencing error
  # rate: sequencing error rate per read base
  # genomes: number of individual haploid genomes
  # depth: average sequencing depth
  # min_reads: minimum number of reads supporting the variant
  # min_variants: required minimum number of genomes having the variant in a site
  # n_sites: number of sites corresponding to that error rate
  
  # Sequencing error in a nucleotide: X ??? B(1,rate)
  # Assume independence among errors of different sequenced nucleotides corresponding to the same site
  # Number of errors on a site: Y ??? B(depth,rate)
  # Probability of calling a variant in an haplotype because of sequencing errors: var = P(Y???min_reads) 
  
  var = sum(dbinom(min_reads:depth, size = depth, prob = rate))
  
  # Assume independence among sequencing errors of different haploid genomes
  # Number of haplotypes with a variant in a site due to sequencing errors: V ??? B(genomes, var)
  # Probability of calling a variant in a site due to sequencing errors: P(V???min_variants)
  prob_variant = sum(dbinom(min_variants:genomes, size = genomes, prob = var))
  
  # Expected number of variants in the considered genomic region
  exp_n_variants = prob_variant * n_sites
  
  return(list(prob_variant = prob_variant, exp_n_variants = exp_n_variants))
}


# Approach without simulations

#####

# Re-calculation of error rates with overlapping motifs considered:
#####
Illumina_curated %>%
  filter(source == "motif") %>% 
  group_by(type, overlap_with_other_nonB ) %>%
  summarize(raw_smm_rate = sum(smm)/sum(depth),
            raw_ins_rate = sum(ins)/sum(depth),
            raw_del_rate = sum(del)/sum(depth),
            per_motif_smm_rate = mean(smm_rate, na.rm=T),
            per_motif_ins_rate = mean(ins_rate, na.rm = T), 
            per_motif_del_rate = mean(del_rate, na.rm = T),
            total_length = sum(length)
  ) %>%
  data.frame() %>%
  mutate(tech = "Illumina",
         filter_level = "moderate") -> raw_error_rates_Illumina_moderate

HiFi_curated %>%
  filter(source == "motif") %>% 
  group_by(type, overlap_with_other_nonB ) %>%
  summarize(raw_smm_rate = sum(smm)/sum(depth),
            raw_ins_rate = sum(ins)/sum(depth),
            raw_del_rate = sum(del)/sum(depth),
            per_motif_smm_rate = mean(smm_rate, na.rm=T),
            per_motif_ins_rate = mean(ins_rate, na.rm = T), 
            per_motif_del_rate = mean(del_rate, na.rm = T),
            total_length = sum(length)
  ) %>%
  data.frame() %>%
  mutate(tech = "HiFi",
         filter_level = "moderate") -> raw_error_rates_HiFi_moderate

ONT_curated %>%
  filter(source == "motif") %>% 
  group_by(type, overlap_with_other_nonB ) %>%
  summarize(raw_smm_rate = sum(smm)/sum(depth),
            raw_ins_rate = sum(ins)/sum(depth),
            raw_del_rate = sum(del)/sum(depth),
            per_motif_smm_rate = mean(smm_rate, na.rm=T),
            per_motif_ins_rate = mean(ins_rate, na.rm = T), 
            per_motif_del_rate = mean(del_rate, na.rm = T),
            total_length = sum(length)
  ) %>%
  data.frame() %>%
  mutate(tech = "ONT",
         filter_level = "moderate") -> raw_error_rates_ONT_moderate

rbind.data.frame(raw_error_rates_Illumina_moderate,
                 raw_error_rates_HiFi_moderate,
                 raw_error_rates_ONT_moderate) -> error_rates_filter_level_combined
#####

# Input preparation

#rbind.data.frame(raw_error_rates_Illumina_moderate,
#                 raw_error_rates_HiFi_moderate,
#                 raw_error_rates_ONT_moderate)  %>% 
#  write.table(., file = "error_rates_filter_level_combined.txt", sep = "\t", col.names = T, row.names = F, quote=F)

error_rates_filter_level_combined <- read.table("error_rates_filter_level_combined.txt", header = T)

# add 'B-DNA' values:

# Total length non-B DNA:
error_rates_filter_level_combined %>%  
  filter(  tech == "Illumina" & filter_level == "moderate" & overlap_with_other_nonB == 0) %>% 
  select(total_length) %>% 
  sum()
# 137051792

# Average error rates control for Technologies moderate filter level:
error_rates_filter_level_combined %>% 
  filter(source == "control" & filter_level == "moderate") %>% 
  group_by(tech) %>% 
  summarize(mean_control_smm_rate = mean(raw_smm_rate, na.rm= T))
# HiFi                   0.00181
# Illumina               0.00219
# ONT                    0.0175 

type <- c("BDNA", "BDNA", "BDNA")
source <- c("BDNA", "BDNA", "BDNA")
raw_smm_rate <- c(0.00181, 0.00219, 0.0175 )
total_length <- c(109619844, 109619844, 109619844)
tech <- c("HiFi", "Illumina", "ONT")
filter_level<-c("moderate", "moderate", "moderate")

cbind.data.frame(type, source, raw_smm_rate, total_length, tech, filter_level) %>% 
  mutate(overlap_with_other_nonB = 0) -> temp
error_rates_filter_level_combined %>% 
  mutate(source = "NA") %>% 
  select(type, source, overlap_with_other_nonB, raw_smm_rate, total_length, tech, filter_level) %>% 
  rbind.data.frame(., temp) %>% 
  merge(., seq(from = 3, to = 30, by =3)) %>% 
  rename(depth = y) %>% 
  mutate(min_reads = ceiling(depth/3)) -> probabilistic_calculations_input

probabilistic_calculations_input %>% 
  write.table(., file = "probabilistic_calculations_input.txt", sep = "\t", quote = FALSE, row.names = FALSE)

probabilistic_calculations_input <- read.table("probabilistic_calculations_input.txt", header = T)

# prepare input

probabilistic_calculations_input %>% 
  filter(overlap_with_other_nonB == 0) %>%  
  merge(., c(200, 1000, 10000)) %>%  
  rename(haploid_genome_number = y) %>%  
  merge(., c(1, 3, 0.01)) %>% 
  rename(MAF = y) %>%  
  mutate(input_min_variants = ifelse(MAF < 1, ceiling(MAF*haploid_genome_number), MAF)) %>%  
  select(-source, -overlap_with_other_nonB,  -filter_level) -> input_simplified

variants = mapply(variantsDueToErrors, 
                  input_simplified$raw_smm_rate, 
                  input_simplified$total_length,
                  genomes = input_simplified$haploid_genome_number, 
                  depth = input_simplified$depth, 
                  min_reads = input_simplified$min_reads,
                  min_variants = input_simplified$input_min_variants, SIMPLIFY = FALSE)
variants = matrix(unlist(variants), ncol = 2, byrow = TRUE)
colnames(variants) = c("prob_variant", "exp_n_variants")

# Comparison across technologies

data.frame(variants) %>%
  cbind.data.frame(input_simplified, .) %>% 
  filter( type =="BDNA" ) %>%
  group_by(tech) %>% 
  mutate(exp_n_variants_log = exp_n_variants + 1) -> BDNA_data 

data.frame(variants) %>%
  cbind.data.frame(input_simplified, .) %>% 
  group_by(depth, tech, MAF, haploid_genome_number) %>%
  mutate(total_n_variants = sum(exp_n_variants)) %>%  
  mutate(log_variants = exp_n_variants+1,
         total_n_variants_log = total_n_variants+1) -> nonB_data

nonB_data  %>% 
  ggplot() + 
  geom_line(aes(x = depth, y = total_n_variants_log, color =factor(MAF))) +
  geom_line(data = BDNA_data, aes(x = depth, y = exp_n_variants_log, color = factor(MAF)), linetype ="dashed") + 
  scale_y_continuous(trans = "log10", labels = point) + 
  scale_color_manual(values = c("1" = "#33a02c", "3" =  "#FFC20A", "0.01" =  "#40B0A6"),
                     labels = c(  "1" = "Singletons",   "3" = "Tripletons" ,"0.01" = "MVF = 0.01"),
                     name = "") +
  ylab("Expected false-positive SNVs") + 
  xlim(3,30) + 
  xlab("Read depth") + 
  theme_light(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank()) +
  facet_grid(factor(tech, levels = c("Illumina", "HiFi", "ONT"))~haploid_genome_number,
             labeller = as_labeller(c("200" = "100 Individuals", "1000" = "1,000 individuals", "10000" = "10,000 individuals",
                                      "Illumina" = "Illumina", "HiFi" = "HiFi" , "ONT" = "ONT"
             )))

ggsave(paste0("FP_SNVs_across_tech.", format(Sys.time(), "%Y-%m-%d"), ".png" ), height= 120, width = 185,units = "mm")

# Supplementary Figures: Non-B types as line plots

nonB_data  %>%
  filter(MAF == 1 & type !="BDNA") %>% 
  ggplot() + 
  geom_line(aes(x = depth, y = log_variants, color =type)) +
  scale_y_continuous(trans = "log10", labels = point) + 
  ylab("Expected false-positive SNVs") + 
  xlim(3,30) + 
  xlab("Read depth") + 
  theme_light(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        axis.title.x = element_blank(),
        legend.title = element_blank()) +
  ggtitle("Singletons") + 
  facet_grid(factor(tech, levels = c("Illumina", "HiFi", "ONT"))~haploid_genome_number,
             labeller = as_labeller(c("200" = "100 Individuals", "1000" = "1,000 individuals", "10000" = "10,000 individuals",
                                      "Illumina" = "Illumina", "HiFi" = "HiFi" , "ONT" = "ONT"
             )))

ggsave(paste0("FP_SNVs_across_tech_nonB_types_separately_singletons.", format(Sys.time(), "%Y-%m-%d"), ".png" ), height= 120, width = 185,units = "mm")

nonB_data  %>%
  filter(MAF == 3 & type !="BDNA") %>% 
  ggplot() + 
  geom_line(aes(x = depth, y = log_variants, color =type)) +
  scale_y_continuous(trans = "log10", labels = point) + 
  ylab("Expected false-positive SNVs") + 
  xlim(3,30) + 
  xlab("Read depth") + 
  theme_light(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        axis.title.x = element_blank(),
        legend.title = element_blank()) +
  ggtitle("Tripletons") + 
  facet_grid(factor(tech, levels = c("Illumina", "HiFi", "ONT"))~haploid_genome_number,
             labeller = as_labeller(c("200" = "100 Individuals", "1000" = "1,000 individuals", "10000" = "10,000 individuals",
                                      "Illumina" = "Illumina", "HiFi" = "HiFi" , "ONT" = "ONT"
             )))

ggsave(paste0("FP_SNVs_across_tech_nonB_types_separately_tripletons.", format(Sys.time(), "%Y-%m-%d"), ".png" ), height= 120, width = 185,units = "mm")

nonB_data  %>%
  filter(MAF == 0.01 & type !="BDNA") %>% 
  ggplot() + 
  geom_line(aes(x = depth, y = log_variants, color =type)) +
  scale_y_continuous(trans = "log10", labels = point) + 
  ylab("Expected false-positive SNVs") + 
  xlim(3,30) + 
  xlab("Read depth") + 
  theme_light(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        axis.title.x = element_blank(),
        legend.title = element_blank()) +
  ggtitle("MVF = 0.01") + 
  facet_grid(factor(tech, levels = c("Illumina", "HiFi", "ONT"))~haploid_genome_number,
             labeller = as_labeller(c("200" = "100 Individuals", "1000" = "1,000 individuals", "10000" = "10,000 individuals",
                                      "Illumina" = "Illumina", "HiFi" = "HiFi" , "ONT" = "ONT"
             )))

ggsave(paste0("FP_SNVs_across_tech_nonB_types_separately_MAF_0.01.", format(Sys.time(), "%Y-%m-%d"), ".png" ), height= 120, width = 185,units = "mm")


# Supplementary Figure: Doubletons


probabilistic_calculations_input %>% 
  filter(overlap_with_other_nonB == 0) %>%  
  merge(., c( 1000, 10000)) %>%  
  rename(haploid_genome_number = y) %>%  
  merge(., c(2)) %>% 
  rename(MAF = y) %>%  
  mutate(input_min_variants = ifelse(MAF < 1, ceiling(MAF*haploid_genome_number), MAF)) %>%  
  select(-source, -overlap_with_other_nonB,  -filter_level) -> input_simplified

variants = mapply(variantsDueToErrors, 
                  input_simplified$raw_smm_rate, 
                  input_simplified$total_length,
                  genomes = input_simplified$haploid_genome_number, 
                  depth = input_simplified$depth, 
                  min_reads = input_simplified$min_reads,
                  min_variants = input_simplified$input_min_variants, SIMPLIFY = FALSE)
variants = matrix(unlist(variants), ncol = 2, byrow = TRUE)
colnames(variants) = c("prob_variant", "exp_n_variants")

# Comparison across technologies

data.frame(variants) %>%
  cbind.data.frame(input_simplified, .) %>% 
  filter( type =="BDNA" ) %>%
  group_by(tech) %>% 
  mutate(exp_n_variants_log = exp_n_variants + 1) -> BDNA_data 

data.frame(variants) %>%
  cbind.data.frame(input_simplified, .) %>% 
  group_by(depth, tech, MAF, haploid_genome_number) %>%
  mutate(total_n_variants = sum(exp_n_variants)) %>%  
  mutate(log_variants = exp_n_variants+1,
         total_n_variants_log = total_n_variants+1) -> nonB_data

nonB_data  %>% 
  ggplot() + 
  geom_line(aes(x = depth, y = total_n_variants_log, color =factor(MAF))) +
  geom_line(data = BDNA_data, aes(x = depth, y = exp_n_variants_log, color = factor(MAF)), linetype ="dashed") + 
  scale_y_continuous(trans = "log10", labels = point) + 
  scale_color_manual(values = c("1" = "#33a02c", "3" =  "#FFC20A", "0.01" =  "#40B0A6"),
                     labels = c(  "1" = "Singletons",   "3" = "Tripletons" ,"0.01" = "MVF = 0.01"),
                     name = "") +
  ylab("Expected false-positive SNVs") + 
  xlim(3,30) + 
  xlab("Read depth") + 
  theme_light(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank()) +
  facet_grid(factor(tech, levels = c("Illumina", "HiFi", "ONT"))~haploid_genome_number,
             labeller = as_labeller(c("200" = "100 Individuals", "1000" = "1,000 individuals", "10000" = "10,000 individuals",
                                      "Illumina" = "Illumina", "HiFi" = "HiFi" , "ONT" = "ONT"
             )))

ggsave(paste0("FP_SNVs_DOUBLETONS_across_tech.", format(Sys.time(), "%Y-%m-%d"), ".png" ), height= 120, width = 185,units = "mm")


# Middle G in G-triplets

# Number of middle Gs in G4motifs: 1715082, error rate of middle Gs: 0.035

depth <- seq(from = 3, to = 30, by = 3)
length <- rep(1715082, 10)
rate <- rep(0.035, 10)
exp_FP_singleton <- rep(NA, 10)
exp_FP_tripleton <- rep(NA, 10)
exp_FP_MAF_0.01 <- rep(NA, 10)

cbind.data.frame(depth, length, rate,
                 exp_FP_singleton, 
                 exp_FP_tripleton, 
                 exp_FP_MAF_0.01
) %>% 
  mutate(min_reads = depth/3) -> results

n_genomes <- 200
for (i in 1:10){
  results[i, 4] <- variantsDueToErrors(results[i, 3], n_genomes, results[i,1], results[i, 7], 1, results[i, 2])[[2]] 
  results[i, 5] <- variantsDueToErrors(results[i, 3], n_genomes, results[i,1], results[i, 7], 3, results[i, 2])[[2]] 
  results[i, 6] <- variantsDueToErrors(results[i, 3], n_genomes, results[i,1], results[i, 7], ceiling(0.01*200), results[i, 2])[[2]] 
}

results %>%  
  mutate(n_haploid_genomes = 200) -> temp

n_genomes <- 2000
for (i in 1:10){
  results[i, 4] <- variantsDueToErrors(results[i, 3], n_genomes, results[i,1], results[i, 7], 1, results[i, 2])[[2]] 
  results[i, 5] <- variantsDueToErrors(results[i, 3], n_genomes, results[i,1], results[i, 7], 3, results[i, 2])[[2]] 
  results[i, 6] <- variantsDueToErrors(results[i, 3], n_genomes, results[i,1], results[i, 7], ceiling(0.01*2000), results[i, 2])[[2]] 
}

results %>%  
  mutate(n_haploid_genomes = 2000) %>% 
  rbind.data.frame(., temp)-> temp

n_genomes <- 20000
for (i in 1:10){
  results[i, 4] <- variantsDueToErrors(results[i, 3], n_genomes, results[i,1], results[i, 7], 1, results[i, 2])[[2]] 
  results[i, 5] <- variantsDueToErrors(results[i, 3], n_genomes, results[i,1], results[i, 7], 3, results[i, 2])[[2]] 
  results[i, 6] <- variantsDueToErrors(results[i, 3], n_genomes, results[i,1], results[i, 7], ceiling(0.01*20000), results[i, 2])[[2]] 
}

results %>%  
  mutate(n_haploid_genomes = 20000) %>% 
  rbind.data.frame(., temp)-> temp

temp %>%  
  select(depth, min_reads, n_haploid_genomes, exp_FP_singleton,  exp_FP_tripleton, exp_FP_MAF_0.01, ) %>%  
  melt(., id.vars = c("depth", "min_reads", "n_haploid_genomes")) %>%
  mutate(log_value = value + 1) %>% 
  ggplot() + 
  geom_line(aes(x = depth, y = log_value, color = variable), size = 1) +
  facet_grid(.~n_haploid_genomes, 
             labeller = as_labeller(c("200" = "100 Individuals", "2000" = "1,000 Individuals", "20000" = "10,000 Individuals"))) + 
  xlab("Read Depth") + 
  ylab("Expected number of false-positive SNVs") +
  theme_light(base_size = 10) + 
  scale_y_continuous(trans = "log10", label = point) + 
  scale_color_manual(values = c("exp_FP_singleton" = "#33a02c", "exp_FP_tripleton" = "#FFC20A", "exp_FP_MAF_0.01" = "#40B0A6"),
                     labels = c(  "exp_FP_singleton" = "Singletons",   "exp_FP_tripleton" = "Tripletons" ,"exp_FP_MAF_0.01" = "MVF = 0.01"),
                     name = "") +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) 

ggsave(paste0("FP_GGGs_in_G4s_combined.", format(Sys.time(), "%Y-%m-%d"), ".png" ), height= 80, width = 185,units = "mm")


# Supplementary Figure Doubletons in middle Gs


depth <- seq(from = 3, to = 30, by = 3)
length <- rep(1715082, 10)
rate <- rep(0.035, 10)
exp_FP_doubleton <- rep(NA, 10)

cbind.data.frame(depth, length, rate,
                 exp_FP_doubleton ) %>% 
  mutate(min_reads = depth/3) -> results


n_genomes <- 2000
for (i in 1:10){
  results[i, 4] <- variantsDueToErrors(results[i, 3], n_genomes, results[i,1], results[i, 5], 2, results[i, 2])[[2]] 
}

results %>%  
  mutate(n_haploid_genomes = 2000) -> temp

n_genomes <- 20000
for (i in 1:10){
  results[i, 4] <- variantsDueToErrors(results[i, 3], n_genomes, results[i,1], results[i, 5], 2, results[i, 2])[[2]] 
}

results %>%  
  mutate(n_haploid_genomes = 20000) %>% 
  rbind.data.frame(., temp)-> temp

temp %>%  
  select(depth, min_reads, n_haploid_genomes, exp_FP_doubleton ) %>%  
  melt(., id.vars = c("depth", "min_reads", "n_haploid_genomes")) %>%
  mutate(log_value = value + 1) %>% 
  ggplot() + 
  geom_line(aes(x = depth, y = log_value, color = variable), size = 1) +
  facet_grid(.~n_haploid_genomes, 
             labeller = as_labeller(c( "2000" = "1,000 Individuals", "20000" = "10,000 Individuals"))) + 
  xlab("Read Depth") + 
  ylab("Expected number of false-positive SNVs") +
  theme_light(base_size = 10) + 
  scale_y_continuous(trans = "log10", label = point) + 
  scale_color_manual(values = c("exp_FP_doubleton" = "#33a02c"),
                     labels = c(  "exp_FP_doubleton" = "doubletons"),
                     name = "") +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) 

ggsave(paste0("FP_GGGs_DOUBLETONS_in_G4s_combined.", format(Sys.time(), "%Y-%m-%d"), ".png" ), height= 80, width = 185,units = "mm")


# Supplementary Table 1000G, SGDP, and gnomAD

GGGs_df <- data.frame(list('GGGs', 0.035, 1715082))
names(GGGs_df) <- c("type", "raw_smm_rate", "total_length")

probabilistic_calculations_input %>% 
  filter(overlap_with_other_nonB == 0 & depth == 3 & tech == "Illumina") %>% 
  select(type, raw_smm_rate, total_length) %>% 
  rbind.data.frame(., GGGs_df) %>% 
  merge(., list("haploid_genome_number" = c(5008, 600, 152312), "depth" = c(1,21,15)))  %>% 
  merge(., list("MAF" = c(1, 3, 0.01))) %>% 
  mutate(input_min_variants = ifelse(MAF < 1, ceiling(MAF*haploid_genome_number), MAF),
         min_reads = ceiling(depth / 3)) -> input_simplified

variants = mapply(variantsDueToErrors, 
                  input_simplified$raw_smm_rate, 
                  input_simplified$total_length,
                  genomes = input_simplified$haploid_genome_number, 
                  depth = input_simplified$depth, 
                  min_reads = input_simplified$min_reads,
                  min_variants = input_simplified$input_min_variants, SIMPLIFY = FALSE)
variants = matrix(unlist(variants), ncol = 2, byrow = TRUE)
colnames(variants) = c("prob_variant", "exp_n_variants")


data.frame(variants) %>%
  cbind.data.frame(input_simplified, .) %>% 
  filter( type =="BDNA" ) %>%
  group_by(haploid_genome_number, MAF) -> BDNA_data 

data.frame(variants) %>%
  cbind.data.frame(input_simplified, .) %>% 
  filter(type %not_in% c("BDNA", "GGGs")) %>% 
  group_by(haploid_genome_number, MAF) %>%
  mutate(total_n_variants = sum(exp_n_variants)) %>%  
  select(-exp_n_variants) %>% 
  rename("exp_n_variants" = "total_n_variants") %>% 
  filter(type == "APhasedRepeats") -> nonB_data

data.frame(variants) %>%
  cbind.data.frame(input_simplified, .) %>% 
  filter(type == "GGGs") -> GGGs_data

rbind.data.frame(BDNA_data, nonB_data, GGGs_data) %>% 
  select(type, haploid_genome_number, MAF, exp_n_variants) %>% 
  dcast(haploid_genome_number+MAF~type, value.var = "exp_n_variants")  %>% 
  arrange(match(haploid_genome_number, c("5008", "600", "152312"))) %>% 
  rename("Dataset" = "haploid_genome_number", "NonB" = "APhasedRepeats") %>% 
  arrange(match(MAF, c(1,3,0.01))) %>% 
  write.table(., file = "Supplementary_Table_FP_SNVs_1000G_SGDP_gnomAD.txt", row.names = FALSE, sep = "\t", quote = FALSE)

# Supplementary Figure simulations of probabilistic approach

sim_dat <- read.table("output_simulations_example", header = F)
names(sim_dat) <- c("type", "raw_smm_error_rate", "lenght", "haploid_genome_number", "depth", "singletons", "MAF_0.01", "MAF_0.05")

sim_dat %>% 
  summary()

sim_dat %>%
  filter( type =="BDNA" ) %>% 
  mutate(log_singletons = singletons + 1, 
         group_plot = "BDNA") %>%
  select(group_plot, depth, log_singletons) -> BDNA_data

sim_dat %>% 
  filter(type %not_in% c("BDNA", "GGGs")) %>% 
  mutate(simulation_number = rep(c(1:100), 60)) %>%  
  group_by(depth, simulation_number) %>%
  mutate(total_n_variants = sum(singletons)) %>% 
  filter(type == "APhasedRepeats" ) %>% 
  mutate(log_singletons = total_n_variants + 1,
         group_plot = "nonB") %>%  
  data.frame() %>% 
  select(group_plot, depth, log_singletons) -> nonB_data
  
rbind.data.frame(BDNA_data, nonB_data) %>% 
ggplot() +
  geom_boxplot(aes(x = factor(depth), y = log_singletons,  color = group_plot) , position = "dodge2") +
  scale_y_continuous(trans = "log10", label = point) +
  ylab("Expected number of false-positive SNVs") +
  xlab("Read depth") + 
  theme_light(base_size = 10) + 
  ggtitle("Singletons - 100 Individuals") + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank())

ggsave(paste0("FP_singletons_simulations_200_inds.", format(Sys.time(), "%Y-%m-%d"), ".png" ), height= 100, width = 185,units = "mm")

# Supplementary Figure False-positive SNVs from HG002 Illumina data:

dat <- read.table("FPs_input_data", header = F)

head(dat)
names(dat) <- c("stats", "values", "source", "depth", "total_length", "callset")
dat %>% 
  dcast( source+depth+total_length+callset~stats, value.var="values") -> dat_cast

names(dat_cast) <- c("source", "depth", "total_length", "callset", "fn", "fp", "total_n", "tp")

dat_cast %>%
  filter( depth < 100) %>% 
  mutate(fraction_fp = fp/((tp+fp)/100)) %>% 
  ggplot(aes(x = factor(depth), y = fraction_fp, fill = factor(source, levels=c("motif", "control")))) +
  geom_bar(stat = "identity", position = "dodge2") +
  scale_fill_manual(values = colors, name = "") +
  xlab("Depth") +
  theme_light(base_size = 10) +
  ylab("% false-positives") +
  theme(legend.position = "bottom") +
  facet_wrap(.~callset)

ggsave(paste0("FPs_giab_v3.3.2_v4.2.1.", format(Sys.time(), "%Y-%m-%d"), ".png" ), height= 100, width = 185,units = "mm")

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

Intersected with the G4 motifs from Illumina_curated_fully_filterd (which already has the mappability  = 1 filter) this yields  1,715,082 potential middle Gs in G4 motifs. 


### Empirical false-positive SNV estimation

In addition, we estimated the proportion of false-positive SNVs using the same Illumina data set that was also employed for the error detection analyses. To illustrate the effect of varying read depths on false-positives, we sub-sampled the HG002 100X data set to approximately 3, 9, 15, and 30��, and performed variant calling using the tool FreeBayes with default parameters and a minimum mapping quality filter set to 30.

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
