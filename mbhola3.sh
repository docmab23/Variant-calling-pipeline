#!/bin/bash
realign=0
index=0
gunzip=0                              # some default variables as by default improvement , verbose ,indexing is "off" if user doesn't specify so 
v=0                                      
help=0
output="false"                   # for the vcf file the name should be in "filename.vcf" format
verbose=0
#getopts block showing different flags , some require argument , some just change the default value.
while getopts "a:b:r:eio:zvf:h" arg; 
do   
  case $arg in 
  a)
  reads1=$OPTARG
  ;;
  b)
  reads2=$OPTARG
  ;;    
  r)
  ref=$OPTARG
  ;;
  e)
  realign=1
  ;;
  i)
  index=1
  ;;
  o)
  output=$OPTARG
  ;;
  z)
  gunzip=1
  ;;
  v)
  verbose=1
  ;;
  f)
  millsFile=$OPTARG
  ;;
  h)
  help=1
  ;;
   esac
done
#help flag to give info about the usage of the script with -h flag
if [[ $help -eq 1 ]];
then
   echo "-r Reference file  -a Read file one -b Read file 2 -e Realign and Improvements -i index output bam file
    -o give output name for vcf file(shoudl be in .vcf format) -z Zip the Output vcf file -v verbose -f Mills file path(necessary for realignmnet and improvements) -h help"
    exit

fi
#to change the verbose settings.
if [[ $verbose -eq 1 ]];
then
   set -x
fi
# this segment basically checks for the reference and reads file and exits if not availible.
if [[ ! -f $ref ]]; 
then 
  echo "No refrence file specificed , check help" 
  exit 
fi

if [[ ! -f $reads1  &&  ! -f $reads2 ]];
then
   echo "Both read files missing ,check help"
   exit
elif [[ ! -f $reads1 ]];
then
   echo "reads1 is missing , check help"
exit
elif [[ ! -f $reads2 ]];
then 
  echo "reads2 is missing,check help"
exit
fi
if [[ $output -eq "false" ]];
then 
echo "You need to give an output vcf filename as 'filename.vcf', using -o  "
exit
fi
#this checks for the specified output vcf file and asks the user to overwrite or exit if the zipped version exists.(also applies to the default output name ) 
#and deletes the file if user wants to overwrite as it would be created in the next block.
if [ -f "$output.gz" ];
then 
   echo "The file you specified already exists ,Do you want to overwrite the file or exit , Type overwrite if you want to OVERWRITE it otherwise type EXIT"
   read answer;
   if [[ $answer -eq "EXIT" ]];
   then
   exit
   elif [[ $answer -eq "OVERWRITE" ]];
   then 
    rm -r "$output.gz"  # this is removed as its anyways going to get created in the pipeline now .
   fi
fi
#mills file would be checked in the improvement step as it is required their only.
#burrows wheeler alignment step , this prepares the reference file for mapping 
command bwa index $ref  
# this step maps the reads to the reference file
command bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' $ref $reads1 $reads2 > lane.sam
#conversion of sam to bam files with samtools followed by sorting 
command samtools fixmate -O bam lane.sam lane_fixmate.bam
command samtools sort -O bam -o lane_sorted.bam lane_fixmate.bam
samtools index lane_sorted.bam
#creates dictionary and .fai files.
samtools dict $ref  
samtools faidx $ref
#checks if user wants to perform realignment or not
if [ $realign -eq 1 ];
then
  #checks for mills file.
  if [ ! -f $millsFile ];
  then 
  echo "No mills file specified, so realignment can't happen"
  exit
  fi
  # realignment steps
  command  java -Xmx2g -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R chr17.fa -I lane_sorted.bam -o lane.intervals --known $millsFile 1> mbhola3_1.log 2>&1
  command java -Xmx4g -jar GenomeAnalysisTK.jar -T IndelRealigner -R $ref -I lane_sorted.bam -targetIntervals lane.intervals -o lane_realigned.bam 1 -known $millsFile 1> mbhola3_2.log 2>&1
  #if the user wants to index the realigned bam file
  if [ $index -eq 1 ];
  then  
  samtools index lane_realigned.bam
  fi 
  
  bcftools mpileup -Ob -o inter.bcf -f $ref lane_realigned.bam

  else # for the non-improvement side.
  bcftools mpileup -Ob -o inter.bcf -f $ref lane_sorted.bam
  fi
  
  bcftools call -vmO z -o $output.gz inter.bcf
  gunzip $output.gz
  # conversion of vcf file to bed format and prepartion of sepearte indels and snps files
  sed -n '/^[chr\d]/,$p' $output | awk '{print $1,$2, length($5)-length($4)+$2,length($5)-length($4)}' OFS='\t'> inter.bed
  sed -i 's/chr//g' inter.bed
  awk '{if($4!=0) print;}' inter.bed > tmp && mv tmp indels.txt
  awk '{if($4==0) print;}' inter.bed >snps.txt
  if [ $gunzip=1 ]
  then
    gzip $output
  fi

  









