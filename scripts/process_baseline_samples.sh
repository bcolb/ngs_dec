#!/bin/bash

echo "starting"

#BAM_LIST is a file of absolute paths to each bam file
BAM_LIST=$1

# Replace commented paths below with the files needed to generate your assay specific error profile
# Required input files include: BAM_LIST (command line argument above), intervals file, capture BED file, reference genome
# Final output is the ERROR_BASELINE file, CONSOLIDATED_ERROR_PROFILE is an intermediary file 
VARSCAN=PATH/TO/VARSCAN/VERSION #/mnt/disk3/genetics/venvs/python-2.7.11/tgc-pipe/tgc-pipe-env/bin/VarScan.v2.3.7.jar
READCOUNTS_SCRIPT='ngc_dec.py readcounts_to_error_profile'
INTERVALS_FILE=PATH/TO/INTERVALS/FILE
BEDFILE=PATH/TO/BED/FILE #/mnt/disk2/com/Genomes/BED_Files/OncoPlex_v5_0734011.bed
REF_GENOME=data/refgene/OPX_v5_refGene.txt
CONSOLIDATED_ERROR_PROFILE=/OUTPUT/PATH/FOR/CONSOLIDATED/ERROR/FILE #consolidated_error_profile_08142017_OPXv5.csv
ERROR_BASELINE=/OUTPUT/PATH/FOR/FINAL/ERROR/BASELINE/FILE #error_baseline_08142017_OPXv5.csv

echo "bams"
for BAM in `sed '/^$/d' $BAM_LIST`; do
    SAVEPATH=$(dirname $BAM)
    echo "SAVEPATH: $SAVEPATH"
    BAMNAME=$(basename $BAM)
    echo "BAMNAME: $BAMNAME"
    PFX=${BAMNAME%.*.*}
    echo "PFX: $PFX"
    SAVEPATH='output/error_baseline/'
    mkdir -p $SAVEPATH/$PFX

    echo “Starting Analysis of $PFX” 

    echo "Making mpileups" 
    date +"%D %H:%M" >> $SAVEPATH/$PFX/msi_run_log.txt;
    
    ##### $BAM is the absolute path of the original final.bam 
    samtools mpileup -f $REF_GENOME -d 100000 -A -E $BAM -l $INTERVALS_FILE | awk '{if($4 >= 6) print $0}' > $SAVEPATH/$PFX/$PFX.mpileup2 
    
    echo "Varscan Readcounts start" >> $SAVEPATH/$PFX/msi_run_log.txt;
    cmd="java -Xmx4g -jar $VARSCAN readcounts $SAVEPATH/$PFX/$PFX.mpileup2 --variants-file $INTERVALS_FILE --min-base-qual 15 --output-file $SAVEPATH/$PFX/$PFX.varscan_readcounts &"
    echo $cmd
    $cmd
    wait
    
    echo "Readcounts parsing script" 
    cmd2="python $READCOUNTS_SCRIPT $SAVEPATH/$PFX/$PFX.varscan_readcounts $SAVEPATH/$PFX/$PFX.variant_profile.csv"
    echo $cmd2
    $cmd2

done

# Consolidate individual variant profiles into a single variant profile
cmd="python ngs_dec.py consolidate_error_profiles /home/illumina/hiseq/analysis/170201_NA0174_OncoPlexKAPA135R-OPXv5/output/modified_msings/ $CONSOLIDATED_ERROR_PROFILE"
echo $cmd
$cmd

# Turn consolidated variant profiles into baseline file
cmd="python ngs_dec.py error_profile_to_baseline $CONSOLIDATED_ERROR_PROFILE $ERROR_BASELINE"
echo $cmd
$cmd
