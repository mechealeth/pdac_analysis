#!/bin/bash       
#SBATCH --job-name=SAMPLENAME
#SBATCH --partition=fat-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=100
#SBATCH --error=%j.err    
#SBATCH --output=%j.out    

CURDIR=`pwd`
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID    
NODES=`scontrol show hostnames $SLURM_JOB_NODELIST`    
for i in $NODES    
do    
echo "$i:$SLURM_NTASKS_PER_NODE" >> $CURDIR/nodelist.$SLURM_JOB_ID    
done
echo $SLURM_NPROCS

echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"

echo "start running spaceranger count..."
#setting environment
source /lustre/home/nzhang/software/cellranger/spaceranger-1.2.2/sourceme.bash

spaceranger count --id=SAMPLENAME \
	--transcriptome=/lustre/home/whhou/00.datasets/10x_2021/refdata-gex-GRCh38-2020-A \
	--fastqs=/lustre/home/whhou/00.datasets/2021_Data/process/data/PDAC/ST/SAMPLENAME \
	--sample=SAMPLENAME \
	--image=IMAGEFILE \
	--slide=SLIDE \
	--area=AREA \
	--slidefile=SLIDEFILE

echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID


