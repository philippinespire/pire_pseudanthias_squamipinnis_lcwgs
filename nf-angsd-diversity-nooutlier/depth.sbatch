#!/bin/bash -l

#SBATCH --job-name=angsd_depth
#SBATCH -o angsd_dep-%j.out 
#SBATCH -p main 
#SBATCH -n 1  
#SBATCH -N 1  
#SBATCH --cpus-per-task=40 

enable_lmod
module load container_env angsd
export SINGULARITY_BIND=/home/e1garcia

OUTDIR=$1


cd ${OUTDIR}

crun angsd \
-b /archive/carpenterlab/pire/pire_pseudanthias_squamipinnis_lcwgs/nf-angsd-diversity-nooutlier/results/inputfiles/bamlist.txt \
-ref /archive/carpenterlab/pire/pire_pseudanthias_squamipinnis_lcwgs/2nd_sequencing_run/GenErode/reference/reference.denovoSSL.Psq20k.rename.fasta \
-out depth_all_sites \
-doCounts 1 \
-doDepth 1 \
-dumpCounts 1 \
-GL 0 \
-doMajorMinor 0 \
-doMaf 0 \
-minMapQ 25 \
-minQ 30 \
-uniqueOnly 1 \
-remove_bads 1 \
-P 4
