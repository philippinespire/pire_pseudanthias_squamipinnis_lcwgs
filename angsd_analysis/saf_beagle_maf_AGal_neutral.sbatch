#!/bin/bash -l

#SBATCH --job-name=angsd_saf
#SBATCH -o angsd_saf-%j.out 
#SBATCH -p main 
#SBATCH -n 1  
#SBATCH -N 1  
#SBATCH --cpus-per-task=8 

enable_lmod
module load container_env angsd
export SINGULARITY_BIND=/home/e1garcia

OUTDIR=$1


cd ${OUTDIR}

crun angsd \
        -b bam_list_all_subset_AGal.txt \
        -anc /archive/carpenterlab/pire/pire_pseudanthias_squamipinnis_lcwgs/2nd_sequencing_run/GenErode/reference/reference.denovoSSL.Psq20k.rename.fasta \
        -ref /archive/carpenterlab/pire/pire_pseudanthias_squamipinnis_lcwgs/2nd_sequencing_run/GenErode/reference/reference.denovoSSL.Psq20k.rename.fasta \
        -out agal_sites_neutral \
      -doSaf 1 -noTrans 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -dumpCounts 1  \
        -P 8 \
	-sites global_snp_list_neutral.txt -rf global_snp_list_depth1_15_notrans.chrs
