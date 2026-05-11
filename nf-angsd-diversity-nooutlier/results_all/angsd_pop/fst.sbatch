#!/bin/bash -l

#SBATCH --job-name=angsd_fst
#SBATCH -o angsd_fst-%j.out 
#SBATCH -p main 
#SBATCH -n 1  
#SBATCH -N 1  
#SBATCH --cpus-per-task=8 

enable_lmod
module load container_env angsd

export SINGULARITY_BIND=/home/e1garcia

# --- CONFIGURATION ---
SAFDIR="/archive/carpenterlab/pire/pire_pseudanthias_squamipinnis_lcwgs/nf-angsd-diversity-nooutlier/results_all/angsd_pop/"
THREADS=8

POPS=(
"PsqAGal_historic"
"PsqCVer_modern"
)

WINDOW=50000
STEP=10000
# ---------------------

cd $SAFDIR

echo "🚀 Starting FULL FST pipeline (SFS → FST)..."

for (( i=0; i<${#POPS[@]}; i++ )); do
  for (( j=i+1; j<${#POPS[@]}; j++ )); do

    P1=${POPS[$i]}
    P2=${POPS[$j]}
    PAIR_NAME="${P1}_vs_${P2}"

    echo "===================================="
    echo "Processing: $PAIR_NAME"
    echo "===================================="

    # Check SAF files
    if [[ ! -f ${P1}.saf.idx || ! -f ${P2}.saf.idx ]]; then
      echo "⚠️ Missing SAF index — skipping"
      continue
    fi

    # -----------------------------
    # 1. Compute 2D SFS (REQUIRED)
    # -----------------------------
    if [[ ! -f ${PAIR_NAME}.2dSFS ]]; then
      echo "➡️ Computing FOLDED 2D SFS..."
      crun.angsd realSFS \
        -fold 1 \
        ${P1}.saf.idx ${P2}.saf.idx \
        -P ${THREADS} \
        > ${PAIR_NAME}.2dSFS
    else
      echo "✅ 2D SFS already exists"
    fi
    # -----------------------------
    # 2. Build FST index
    # -----------------------------
    if [[ ! -f ${PAIR_NAME}.fst.idx ]]; then
      echo "➡️ Building FST index..."
      crun.angsd realSFS fst index \
        ${P1}.saf.idx ${P2}.saf.idx \
        -sfs ${PAIR_NAME}.2dSFS \
        -fstout ${PAIR_NAME}
    else
      echo "✅ FST index exists"
    fi

    # -----------------------------
    # 3. Global FST
    # -----------------------------
    echo "➡️ Computing global FST..."
    crun.angsd realSFS fst stats \
      ${PAIR_NAME}.fst.idx \
      > ${PAIR_NAME}.global.fst.txt

    # -----------------------------
    # 4. Windowed FST
    # -----------------------------
    echo "➡️ Computing windowed FST..."
    crun.angsd realSFS fst stats2 \
      ${PAIR_NAME}.fst.idx \
      -win ${WINDOW} \
      -step ${STEP} \
      > ${PAIR_NAME}.windowed.fst.txt

    echo "✅ Done: $PAIR_NAME"
    echo ""

  done
done

echo "🎉 All FST calculations complete!"
