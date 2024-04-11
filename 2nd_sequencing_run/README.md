# Directory Structure

| directory-name | date-created | dataset |
| -------------- | ------------ | ------- |
| 2nd_sequencing_run | 2023-12-05 | Alb & Contemp pair from Verde Island|

# Pseudanthias squamipinnis lcWGS

---

Jem Baldisimo

---
This repository outlines the roadmap we followed to move *Pseudanthias squamipinnis* through the [Low Coverage Whole Genome Sequencing Pipeline](https://github.com/philippinespire/pire_lcwgs_data_processing). This provides steps taken & analysis we did to gain insight on the genetic diversity and historical population demography of this very popular aquarium fish, also known as the sea goldie.

---

## 1. Completed fq.gz pre-processing

The [pire_fq_gz_processing](https://github.com/philippinespire/pire_fq_gz_processing) instructions and scripts (Garcia et al 2021) were used. The purpose of the preprocessing pipeline was to Trim, deduplicate, decontaminate, and repair the raw `fq.gz` files.

<details><summary><i>Checked quality of data using Multi_FASTQC</i></summary>
<p>

Multi_FASTQC [Report](https://github.com/philippinespire/pire_pseudanthias_squamipinnis_lcwgs/blob/main/2nd_sequencing_run/fq_raw/fqc_raw_report.html)
```
Potential issues:
 * duplication - low to high in Albatross; low to high in Contemporary but only due to 1 outlier. The rest of contemp are moderate (10 to 20%)
   * Alb: 5-50.9%, Contemp: 2.9-51.2-%
 * gc content - reasonable
   *39-59%
 *number of reads
   *Albatross - 0.1-123.5 M seq; Contemporary - 0-72.1 M seq
```
</p>
</details>

<details><summary><i>1st trim</i></i></summary>
<p>

 1st FASTP [Report]()
```
Potential issues:
 * duplication - Low to moderate in Albatross, low in Contemp
   * Alb: 1.3 to 22.8%, Contemp: 1.3-15.8%
 * gc content - Reasonable
   * Alb: 34.2-38.8%, Contemp: 37.3-43.9%
 * Passing filter
   * Alb: 92.9-97.2, Contemp: 56.5-97.3%
 * % adapter 
   * Alb: 78.5-94.4%, Contemp: 15.5-63.9%
 * number of reads - Range of reads very erratic for both
   * Albatross - 0.1 - 108 M seq; Contemp - 16,000 to 60 M seq
```
</p>
</details>

 Clumpify [Report]()

```
bash ../../pire_fq_gz_processing/runCLUMPIFY_r1r2_array.bash fq_fp1 fq_fp1_clmp /scratch/jbald004 20
#copied the checkClumpify_EG.R file to the fq_fp1_clmp folder after running clumpify
cp /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_equula/fq_fp1_clmp/checkClumpify_EG.R .
#ran checkClumpify_EG.R in the fq_fp1_clmp folder
cd /home/e1garcia/shotgun_PIRE/pire_lcwgs_data_processing/sphaeramia_nematoptera/fq_fp1_clmp
salloc
enable_lmod
module load container_env mapdamage2
crun R <checkClumpify_EG.R --no-save
```
Clumpify worked well so I moved on to the next step
</p>
</details>

<details><summary><i>2nd Trim</i></summary>
<p>

```
sbatch ../../pire_fq_gz_processing/runFASTP_2_cssl.sbatch fq_fp1_clmp fq_fp1_clmp_fp2
```

Fastp2 [Report]()
```
Potential issues:
  * % duplication - low, with 1-3 outliers
    * Alb: 0.7-34.1% (only 1 individual, the rest are <2%, 1 is 12.1%), Contemp: 0.3-16.5% (16.5 is an outlier)
  * GC content - acceptable range of ~40%
    * Alb: 37.5-52.2%, Contemp: 45.3-49%
  * passing filter - very high
    * Alb: 98.3-99.7%, Contemp: 99.0-99.6%
  * % adapter - low
    * Alb: 0.9-1.8%, Contemp: 0.2-1.7%
  * number of reads - very erratic, no patttern across samples
    * Alb: 1-148 M, Contemp: 1-270 M
```
</p>
</details>

<details><summary><i>Checked Fastqscreen files</i></summary>
<p>

Fastqscreen [Report]()

```
Potential issues 
  * one hit, one genome, no ID -
    * Alb: 78-97.5%, Contemp: 93-96%
  * No one hit, one genome to any potential contaminators (bacteria, virus, human, etc) -
    * Alb: 2.5-22%, Contemp: 4-7%
```
</p>
</details>

<details><summary><i>Revewied Fastqscreen repaired files</i></summary>
<p>

MultiQC [Report]()

```
Potential issues:
 * % duplication - 
   Alb: 0.7-21.5% (only one had 21.5, the rest are only up to 7%), Contemp: 0.7-14.5%
 * GC content -
   Alb: 37-48%, Contemp: 45-48%
 * number of reads -
   Alb: - 0.5-71.7 M (wide range), Contemp: 0.3-131.4 M
```
