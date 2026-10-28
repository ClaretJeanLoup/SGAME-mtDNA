#!/bin/bash
#SBATCH -J UniqueVariants
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 12:00:00

# Filters and removes background mutations shared between lines and merges the unique variants of callers.
##sbatch 06_UniqueVariants.sh -v /path/to/gatk_joint.vcf.gz -m /path/to/mpileup_dir -g /path/to/ref.fa -o /path/to/outdir [-a min_mtDNA_fraction]

minaf=0.02

while getopts v:m:g:o:a: flag
do
case "${flag}" in
                v) VCF=${OPTARG}
                    	;;
                m) MPDIR=${OPTARG}
                        ;;
                g) ref=${OPTARG}
                        ;;
                o) outdir=${OPTARG}
                        ;;
                a) minaf=${OPTARG}
                        ;;
                *) echo "Invalid option: -$flag"
                      	;;
        esac
done

for v in VCF MPDIR ref outdir; do [[ -z "${!v}" ]] && { echo "Missing argument: $v"; exit 1; }; done
[[ -e "$VCF" ]] || { echo "Not found: $VCF"; exit 1; }
[[ -e "$MPDIR/mpileup.MtDNA.vcf.gz" ]] || { echo "Not found: $MPDIR/mpileup.MtDNA.vcf.gz"; exit 1; }
[[ -e "$ref" ]] || { echo "Not found: $ref"; exit 1; }

helper=$(dirname $(readlink -f $0))/06_unique_variants.py
# Under slurm the script runs from a spool copy, so the helper is looked up next to the submitted script as well.
[[ -e "$helper" ]] || helper=$SLURM_SUBMIT_DIR/06_unique_variants.py
[[ -e "$helper" ]] || { echo "Not found: 06_unique_variants.py"; exit 1; }

mkdir -p $outdir
cd $outdir

echo "## ------------------------------------------------ ##"
echo `date` ": GATK calls."
bcftools norm -f $ref -m -any $VCF -Oz -o gatk.split.vcf.gz
bcftools view -t ^MtDNA gatk.split.vcf.gz -Oz -o gatk.nuclear.vcf.gz
python3 $helper gatk.nuclear.vcf.gz GATK gatk.nuclear
# GATK genotypes are diploid; on MtDNA the variant fraction is read from AD as for BCFtools.
bcftools view -t MtDNA gatk.split.vcf.gz -Oz -o gatk.MtDNA.vcf.gz
python3 $helper gatk.MtDNA.vcf.gz GATK gatk.MtDNA --min-af $minaf

echo "## ------------------------------------------------ ##"
echo `date` ": BCFtools calls."
bcftools concat -a `ls $MPDIR/mpileup.*.vcf.gz | grep -v mpileup.MtDNA` -Oz -o bcftools.nuclear.vcf.gz
python3 $helper bcftools.nuclear.vcf.gz BCFtools bcftools.nuclear
python3 $helper $MPDIR/mpileup.MtDNA.vcf.gz BCFtools bcftools.MtDNA --min-af $minaf

echo "## ------------------------------------------------ ##"
echo `date` ": merging."
tables="gatk.nuclear bcftools.nuclear gatk.MtDNA bcftools.MtDNA"
# One row per variant and line, listing the callers that found it (values are those of the last caller).
for t in $tables; do tail -n +2 $t.tsv; done \
    | awk -F'\t' -v OFS='\t' '{k=$2 FS $3 FS $4 FS $5 FS $6; c[k]=(k in c) ? c[k]","$1 : $1; r[k]=$0}
        END {for (k in r) {$0=r[k]; $1=c[k]; print}}' \
    | sort -t$'\t' -k2,2V -k3,3n > body.tsv
cat <(head -1 gatk.nuclear.tsv) body.tsv > unique_variants.tsv && rm body.tsv

for t in $tables; do
    bcftools sort $t.sites.vcf -Oz -o $t.sites.vcf.gz && bcftools index -t $t.sites.vcf.gz
done
bcftools concat -a -D `for t in $tables; do echo $t.sites.vcf.gz; done` -Oz -o unique_variants.sites.vcf.gz
bcftools index -t unique_variants.sites.vcf.gz
echo "Unique variants (variant x line): `tail -n +2 unique_variants.tsv | wc -l`, of which MtDNA: `awk -F'\t' '$2=="MtDNA"' unique_variants.tsv | wc -l`"
