#!/bin/bash
#SBATCH -J Mpileup
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 3-00:00:00

# Calls small variants on MtDNA every observed alternative allele is kept with its read counts so low-frequency variants are not lost.
##sbatch 05_Mpileup.sh -d /path/to/bam_dir -g /path/to/ref.fa -r contig -o /path/to/outdir

while getopts d:g:r:o: flag
do
case "${flag}" in
                d) DIRDATA=${OPTARG}
                    	;;
                g) ref=${OPTARG}
                        ;;
                r) region=${OPTARG}
                        ;;
                o) outdir=${OPTARG}
                        ;;
                *) echo "Invalid option: -$flag"
                      	;;
        esac
done

for v in DIRDATA ref region outdir; do [[ -z "${!v}" ]] && { echo "Missing argument: $v"; exit 1; }; done
[[ -d "$DIRDATA" ]] || { echo "Not found: $DIRDATA"; exit 1; }
ls $DIRDATA/*.mkdup.bam >/dev/null 2>&1 || { echo "No *.mkdup.bam in $DIRDATA"; exit 1; }
[[ -e "$ref" ]] || { echo "Not found: $ref"; exit 1; }

mkdir -p $outdir
ls $DIRDATA/*.mkdup.bam > $outdir/bams.${region}.list

echo "## ------------------------------------------------ ##"
echo `date` ": mpileup on $region."
if [[ "$region" == "MtDNA" ]]; then
    bcftools mpileup -f $ref -r $region -b $outdir/bams.${region}.list -d 100000 -q 50 \
        -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,INFO/AD -Ou \
        | bcftools norm -f $ref -m -any -Ou \
        | bcftools view -e 'ALT="<*>"' -Oz -o $outdir/mpileup.${region}.vcf.gz
else
    bcftools mpileup -f $ref -r $region -b $outdir/bams.${region}.list \
        -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,INFO/AD -Ou \
        | bcftools call -mv -Ou \
        | bcftools norm -f $ref -m -any -Oz -o $outdir/mpileup.${region}.vcf.gz
fi
bcftools index -t $outdir/mpileup.${region}.vcf.gz
