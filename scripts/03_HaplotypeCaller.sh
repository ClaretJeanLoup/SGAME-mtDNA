#!/bin/bash
#SBATCH -J HaplotypeCaller
#SBATCH -c 12
#SBATCH --mem=64G
#SBATCH -t 24:00:00
##sbatch 03_HaplotypeCaller.sh -b /path/to/sample.mkdup.bam -g /path/to/ref.fa -o /path/to/gvcf_outdir

while getopts b:g:o: flag
do
case "${flag}" in
                b) BAM=${OPTARG}
                	;;
                g) ref=${OPTARG}
                        ;;
                o) outdir=${OPTARG}
                        ;;
                *) echo "Invalid option: -$flag"
                	;;
        esac
done

for v in BAM ref outdir; do [[ -z "${!v}" ]] && { echo "Missing argument: $v"; exit 1; }; done
[[ -e "$BAM" ]] || { echo "Not found: $BAM"; exit 1; }
[[ -e "$ref" ]] || { echo "Not found: $ref"; exit 1; }

echo "## ------------------------------------------------ ##"
echo `date` ": sync files."
rsync -ta $BAM $BAM.bai $TMPDIR
rsync -ta ${ref} ${ref}.fai ${ref%.*}.dict $TMPDIR
cd $TMPDIR

genome=`basename $ref`
bam=`basename $BAM`
sample=${bam%.mkdup.bam}
echo "Sample $sample, reference $genome"

echo "## ------------------------------------------------ ##"
echo `date` ": HaplotypeCaller."
gatk --java-options "-Xmx48G" HaplotypeCaller -R $genome -I $bam -ERC GVCF \
    --native-pair-hmm-threads 12 --tmp-dir $TMPDIR -O ${sample}.g.vcf.gz

mkdir -p $outdir
rsync -ta ${sample}.g.vcf.gz* $outdir
