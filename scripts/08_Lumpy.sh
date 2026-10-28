#!/bin/bash
#SBATCH -J Lumpy
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 12:00:00

# Calls structural variants of one bam from split reads and discordant pairs with LUMPY (lumpyexpress).
##sbatch 08_Lumpy.sh -b /path/to/sample.mkdup.bam -o /path/to/outdir

while getopts b:o: flag
do
case "${flag}" in
                b) BAM=${OPTARG}
                	;;
                o) outdir=${OPTARG}
                        ;;
                *) echo "Invalid option: -$flag"
                	;;
        esac
done

for v in BAM outdir; do [[ -z "${!v}" ]] && { echo "Missing argument: $v"; exit 1; }; done
[[ -e "$BAM" ]] || { echo "Not found: $BAM"; exit 1; }

echo "## ------------------------------------------------ ##"
echo `date` ": sync files."
rsync -ta $BAM $BAM.bai $TMPDIR
cd $TMPDIR

bam=`basename $BAM`
sample=${bam%.mkdup.bam}
echo "Sample $sample"

echo "## ------------------------------------------------ ##"
echo `date` ": LUMPY."
lumpyexpress -B $bam -T $TMPDIR/lumpy_${sample} -o ${sample}.lumpy.vcf

mkdir -p $outdir
rsync -ta ${sample}.lumpy.vcf $outdir
