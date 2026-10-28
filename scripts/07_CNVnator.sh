#!/bin/bash
#SBATCH -J CNVnator
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 12:00:00

# CNVnator rujn
##sbatch 07_CNVnator.sh -b /path/to/sample.mkdup.bam -c /path/to/genome/chromosomes -o /path/to/outdir [-w bin_size]

bin=100

while getopts b:c:o:w: flag
do
case "${flag}" in
                b) BAM=${OPTARG}
                	;;
                c) chromdir=${OPTARG}
                        ;;
                o) outdir=${OPTARG}
                        ;;
                w) bin=${OPTARG}
                        ;;
                *) echo "Invalid option: -$flag"
                	;;
        esac
done

for v in BAM chromdir outdir; do [[ -z "${!v}" ]] && { echo "Missing argument: $v"; exit 1; }; done
[[ -e "$BAM" ]] || { echo "Not found: $BAM"; exit 1; }
ls $chromdir/*.fa >/dev/null 2>&1 || { echo "No per-chromosome fasta in $chromdir (run 00_Reference.sh)"; exit 1; }

echo "## ------------------------------------------------ ##"
echo `date` ": sync files."
rsync -ta $BAM $BAM.bai $TMPDIR
rsync -ta $chromdir $TMPDIR
cd $TMPDIR

bam=`basename $BAM`
sample=${bam%.mkdup.bam}
chroms=`ls chromosomes/*.fa | xargs -n1 basename | sed 's/\.fa$//' | tr '\n' ' '`
echo "Sample $sample, bin $bin bp"

echo "## ------------------------------------------------ ##"
echo `date` ": CNVnator."
cnvnator -root ${sample}.root -tree $bam -chrom $chroms
cnvnator -root ${sample}.root -his $bin -d chromosomes
cnvnator -root ${sample}.root -stat $bin
cnvnator -root ${sample}.root -partition $bin
cnvnator -root ${sample}.root -call $bin > ${sample}.cnvnator.txt
cnvnator2VCF.pl -prefix $sample ${sample}.cnvnator.txt chromosomes > ${sample}.cnvnator.vcf

mkdir -p $outdir
rsync -ta ${sample}.cnvnator.txt ${sample}.cnvnator.vcf $outdir
