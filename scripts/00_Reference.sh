#!/bin/bash
#SBATCH -J Reference
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 2:00:00

# WormBase WS295 N2 reference bwa, samtools and GATK
##sbatch 00_Reference.sh -o /path/to/genome_outdir

while getopts o: flag
do
case "${flag}" in
                o) outdir=${OPTARG}
                        ;;
                *) echo "Invalid option: -$flag"
                	;;
        esac
done

[[ -z "$outdir" ]] && { echo "Missing argument: outdir"; exit 1; }

mkdir -p $outdir
cd $outdir

genome=c_elegans.PRJNA13758.WS295.genomic.fa
echo "## ------------------------------------------------ ##"
gunzip -k ${genome}.gz

echo "## ------------------------------------------------ ##"
echo `date` ": indexing."
bwa index $genome
samtools faidx $genome
gatk CreateSequenceDictionary -R $genome
mkdir -p chromosomes
for c in `cut -f1 ${genome}.fai`; do samtools faidx $genome $c > chromosomes/${c}.fa; done
