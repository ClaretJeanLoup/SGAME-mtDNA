#!/bin/bash
#SBATCH -J Trimming
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 4:00:00

# Checks read quality of one paired-end sample with FastQC, trims adapters and low-quality ends with Trimmomaticcchecks the trimmed reads.
##sbatch 01_Trimming.sh -1 /path/to/sample_R1.fastq.gz -o /path/to/trimmed_outdir [-a adapters.fa]

adapters=$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE-2.fa

while getopts 1:o:a: flag
doß
case "${flag}" in
                1) R1=${OPTARG}
                	;;
                o) outdir=${OPTARG}
                        ;;
                a) adapters=${OPTARG}
                        ;;
                *) echo "Invalid option: -$flag"
                	;;
        esac
done

for v in R1 outdir; do [[ -z "${!v}" ]] && { echo "Missing argument: $v"; exit 1; }; done
[[ -e "$R1" ]] || { echo "Not found: $R1"; exit 1; }
[[ -e "$adapters" ]] || { echo "Not found: $adapters (give the adapter fasta with -a)"; exit 1; }

R2=${R1/_R1/_R2}
[[ -e "$R2" ]] || { echo "Not found: $R2"; exit 1; }

echo "## ------------------------------------------------ ##"
echo `date` ": sync files."
rsync -ta $R1 $R2 $TMPDIR
cd $TMPDIR

r1=`basename $R1`
r2=`basename $R2`
sample=${r1%%_R1*}
echo "Sample $sample"
mkdir -p qc

echo "## ------------------------------------------------ ##"
echo `date` ": quality of raw reads."
fastqc -t 2 -o qc $r1 $r2

echo "## ------------------------------------------------ ##"
echo `date` ": trimming."
#trailing bases below quality 3 are cut and reads shorter than 75 bp are removed
trimmomatic PE -threads 8 $r1 $r2 \
    ${sample}_R1.trim.fq.gz ${sample}_R1.unpaired.fq.gz ${sample}_R2.trim.fq.gz ${sample}_R2.unpaired.fq.gz \
    ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 MINLEN:75 2> qc/${sample}.trimmomatic.log

echo "## ------------------------------------------------ ##"
echo `date` ": quality of trimmed reads."
fastqc -t 2 -o qc ${sample}_R1.trim.fq.gz ${sample}_R2.trim.fq.gz

mkdir -p $outdir/qc
rsync -ta ${sample}_R1.trim.fq.gz ${sample}_R2.trim.fq.gz $outdir
rsync -ta qc/ $outdir/qc
