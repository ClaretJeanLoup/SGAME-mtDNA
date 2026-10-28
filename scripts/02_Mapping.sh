#!/bin/bash
#SBATCH -J Mapping
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 12:00:00

# Maps one trimmed sample with BWA-MEM, sorts and indexes the bam and tags PCR and optical duplicates with Picard MarkDuplicates.
##sbatch 02_Mapping.sh -1 /path/to/sample_R1.trim.fq.gz -g /path/to/ref.fa -o /path/to/bam_outdir

while getopts 1:g:o: flag
do
case "${flag}" in
                1) R1=${OPTARG}
                	;;
                g) ref=${OPTARG}
                        ;;
                o) outdir=${OPTARG}
                        ;;
                *) echo "Invalid option: -$flag"
                	;;
        esac
done

for v in R1 ref outdir; do [[ -z "${!v}" ]] && { echo "Missing argument: $v"; exit 1; }; done
[[ -e "$R1" ]] || { echo "Not found: $R1"; exit 1; }
[[ -e "$ref" ]] || { echo "Not found: $ref"; exit 1; }
[[ -e "${ref}.bwt" ]] || { echo "No bwa index for $ref (run 00_Reference.sh)"; exit 1; }

R2=${R1/_R1.trim/_R2.trim}

echo "## ------------------------------------------------ ##"
echo `date` ": sync files."
rsync -ta $R1 $R2 $TMPDIR
rsync -ta ${ref}* ${ref%.*}.dict $TMPDIR
cd $TMPDIR

genome=`basename $ref`
r1=`basename $R1`
r2=`basename $R2`
sample=${r1%%_R1*}
pu=`zcat $r1 | head -1 | cut -d: -f3,4 | tr ':' '.'`
echo "Sample $sample, PU $pu, reference $genome"

echo "## ------------------------------------------------ ##"
echo `date` ": mapping."
bwa mem -t 16 -R "@RG\tID:${sample}\tSM:${sample}\tLB:lib_${sample}\tPL:ILLUMINA\tPU:${pu}" $genome $r1 $r2 \
    | samtools sort -@ 4 -m 2G -T $TMPDIR/${sample}.sort -o ${sample}.sorted.bam
samtools index ${sample}.sorted.bam

echo "## ------------------------------------------------ ##"
echo `date` ": marking duplicates."
# Duplicates are tagged, not removed.
picard -Xmx32G MarkDuplicates I=${sample}.sorted.bam O=${sample}.mkdup.bam M=${sample}.mkdup.metrics.txt \
    TMP_DIR=$TMPDIR
samtools index ${sample}.mkdup.bam
samtools quickcheck ${sample}.mkdup.bam && samtools flagstat ${sample}.mkdup.bam > ${sample}.flagstat.txt

mkdir -p $outdir/qc
rsync -ta ${sample}.mkdup.bam ${sample}.mkdup.bam.bai $outdir
rsync -ta ${sample}.mkdup.metrics.txt ${sample}.flagstat.txt $outdir/qc
