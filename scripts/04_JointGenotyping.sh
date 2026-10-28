#!/bin/bash
#SBATCH -J JointGenotyping
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 3-00:00:00

#  genotypes all lines with GATK GenotypeGVCFs.
##sbatch 04_JointGenotyping.sh -d /path/to/gvcf_dir -g /path/to/ref.fa -n outname -o /path/to/outdir

while getopts d:o:n:g: flag
do
case "${flag}" in
                d) DIRDATA=${OPTARG}
                    	;;
                o) outdir=${OPTARG}
                        ;;
                n) name=${OPTARG}
			            ;;
                g) ref=${OPTARG}
			            ;;
                *) echo "Invalid option: -$flag"
                      	;;
        esac
done

for v in DIRDATA ref name outdir; do [[ -z "${!v}" ]] && { echo "Missing argument: $v"; exit 1; }; done
[[ -d "$DIRDATA" ]] || { echo "Not found: $DIRDATA"; exit 1; }
ls $DIRDATA/*.g.vcf.gz >/dev/null 2>&1 || { echo "No *.g.vcf.gz in $DIRDATA"; exit 1; }
[[ -e "$ref" ]] || { echo "Not found: $ref"; exit 1; }

echo "## ------------------------------------------------ ##"
echo `date` ": sync files."
rsync -ta $DIRDATA/*.g.vcf.gz* $TMPDIR
rsync -ta ${ref} ${ref}.fai ${ref%.*}.dict $TMPDIR
cd $TMPDIR

genome=`basename $ref`

rm -f cohort.sample_map
for i in *.g.vcf.gz; do
    if [[ ! -f "${i}.tbi" ]]; then gatk IndexFeatureFile -I $i; fi
    echo -e ${i%.g.vcf.gz}'\t'./${i} >> cohort.sample_map
done

# One interval per contig, taken from the fasta index.
awk '{print $1"\t0\t"$2}' ${genome}.fai > contigs.bed

echo "## ------------------------------------------------ ##"
echo `date` ": building GenomicsDB workspace."
gatk --java-options "-Xmx100G" GenomicsDBImport --genomicsdb-workspace-path ./my_database \
    --batch-size 20 --sample-name-map cohort.sample_map -L contigs.bed --tmp-dir $TMPDIR

echo "## ------------------------------------------------ ##"
echo `date` ": genotyping GenomicsDB workspace."
gatk --java-options "-Xmx100G" GenotypeGVCFs -V gendb://my_database -R $genome -O ${name}.vcf.gz --tmp-dir $TMPDIR

cd -
mkdir -p $outdir
rsync -ta $TMPDIR/${name}.vcf.gz* $outdir
