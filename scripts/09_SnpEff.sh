#!/bin/bash
#SBATCH -J SnpEff
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 2:00:00
##sbatch 09_SnpEff.sh -v /path/to/unique_variants.sites.vcf.gz -o /path/to/outdir [-d snpeff_database]

db=WBcel235.105

while getopts v:o:d: flag
do
case "${flag}" in
                v) VCF=${OPTARG}
                    	;;
                o) outdir=${OPTARG}
                        ;;
                d) db=${OPTARG}
                        ;;
                *) echo "Invalid option: -$flag"
                      	;;
        esac
done

for v in VCF outdir; do [[ -z "${!v}" ]] && { echo "Missing argument: $v"; exit 1; }; done
[[ -e "$VCF" ]] || { echo "Not found: $VCF"; exit 1; }

mkdir -p $outdir
cd $outdir

echo "## ------------------------------------------------ ##"
echo `date` ": SnpEff on $db."
snpEff download -v $db
# MtDNA genes  table 5
snpEff -Xmx8g -v -configOption ${db}.MtDNA.codonTable=Invertebrate_Mitochondrial -stats unique_variants.snpeff.html \
    $db $VCF | bgzip > unique_variants.snpeff.vcf.gz
bcftools index -t unique_variants.snpeff.vcf.gz

#First ANN entry is the most severe effect
echo -e "CHROM\tPOS\tREF\tALT\tLINE\tAF_LINE\tEFFECT\tIMPACT\tGENE\tHGVS_C\tHGVS_P" > unique_variants.snpeff.tsv
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/LINE\t%INFO/AF_LINE\t%INFO/ANN\n' unique_variants.snpeff.vcf.gz \
    | awk -F'\t' -v OFS='\t' '{split($7, a, ","); split(a[1], f, "|"); print $1, $2, $3, $4, $5, $6, f[2], f[3], f[4], f[10], f[11]}' \
    >> unique_variants.snpeff.tsv
echo "Annotated: `tail -n +2 unique_variants.snpeff.tsv | wc -l` variants"
