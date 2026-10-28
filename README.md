# SGAME mtDNA adaptation

Code for the variant analyses of

> Katju *et al.* No passive sidekick: mitochondrial compensatory evolution during adaptation of ETC mutant lines in *Caenorhabditis elegans*. *Communications Biology* (in revision).

Four electron transport chain (ETC) mutants of *C. elegans* (*cox-1*, *ctb-1*, *gas-1*, *isp-1*), an *isp-1; ctb-1* double mutant and mutant-free controls were evolved for 60 generations at large population size under three breeding systems (N, facultative outcrossing, wild-type N2; F, obligate outcrossing, *fog-2*; X, obligate selfing, *xol-1*). Pools of worms from each line were sequenced (Illumina NovaSeq 6000, 2 x 150 bp). These scripts map the reads, call small variants (nuclear and mitochondrial, down to 2% frequency in the mtDNA) and copy-number variants, remove mutations shared between lines, and annotate the variants unique to one line.


## Pipeline

| Step | Script | Tools | Settings |
|---|---|---|---|
| 0. Reference | [00_Reference.sh](scripts/00_Reference.sh) | bwa, samtools, GATK | WormBase WS295 N2 genome (WBcel235 assembly) |
| 1. Read quality and trimming | [01_Trimming.sh](scripts/01_Trimming.sh) | FastQC 0.11.9, Trimmomatic 0.39 | adapter clipping, leading and trailing low-quality bases removed, reads < 75 bp removed |
| 2. Mapping | [02_Mapping.sh](scripts/02_Mapping.sh) | BWA-MEM 0.7.18, samtools 1.20, Picard 3.1.1 | BWA defaults; PCR and optical duplicates tagged with MarkDuplicates |
| 3. Small variants, GATK | [03_HaplotypeCaller.sh](scripts/03_HaplotypeCaller.sh) | GATK 4 HaplotypeCaller | per-sample gvcf |
| 4. Joint genotyping | [04_JointGenotyping.sh](scripts/04_JointGenotyping.sh) | GATK 4 GenomicsDBImport, GenotypeGVCFs | all lines together |
| 5. Small variants, BCFtools | [05_Mpileup.sh](scripts/05_Mpileup.sh) | BCFtools 1.20 mpileup, call | one job per contig; on MtDNA all observed alleles kept with their read counts |
| 6. Filtering and background removal | [06_UniqueVariants.sh](scripts/06_UniqueVariants.sh), [06_unique_variants.py](scripts/06_unique_variants.py) | BCFtools, Python 3 | MQ ≥ 50, FS ≤ 20, SOR ≤ 3, per-sample DP ≥ 5; mtDNA variants down to 2% of reads; only variants carried by a single line are kept |
| 7. Copy-number variants, depth | [07_CNVnator.sh](scripts/07_CNVnator.sh) | CNVnator 0.4.1 | per sample |
| 8. Structural variants, read mapping | [08_Lumpy.sh](scripts/08_Lumpy.sh) | LUMPY (lumpyexpress) | per sample |
| 9. Annotation | [09_SnpEff.sh](scripts/09_SnpEff.sh) | SnpEff 5.2c | WBcel235 database; MtDNA with the invertebrate mitochondrial code |

The unique variants (step 6) and copy-number variants (steps 7 and 8) were then checked by eye in IGV. Only variants that passed this check were kept in the paper.

### How step 6 works

- **Site filters.** MQ, FS and SOR are read from the GATK INFO field. For BCFtools calls, FS and SOR are computed from the forward and reverse read counts summed over individual lines. On MtDNA, where bcftools makes no genotype call, the MQ ≥ 50 threshold is applied to the reads counts from in step 5.
- **Background removal.** A variant carried by more than one line is treated as standing variation of the ancestral strains or a recurrent artefact, and is removed. 
## Running

```bash
conda env create -f envs/main.yml          
conda env create -f envs/cnvnator.yml      
conda env create -f envs/lumpy.yml         

P=/path/to/project                        
READS=/path/to/fastq                      
REF=$P/genome/c_elegans.PRJNA13758.WS295.genomic.fa
cd scripts && mkdir -p logs
conda activate sgame

ref_job=$(sbatch --parsable -o logs/00_Reference.%j.out 00_Reference.sh -o $P/genome)

map_jobs=""
for r1 in $READS/*_R1.fastq.gz; do
    s=$(basename $r1 _R1.fastq.gz)
    t=$(sbatch --parsable -o logs/01_Trimming_$s.%j.out 01_Trimming.sh -1 $r1 -o $P/reads_trimmed)
    m=$(sbatch --parsable --dependency=afterok:$t:$ref_job -o logs/02_Mapping_$s.%j.out \
        02_Mapping.sh -1 $P/reads_trimmed/${s}_R1.trim.fq.gz -g $REF -o $P/bams)
    map_jobs="$map_jobs:$m"
done
```


```bash
hc_jobs=""
for bam in $P/bams/*.mkdup.bam; do
    s=$(basename $bam .mkdup.bam)
    j=$(sbatch --parsable -o logs/03_HaplotypeCaller_$s.%j.out 03_HaplotypeCaller.sh -b $bam -g $REF -o $P/gvcf)
    hc_jobs="$hc_jobs:$j"
done
jg_job=$(sbatch --parsable --dependency=afterok$hc_jobs -o logs/04_JointGenotyping.%j.out \
    04_JointGenotyping.sh -d $P/gvcf -g $REF -n SGAME_joint -o $P/gatk)

mp_jobs=""
for c in I II III IV V X MtDNA; do
    j=$(sbatch --parsable -o logs/05_Mpileup_$c.%j.out 05_Mpileup.sh -d $P/bams -g $REF -r $c -o $P/bcftools)
    mp_jobs="$mp_jobs:$j"
done

uv_job=$(sbatch --parsable --dependency=afterok:$jg_job$mp_jobs -o logs/06_UniqueVariants.%j.out \
    06_UniqueVariants.sh -v $P/gatk/SGAME_joint.vcf.gz -m $P/bcftools -g $REF -o $P/unique)
sbatch --dependency=afterok:$uv_job -o logs/09_SnpEff.%j.out 09_SnpEff.sh -v $P/unique/unique_variants.sites.vcf.gz -o $P/annotation
```

```bash
conda activate sgame-cnvnator
for bam in $P/bams/*.mkdup.bam; do
    s=$(basename $bam .mkdup.bam)
    sbatch -o logs/07_CNVnator_$s.%j.out 07_CNVnator.sh -b $bam -c $P/genome/chromosomes -o $P/cnv/cnvnator
done
conda activate sgame-lumpy
for bam in $P/bams/*.mkdup.bam; do
    s=$(basename $bam .mkdup.bam)
    sbatch -o logs/08_Lumpy_$s.%j.out 08_Lumpy.sh -b $bam -o $P/cnv/lumpy
done
```

## Main outputs

| File | Content |
|---|---|
| `unique/unique_variants.tsv` | one row per variant and line: callers, position, alleles, line, depth, alternative read count and fraction, MQ, FS, SOR |
| `unique/{gatk,bcftools}.{nuclear,MtDNA}.tsv` | the same per caller and genome, with the counts of sites removed at each filter in the slurm log |
| `annotation/unique_variants.snpeff.tsv` | most severe SnpEff effect f unique variant |
| `cnv/cnvnator/*.cnvnator.vcf`, `cnv/lumpy/*.lumpy.vcf` | copy-number and SV calls per line |

## Not scripted here

- Manual inspection of candidate variants in IGV
- 3D structural mapping of the nonsynonymous mtDNA variants with AlphaFold3
- The list of nuclear-encoded ETC subunits (Supplementary File S2), compiled from MitoCarta 3.0, OrthoList 2, WormBase WS292 and BioGRID


