# Assemble and Identify
We will assemble raw reads into contigs with metaSPAdes and build report with quast and multiQC. Then we will identify viral and plasmid contigs using geNomad, and prokaryotic and eukaryotic contigs using whokaryote

## Assembly
Assemble with SPAdes, then remove contigs shorter than 1000bp and gzip all the fasta files.
```
mamba activate Assemble_Bin
mkdir 03_ASSEMBLIES

SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    N=${SAMPLE}_metaSPAdes

    F_cp=00_CLEANREADS/hocort/${SAMPLE}/${SAMPLE}_clean_1.fq.gz
    R_cp=00_CLEANREADS/hocort/${SAMPLE}/${SAMPLE}_clean_2.fq.gz

    spades.py --meta -1 ${F_cp} -2 ${R_cp} -o 03_ASSEMBLIES/${SAMPLE} -m 250 -t 32 &&\ 
    seqkit seq -m 1000 03_ASSEMBLIES/${SAMPLE}/scaffolds.fasta -o 03_ASSEMBLIES/${SAMPLE}/cleaned_scaffolds.fasta.gz &&\
    gunzip 03_ASSEMBLIES/${SAMPLE}/*.fasta &&\
    gunzip 03_ASSEMBLIES/${SAMPLE}/*/*.fasta

done
```

### Clean up SPAdes output
Spades produces a lot of large files we dont need to keep, lets get rid of them.
```
SAMPLES=`cut -f 1 samples.txt | sed '1d'`

for SAMPLE in $SAMPLES; do
    rm -r 03_ASSEMBLIES/${SAMPLE}/{corrected,K55,K33,K21,misc,pipeline_state,tmp} 
    rm 03_ASSEMBLIES/${SAMPLE}/{dataset.info,*.yaml,contigs.*,first_pe_contigs.*,before_rr.*,assembly_graph_after_simplification.gfa,strain_graph.gfa}
    gzip 03_ASSEMBLIES/${SAMPLE}/assembly_graph.fastg
done
```

## Assess Assemblies with Quast
```
mkdir 01_QC/QUAST
mamba activate Assemble_Bin

SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    N=${SAMPLE}_quast

        metaquast 03_ASSEMBLIES/${SAMPLE}/cleaned_scaffolds.fasta.gz \
        --output-dir 01_QC/QUAST/${SAMPLE} \
        --max-ref-number 0 -L \
        --threads 8 --fast --silent

done
```

## Identify viral and plasmid contigs with geNomad
```
mamba activate genomad
mkdir 04_geNomad

SAMPLES=`cut -f 1 samples_1.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    N=${SAMPLE}_geNomad
    
    genomad end-to-end --cleanup --disable-find-proviruses --quiet --enable-score-calibration --composition metagenome --threads 16 \
    --min-score 0.8 --min-virus-hallmarks 3 --min-plasmid-hallmarks 3 --min-plasmid-hallmarks-short-seqs 3 --min-virus-hallmarks-short-seqs 3 \
    03_ASSEMBLIES/${SAMPLE}/cleaned_scaffolds.fasta.gz 04_geNomad/${SAMPLE} /blue/hlaughinghouse/flefler/genomad_db

done
```
### Extract viral and plasmid contigs
```
SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    seqkit grep -v -f <(cat 04_geNomad/${SAMPLE}/cleaned_scaffolds_summary/cleaned_scaffolds_{virus,plasmid}_summary.tsv \
    | cut -f 1 | grep -v "seq_name" | sort | uniq) 03_ASSEMBLIES/${SAMPLE}/cleaned_scaffolds.fasta \
    -o 03_ASSEMBLIES/${SAMPLE}/${SAMPLE}_non_vir_plas_contigs.fasta
    
done
```

## Identify eukaryotic and prokaryotic contigs with whokaryote
```
mkdir 06_WHOKARYOTE

conda activate whokaryote

SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    N=WHOKARYOTE_${SAMPLE}
    whokaryote.py --contigs 03_ASSEMBLIES/${SAMPLE}/non_viral_non_plasmid_contigs.fasta --outdir 06_WHOKARYOTE/${SAMPLE} --minsize 1000 --threads 16

done
```

### Extract eukaryotic contigs
```
mamba activate reassemble
SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    seqkit grep -v -f 06_WHOKARYOTE/${SAMPLE}/eukaryote_contig_headers.txt 03_ASSEMBLIES/${SAMPLE}/non_viral_non_plasmid_contigs.fasta.gz \
    -o 03_ASSEMBLIES/${SAMPLE}/non_vir_plas_euk_contigs.fasta.gz
  
done
```

## Now we are left with only prokaryotic contigs in the ```non_vir_plas_euk_contigs.fasta.gz``` file
