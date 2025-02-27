# Binning prokaryote MAGs
Were going to do this in three iterations, only keeping medium and high quality MAGs after each iteration and re-binning all unbinned contigs again.

# BIN 1
## Make BAM files
```
mamba activate reassemble
SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    mkdir 05_BIN/${SAMPLE}

    N=${SAMPLE}_BAM

    CONTIGS=03_ASSEMBLIES/${SAMPLE}/${SAMPLE}_non_vir_plas_euk_contigs.fasta

    F=00_CLEANREADS/hocort/${SAMPLE}/${SAMPLE}_clean_1.fq.gz
    R=00_CLEANREADS/hocort/${SAMPLE}/${SAMPLE}_clean_2.fq.gz

    coverm make -r ${CONTIGS} -1 ${F} -2 ${R} -o 05_BIN/${SAMPLE} -t 16 &&\
    mv 05_BIN/${SAMPLE}/*.bam 05_BIN/BAM_1/${SAMPLE}.bam

done
```

## METABAT2
```
mamba activate metabat2

mkdir 05_BIN/BIN_1

SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in ${SAMPLE} ; do

    mkdir -p 05_BIN/BIN_1/${SAMPLE}/METABAT

    N=METABAT

    CONTIGS=03_ASSEMBLIES/${SAMPLE}/${SAMPLE}_non_vir_plas_euk_contigs.fasta

    OUT_DIR=05_BIN/BIN_1/${SAMPLE}/METABAT

    jgi_summarize_bam_contig_depths \
    05_BIN/BAM_1/${SAMPLE}.bam \
    --outputDepth 05_BIN/BIN_1/${SAMPLE}/METABAT/depth.txt &&\
    metabat2 -i ${CONTIGS} -a ${OUT_DIR}/depth.txt -o ${OUT_DIR}/BINS &&\
    gzip ${OUT_DIR}/*.fa

done
```

## Semibin
```
mamba activate Assemble_Bin

SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in ${SAMPLE} ; do

    mkdir -p 05_BIN/BIN_1/${SAMPLE}/Semibin

    N=${SAMPLE}_SemiBin

    CONTIGS=03_ASSEMBLIES/${SAMPLE}/non_viral_non_plasmid_contigs.fasta
    BAM=05_BIN/${SAMPLE}/*.bam

    SemiBin2 single_easy_bin \
    --input-fasta ${CONTIGS} --input-bam ${BAM} \
    --environment wastewater --quiet --output 05_BIN/${SAMPLE}/BIN_1/SEMIBIN

done
```

## comebin 
```
module load comebin
for SAMPLE in ${SAMPLE} ; do

    mkdir -p 05_BIN/BIN_1/${SAMPLE}/comebin

    N=comebin

    CONTIGS=03_ASSEMBLIES/${SAMPLE}/${SAMPLE}_non_vir_plas_euk_contigs.fasta

    OUT_DIR=05_BIN/BIN_1/${SAMPLE}/comebin

    run_comebin.sh \
    -a ${CONTIGS} \
    -o ${OUT_DIR} \
    -p 05_BIN/${SAMPLE} -t 32

done
```

## METADECODER
```
mamba activate metadecoder2

for SAMPLE in ${SAMPLE} ; do

    mkdir 05_BIN/BIN_1/${SAMPLE}/METADECODER
    mkdir 05_BIN/BIN_1/${SAMPLE}/METADECODER/FILES
    mkdir 05_BIN/BIN_1/${SAMPLE}/METADECODER/BINS

    CONTIGS=03_ASSEMBLIES/${SAMPLE}/${SAMPLE}_non_vir_plas_euk_contigs.fasta

    FILES=05_BIN/BIN_1/${SAMPLE}/METADECODER

    BAMPATH=05_BIN/BAM_1

    N=metadecoder

    cd ${FILES} &&\
    metadecoder coverage \
    --threads 32 \
    -b ${BAMPATH}/${SAMPLE}.bam \
    -o FILES/METADECODER.COVERAGE &&\
    metadecoder seed \
    --threads 32 \
    -f ${CONTIGS} \
    -o FILES/METADECODER.SEED &&\
    metadecoder cluster \
    -f ${CONTIGS} \
    -c FILES/METADECODER.COVERAGE \
    -s FILES/METADECODER.SEED \
    -o ${SAMPLE} &&
    gzip *.fasta && mv *.fasta.gz BINS/ && mv *.{dpgmm,kmers} FILES/

done
```

## metacoag
Prepare the bam files, this need the raw contigs.
```
mamba activate reassemble
SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    gunzip 03_ASSEMBLIES/${SAMPLE}/scaffolds.fasta.gz

    N=${SAMPLE}_covermcontig
    F_cp=00_CLEANREADS/hocort/${SAMPLE}/${SAMPLE}_clean_1.fq.gz
    R_cp=00_CLEANREADS/hocort/${SAMPLE}/${SAMPLE}_clean_2.fq.gz
    CONTIGS=03_ASSEMBLIES/${SAMPLE}/scaffolds.fasta
    ABUNDANCE=05_BIN/${SAMPLE}/abundance.tsv

    coverm contig -1 ${F_cp} -2 ${R_cp} -r ${CONTIGS} -o ${ABUNDANCE} -t 8 
    sed -i '1d' ${ABUNDANCE}

done
```
Run the binner
```
mamba activate gbintk
SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    #mkdir 05_BIN/${SAMPLE}/metacoag22222

    gunzip 03_ASSEMBLIES/${SAMPLE}/scaffolds.fasta.gz

    N=${SAMPLE}_metacoag
    CONTIGS=03_ASSEMBLIES/${SAMPLE}/scaffolds.fasta
    GRAPH=03_ASSEMBLIES/${SAMPLE}/assembly_graph_with_scaffolds.gfa
    PATHS=03_ASSEMBLIES/${SAMPLE}/scaffolds.paths
    ABUNDANCE=05_BIN/${SAMPLE}/abundance.tsv

    gbintk metacoag --assembler spades \
    --graph ${GRAPH} \
    --contigs ${CONTIGS} \
    --paths ${PATHS} \
    --abundance ${ABUNDANCE} \
    --nthreads 32 \
    --min_length 2000 \
    --output 05_BIN/${SAMPLE}/metacoag22222

done
```

## binette
```
conda activate binette
export CHECKM2DB="/blue/hlaughinghouse/flefler/databases/CheckM2_database/uniref100.KO.1.dmnd"

SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    N=${SAMPLE}_binette

    CONTIGS=03_ASSEMBLIES/${SAMPLE}/scaffolds.fasta

    binette --bin_dirs \
    05_BIN/${SAMPLE}/metacoag22222
    05_BIN/${SAMPLE}/BIN_1/SEMIBIN/output_bins \
    05_BIN/${SAMPLE}/BIN_1/METABAT \
    05_BIN/${SAMPLE}/BIN_1/metadecoder \
    05_BIN/${SAMPLE}/BIN_1/comebin/comebin_res/comebin_res_bins \
    -c ${CONTIGS} \
    -m 50 -t 8 \
    -o 05_BIN/${SAMPLE}/binette

done
```

### Quarntine the high-contam bins
```
mamba activate reassemble
for SAMPLE in ${SAMPLE} ; do

    mkdir 05_BIN/BIN_1/${SAMPLE}/binette/quarantine

    awk -F '\t' 'NR>1 && $5 > 10 {print "bin_" $1}' 05_BIN/BIN_1/${SAMPLE}/binette/final_bins_quality_reports.tsv |\
    xargs -I {} sh -c "mv 05_BIN/BIN_1/${SAMPLE}/binette/final_bins/{}.fa.gz 05_BIN/BIN_1/${SAMPLE}/binette/quarantine"

done 
```

# BIN 2
## wrangle the contigs
```
mamba activate reassemble

for SAMPLE in ${SAMPLE} ; do

    zgrep ">" 05_BIN/BIN_1/${SAMPLE}/binette/final_bins/*.fa.gz | sed 's/^.*fa.gz://' | sed 's/>//' | \
    seqkit grep -v -f - 05_WHOKARYOTE/${SAMPLE}/{archaea,bacteria,prokarya,unknown}_${SAMPLE}_non_vir_plas_euk_contigs.fasta.gz | \
    gzip > 03_ASSEMBLIES/${SAMPLE}/unbinned_1.fasta.gz


done
```

## Make BAM files
```
mamba activate reassemble
SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    mkdir 05_BIN/${SAMPLE}

    N=${SAMPLE}_BAM

    CONTIGS=03_ASSEMBLIES/${SAMPLE}/unbinned_1.fasta.gz

    F=00_CLEANREADS/hocort/${SAMPLE}/${SAMPLE}_clean_1.fq.gz
    R=00_CLEANREADS/hocort/${SAMPLE}/${SAMPLE}_clean_2.fq.gz

    coverm make -r ${CONTIGS} -1 ${F} -2 ${R} -o 05_BIN/${SAMPLE} -t 16 &&\
    mv 05_BIN/${SAMPLE}/*.bam 05_BIN/BAM_2/${SAMPLE}.bam

done
```

## METABAT2
```
mamba activate metabat2

mkdir 05_BIN/BIN_2

SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in ${SAMPLE} ; do

    mkdir -p 05_BIN/BIN_2/${SAMPLE}/METABAT

    N=METABAT

    CONTIGS=03_ASSEMBLIES/${SAMPLE}/unbinned_1.fasta.gz

    OUT_DIR=05_BIN/BIN_2/${SAMPLE}/METABAT

    jgi_summarize_bam_contig_depths \
    05_BIN/BAM_2/${SAMPLE}.bam \
    --outputDepth 05_BIN/BIN_2/${SAMPLE}/METABAT/depth.txt &&\
    metabat2 -i ${CONTIGS} -a ${OUT_DIR}/depth.txt -o ${OUT_DIR}/BINS &&\
    gzip ${OUT_DIR}/*.fa

done
```

## Semibin
```
mamba activate Assemble_Bin

SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in ${SAMPLE} ; do

    mkdir -p 05_BIN/BIN_2/${SAMPLE}/Semibin

    N=${SAMPLE}_SemiBin

    CONTIGS=03_ASSEMBLIES/${SAMPLE}/non_viral_non_plasmid_contigs.fasta
    BAM=05_BIN/${SAMPLE}/*.bam

    SemiBin2 single_easy_bin \
    --input-fasta ${CONTIGS} --input-bam ${BAM} \
    --environment wastewater --quiet --output 05_BIN/${SAMPLE}/BIN_2/SEMIBIN

done
```

## comebin 
```
module load comebin
for SAMPLE in ${SAMPLE} ; do

    mkdir -p 05_BIN/BIN_2/${SAMPLE}/comebin

    N=comebin

    CONTIGS=03_ASSEMBLIES/${SAMPLE}/unbinned_1.fasta.gz

    OUT_DIR=05_BIN/BIN_2/${SAMPLE}/comebin

    run_comebin.sh \
    -a ${CONTIGS} \
    -o ${OUT_DIR} \
    -p 05_BIN/${SAMPLE} -t 32

done
```

## METADECODER
```
mamba activate metadecoder2

for SAMPLE in ${SAMPLE} ; do

    mkdir 05_BIN/BIN_2/${SAMPLE}/METADECODER
    mkdir 05_BIN/BIN_2/${SAMPLE}/METADECODER/FILES
    mkdir 05_BIN/BIN_2/${SAMPLE}/METADECODER/BINS

    CONTIGS=03_ASSEMBLIES/${SAMPLE}/unbinned_1.fasta.gz

    FILES=05_BIN/BIN_2/${SAMPLE}/METADECODER

    BAMPATH=05_BIN/BAM_2

    N=metadecoder

    cd ${FILES} &&\
    metadecoder coverage \
    --threads 32 \
    -b ${BAMPATH}/${SAMPLE}.bam \
    -o FILES/METADECODER.COVERAGE &&\
    metadecoder seed \
    --threads 32 \
    -f ${CONTIGS} \
    -o FILES/METADECODER.SEED &&\
    metadecoder cluster \
    -f ${CONTIGS} \
    -c FILES/METADECODER.COVERAGE \
    -s FILES/METADECODER.SEED \
    -o ${SAMPLE} &&
    gzip *.fasta && mv *.fasta.gz BINS/ && mv *.{dpgmm,kmers} FILES/

done
```

## binette
```
conda activate binette
export CHECKM2DB="/blue/hlaughinghouse/flefler/databases/CheckM2_database/uniref100.KO.1.dmnd"

SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    N=${SAMPLE}_binette

    CONTIGS=03_ASSEMBLIES/${SAMPLE}/scaffolds.fasta

    binette --bin_dirs \
    05_BIN/${SAMPLE}/BIN_2/SEMIBIN/output_bins \
    05_BIN/${SAMPLE}/BIN_2/METABAT \
    05_BIN/${SAMPLE}/BIN_2/metadecoder \
    05_BIN/${SAMPLE}/BIN_2/comebin/comebin_res/comebin_res_bins \
    -c ${CONTIGS} \
    -m 50 -t 8 \
    -o 05_BIN/${SAMPLE}/binette

done
```

### Quarntine the high-contam bins
```
mamba activate reassemble
for SAMPLE in ${SAMPLE} ; do

    mkdir 05_BIN/BIN_2/${SAMPLE}/binette/quarantine

    awk -F '\t' 'NR>1 && $5 > 10 {print "bin_" $1}' 05_BIN/BIN_2/${SAMPLE}/binette/final_bins_quality_reports.tsv |\
    xargs -I {} sh -c "mv 05_BIN/BIN_2/${SAMPLE}/binette/final_bins/{}.fa.gz 05_BIN/BIN_2/${SAMPLE}/binette/quarantine"

done 
```

# BIN 3
## wrangle the contigs
```
mamba activate reassemble

for SAMPLE in ${SAMPLE} ; do

    zgrep ">" 05_BIN/BIN_2/${SAMPLE}/binette/final_bins/*.fa.gz | sed 's/^.*fa.gz://' | sed 's/>//' | \
    seqkit grep -v -f - 05_WHOKARYOTE/${SAMPLE}/{archaea,bacteria,prokarya,unknown}_${SAMPLE}_non_vir_plas_euk_contigs.fasta.gz | \
    gzip > 03_ASSEMBLIES/${SAMPLE}/unbinned_2.fasta.gz


done
```

## Make BAM files
```
mamba activate reassemble
SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    mkdir 05_BIN/${SAMPLE}

    N=${SAMPLE}_BAM

    CONTIGS=03_ASSEMBLIES/${SAMPLE}/unbinned_2.fasta.gz

    F=00_CLEANREADS/hocort/${SAMPLE}/${SAMPLE}_clean_1.fq.gz
    R=00_CLEANREADS/hocort/${SAMPLE}/${SAMPLE}_clean_2.fq.gz

    coverm make -r ${CONTIGS} -1 ${F} -2 ${R} -o 05_BIN/${SAMPLE} -t 16 &&\
    mv 05_BIN/${SAMPLE}/*.bam 05_BIN/BAM_3/${SAMPLE}.bam

done
```

## METABAT2
```
mamba activate metabat2

mkdir 05_BIN/BIN_3

SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in ${SAMPLE} ; do

    mkdir -p 05_BIN/BIN_3/${SAMPLE}/METABAT

    N=METABAT

    CONTIGS=03_ASSEMBLIES/${SAMPLE}/unbinned_2.fasta.gz

    OUT_DIR=05_BIN/BIN_3/${SAMPLE}/METABAT

    jgi_summarize_bam_contig_depths \
    05_BIN/BAM_3/${SAMPLE}.bam \
    --outputDepth 05_BIN/BIN_3/${SAMPLE}/METABAT/depth.txt &&\
    metabat2 -i ${CONTIGS} -a ${OUT_DIR}/depth.txt -o ${OUT_DIR}/BINS &&\
    gzip ${OUT_DIR}/*.fa

done
```

## Semibin
```
mamba activate Assemble_Bin

SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in ${SAMPLE} ; do

    mkdir -p 05_BIN/BIN_3/${SAMPLE}/Semibin

    N=${SAMPLE}_SemiBin

    CONTIGS=03_ASSEMBLIES/${SAMPLE}/non_viral_non_plasmid_contigs.fasta
    BAM=05_BIN/${SAMPLE}/*.bam

    SemiBin2 single_easy_bin \
    --input-fasta ${CONTIGS} --input-bam ${BAM} \
    --environment wastewater --quiet --output 05_BIN/${SAMPLE}/BIN_3/SEMIBIN

done
```

## comebin 
```
module load comebin
for SAMPLE in ${SAMPLE} ; do

    mkdir -p 05_BIN/BIN_3/${SAMPLE}/comebin

    N=comebin

    CONTIGS=03_ASSEMBLIES/${SAMPLE}/unbinned_2.fasta.gz

    OUT_DIR=05_BIN/BIN_3/${SAMPLE}/comebin

    run_comebin.sh \
    -a ${CONTIGS} \
    -o ${OUT_DIR} \
    -p 05_BIN/${SAMPLE} -t 32

done
```

## METADECODER
```
mamba activate metadecoder2

for SAMPLE in ${SAMPLE} ; do

    mkdir 05_BIN/BIN_3/${SAMPLE}/METADECODER
    mkdir 05_BIN/BIN_3/${SAMPLE}/METADECODER/FILES
    mkdir 05_BIN/BIN_3/${SAMPLE}/METADECODER/BINS

    CONTIGS=03_ASSEMBLIES/${SAMPLE}/unbinned_2.fasta.gz

    FILES=05_BIN/BIN_3/${SAMPLE}/METADECODER

    BAMPATH=05_BIN/BAM_3

    N=metadecoder

    cd ${FILES} &&\
    metadecoder coverage \
    --threads 32 \
    -b ${BAMPATH}/${SAMPLE}.bam \
    -o FILES/METADECODER.COVERAGE &&\
    metadecoder seed \
    --threads 32 \
    -f ${CONTIGS} \
    -o FILES/METADECODER.SEED &&\
    metadecoder cluster \
    -f ${CONTIGS} \
    -c FILES/METADECODER.COVERAGE \
    -s FILES/METADECODER.SEED \
    -o ${SAMPLE} &&
    gzip *.fasta && mv *.fasta.gz BINS/ && mv *.{dpgmm,kmers} FILES/

done
```

## binette
```
conda activate binette
export CHECKM2DB="/blue/hlaughinghouse/flefler/databases/CheckM2_database/uniref100.KO.1.dmnd"

SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    N=${SAMPLE}_binette

    CONTIGS=03_ASSEMBLIES/${SAMPLE}/scaffolds.fasta

    binette --bin_dirs \
    05_BIN/${SAMPLE}/BIN_3/SEMIBIN/output_bins \
    05_BIN/${SAMPLE}/BIN_3/METABAT \
    05_BIN/${SAMPLE}/BIN_3/metadecoder \
    05_BIN/${SAMPLE}/BIN_3/comebin/comebin_res/comebin_res_bins \
    -c ${CONTIGS} \
    -m 50 -t 8 \
    -o 05_BIN/${SAMPLE}/binette

done
```

### Quarntine the high-contam bins
```
mamba activate reassemble
for SAMPLE in ${SAMPLE} ; do

    mkdir 05_BIN/BIN_3/${SAMPLE}/binette/quarantine

    awk -F '\t' 'NR>1 && $5 > 10 {print "bin_" $1}' 05_BIN/BIN_3/${SAMPLE}/binette/final_bins_quality_reports.tsv |\
    xargs -I {} sh -c "mv 05_BIN/BIN_3/${SAMPLE}/binette/final_bins/{}.fa.gz 05_BIN/BIN_3/${SAMPLE}/binette/quarantine"

done 
```
