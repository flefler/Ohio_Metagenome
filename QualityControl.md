# Quality Control Script
Last updated Feb 27 2025 by Forrest W. Lefler
Sample names are in samples.txt

## fastQC
Run fastQC on raw data
```
SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in $SAMPLES; do
    fastqc --quiet 00_Reads/${SAMPLE}/*.fq.gz -o 01_QC/FASTQC
done
```

## fastP
We want to deduplicate the reads to reduce redundancy in assembly, less is more. we also want the -g option as these were sequenced with NovaSeq
```
mamba activate hocort
SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    N=${SAMPLE}_fastp
    mkdir 00_CLEANREADS/fastP
    mkdir 00_CLEANREADS/fastP/HTML
    mkdir 00_CLEANREADS/fastP/JSON
    mkdir 00_CLEANREADS/fastP/${SAMPLE}

    F=00_Reads/${SAMPLE}/*_1.fq.gz
    R=00_Reads/${SAMPLE}/*_2.fq.gz
    F_cp=00_CLEANREADS/fastP/${SAMPLE}/${SAMPLE}_clean_1.fq.gz
    R_cp=00_CLEANREADS/fastP/${SAMPLE}/${SAMPLE}_clean_2.fq.gz
    HTML=01_QC/fastP/HTML/${SAMPLE}.html
    JSON=01_QC/fastP/JSON/${SAMPLE}.json

    fastp -i ${F} -I ${R} -o ${F_cp} -O ${R_cp} -D -g -j ${JSON} -h ${HTML} -w 16

done
```

## hocort
Remove human, mouse, and phix reads with bowtie2, remove low quality and duplicate reads with fastP
fastP doesnt need to be run twice, i just wasted time and created more files for no good reason.
```
mamba activate hocort
SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    mkdir 00_CLEANREADS/hocort
    mkdir 00_CLEANREADS/hocort/${SAMPLE}

    N=${SAMPLE}_decontam
    F_cp=00_CLEANREADS/fastP/${SAMPLE}/${SAMPLE}_clean_1.fq.gz
    R_cp=00_CLEANREADS/fastP/${SAMPLE}/${SAMPLE}_clean_2.fq.gz
    REF=/blue/hlaughinghouse/flefler/contam_databases/homo_mus_phix/homo_mus_phix
    F_ch=00_CLEANREADS/hocort/${SAMPLE}/${SAMPLE}_clean_1.fq.gz
    R_ch=00_CLEANREADS/hocort/${SAMPLE}/${SAMPLE}_clean_2.fq.gz

    hocort map bowtie2 -x ${REF} -i ${F_cp} ${R_cp} -o ${F_ch} ${R_ch} -c=--very-fast
done
```

## fastQC
On the quality controled data
```
mamba activate Assemble_Bin
mkdir 01_QC/FASTQC

SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    N=${SAMPLE}_FASTQC

    fastqc --quiet 00_CLEANREADS/hocort/${SAMPLE}/*.fq.gz -o 01_QC/FASTQC

done
```
## MultiQC
Visualize everything
```
multiqc 01_QC -o 01_QC --interactive
```
