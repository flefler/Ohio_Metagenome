# Prokaryotic bins

## Move the good Prok Bins
```
mkdir 07_MAGs
SAMPLES=`cut -f 1 samples_1.txt | sed '1d'`
for SAMPLE in $SAMPLES; do 
    mkdir -p 07_MAGs/${SAMPLE}/MAGs
    for file in 05_BIN/${SAMPLE}/{BIN_1,BIN_2,BIN_3}/binette/final_bins/*.fa.gz; do 
        cp "$file" "07_MAGs/${SAMPLE}/MAGs/${SAMPLE}_$(basename "$file" .fa.gz).fasta.gz"; 
    done
done
```

## Simple stats
```
SAMPLES=`cut -f 1 samples_1.txt | sed '1d'`
mamba activate reassemble
for SAMPLE in $SAMPLES; do 
    mkdir 07_MAGs/${SAMPLE}/DATA
    seqkit stats -a 07_MAGs/${SAMPLE}/MAGs/*.fasta.gz -T | csvtk tab2csv | csvtk cut -f -2,-3,-9,-10,-11,-12,-15,-16,-17 -o 07_MAGs/${SAMPLE}/DATA/${SAMPLE}_INFO.csv
done
```

## compareM2
```
mamba activate comparem2
SAMPLES=`cut -f 1 samples_1.txt | sed '1d'`
for SAMPLE in $SAMPLES; do 

    N=${SAMPLE}_comparem2

    cd 07_MAGs/${SAMPLE} &&/
    comparem2 --config input_genomes=MAGs/*.fasta.gz \
    output_directory=DATA/COMPAREM2 \
    --until sequence_lengths checkm2 gtdbtk

done
```

## CoverM
```
mamba activate reassemble
SAMPLES=`cut -f 1 samples_1.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    N=${SAMPLE}_coverM

    coverm genome \
    --read1 00_CLEANREADS/hocort/${SAMPLE}/${SAMPLE}_clean_1.fq.gz \
    --read2 00_CLEANREADS/hocort/${SAMPLE}/${SAMPLE}_clean_2.fq.gz \
    --genome-fasta-directory 07_MAGs/${SAMPLE}/MAGs \
    -x .gz --threads 8 --methods relative_abundance \
    -o 07_MAGs/${SAMPLE}/DATA/${SAMPLE}_coverM.tsv

done
```

# Viral Bins
I made a list, manually, of the viral contigs which were circular based on geNomad called ```vMAGs_2.txt```

## Move and rename those circular viral bins
```
mkdir 08_vMAGs
seqkit grep -f vMAGs_2.txt 03_ASSEMBLIES/${SAMPLE}/cleaned_scaffolds.fasta.gz | \
seqkit replace -p "^(\S+)_(\S+).*" -r "${SAMPLE}_vMAG_\${1}_\${2}" | \
seqkit split --by-id -O 08_vMAGs | \
sed -i 's/^>.*/>LE_vMAG_NODE_/' 08_vMAGs/*.fasta

for file in 08_vMAGs/*.fasta; do
    new_name=$(head -n1 "$file" | cut -d' ' -f1 | sed 's/>//' | cut -d'_' -f1-4).fasta.gz
    gzip -c "$file" > "08_vMAGs/$new_name"
    rm "$file"
done
```

## Run CheckV
```
module load checkv
for file in 08_vMAGs/*.fasta.gz; do
    checkv end_to_end ${file} 03_INFO/${file}
done
```

# Annotate all MAGs
```
mkdir 09_ANNOTATE
cp {07_MAGs,08_vMAGs}/* 04_ANNOTATE/
gunzip 09_ANNOTATE/*.fasta.gz

mamba activate metacerberus

metacerberus.py \
--illumina \
09_ANNOTATE/* \
--dir_out 09_ANNOTATE/METACERBERUS \
--hmm ALL
```

