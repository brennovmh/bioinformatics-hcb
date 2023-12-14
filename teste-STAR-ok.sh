#!/bin/bash

### LEMBRAR DE DESCOMPACTAR ARQUIVOS

cd /home/hcb/rawData/HCBagilent1/FASTQcat/teste/paired
gunzip *fastq.gz

diretorio="./STARresults"

# Verifica se o diretório existe
if [ -d "$diretorio" ]; then
    echo "O diretório $diretorio já existe."
else
    # Cria o diretório se não existir
    mkdir -p "$diretorio"
    echo "Diretório $diretorio criado com sucesso."
fi


# Path to the genome index directory
genome_index="/home/hcb/tools/ref/ref_STAR/"

# Path to the directory containing all the datasets
data_dir="/home/hcb/rawData/HCBagilent1/FASTQcat/teste/paired"

# Output directory
output_dir="/home/hcb/rawData/HCBagilent1/FASTQcat/teste/paired/STARresults"

# Iterate over each dataset in the data directory
# Iterate over each R1 file
for r1_file in "${data_dir}"/*_R1.fastq; do
    if [ -e "$r1_file" ]; then
        # Construct the corresponding R2 filename based on the R1 filename
        r2_file="${r1_file/_R1/_R2}"

        # Check if the corresponding R2 file exists
        if [ -e "$r2_file" ]; then
            # Extract dataset name from the R1 filename
            dataset_name=$(basename "$r1_file" _R1.fq)

            # Run STAR aligner for the current pair
            STAR \
                --runMode alignReads \
                --outSAMtype BAM SortedByCoordinate \
                --genomeDir ${genome_index} \
                --readFilesIn ${r1_file} ${r2_file} \
                --outFileNamePrefix ${output_dir}/${dataset_name}_
        else
            echo "Corresponding R2 file not found for $r1_file"
        fi
    fi
done




