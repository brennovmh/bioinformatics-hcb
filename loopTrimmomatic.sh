#!/bin/bash

# Diretorio onde estão os arquivos concatenados
cd /home/hcb/rawData/HCBagilent1/FASTQcat/teste

diretorio1="./paired"
diretorio2="./unpaired"

# Verifica se o diretório existe
if [ -d "$diretorio1" ]; then
    echo "O diretório $diretorio já existe."
else
    # Cria o diretório se não existir
    mkdir -p "$diretorio1"
    echo "Diretório $diretorio criado com sucesso."
fi

# Verifica se o diretório existe
if [ -d "$diretorio2" ]; then
    echo "O diretório $diretorio já existe."
else
    # Cria o diretório se não existir
    mkdir -p "$diretorio2"
    echo "Diretório $diretorio criado com sucesso."
fi


# Path to the dataset directory
dataset_dir="/home/hcb/rawData/HCBagilent1/FASTQcat/teste"

# Iterate over each R1 file in the dataset directory
for r1_file in "${dataset_dir}"/*_R1.fastq.gz; do
    # Construct the corresponding R2 filename based on the R1 filename
    r2_file="${r1_file/_R1/_R2}"

    # Run Trimmomatic for the current pair
    trimmomatic PE -threads 6 -phred33 \
        "${r1_file}" "${r2_file}" \
        "${diretorio1}/$(basename ${r1_file})" \
        "${diretorio2}/$(basename ${r1_file})" \
        "${diretorio1}/$(basename ${r2_file})" \
        "${diretorio2}/$(basename ${r2_file})" \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
done

cd paired/

fastqc *

diretorio3="./fastqc"

# Verifica se o diretório existe
if [ -d "$diretorio3" ]; then
    echo "O diretório $diretorio já existe."
else
    # Cria o diretório se não existir
    mkdir -p "$diretorio3"
    echo "Diretório $diretorio criado com sucesso."
fi

mv *_fastqc* ./fastqc

## Os arquivos para trabalho estão na pasta paired
