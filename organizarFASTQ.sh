#!/bin/bash

## O objetivo desse script é concatenar os arquivos FASTQ de diversas linhas de sequenciamento
## Ao final da execução os haverá apenas um arquivo R1 e um arquivo R2 para cada amostra
## Os arquivos que serão utilizados para análise estão nas pasta FASTQcat

cd /home/hcb/rawData/HCBagilent1/
mkdir FASTQcat


new_directory="/home/hcb/rawData/HCBagilent1/FASTQcat"


for d in /home/hcb/rawData/HCBagilent1/*
do
    # Get the base name of the directory
    fileName=$(basename $d)

    # Check if the current directory is the new_directory
    if [ "$d" != "$new_directory" ]; then
        (cd $d && cat *L00{1,2,3,4}*R1*.fastq.gz > ${fileName}_R1.fastq.gz && \
        cat *L00{1,2,3,4}*R2*.fastq.gz > ${fileName}_R2.fastq.gz)

        # Move the recently created files to the new directory
        mv ${d}/${fileName}_R1.fastq.gz $new_directory
        mv ${d}/${fileName}_R2.fastq.gz $new_directory
    else
        echo "Skipping processing of new_directory: $d"
    fi
done

echo "Arquivos transferidos para a pasta FASTQcat com sucesso!"



