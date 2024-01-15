#!/bin/bash

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
genome_index="/home/hcb-info/RefSeq/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/"

# Path to the directory containing all the datasets
data_dir="/home/hcb-info/teste_STAR"

# Output directory
output_dir="/home/hcb-info/teste_STAR/STARresults"

# Iterate over each dataset in the data directory
# Iterate over each R1 file
for r1_file in "${data_dir}"/*_R1.fastq; do
    if [ -e "$r1_file" ]; then
        # Construct the corresponding R2 filename based on the R1 filename
        r2_file="${r1_file/_R1/_R2}"

        # Check if the corresponding R2 file exists
        if [ -e "$r2_file" ]; then
            # Extract dataset name from the R1 filename
            dataset_name=$(basename "$r1_file")

            # Run STAR aligner for the current pair
            /home/hcb-info/tools/STAR-Fusion-v1.12.0/STAR-Fusion \
                --genome_lib_dir $genome_index \
                --left_fq ${r1_file} \
                --right_fq ${r2_file} \
		--output_dir ${output_dir}/${dataset_name}
        else
            echo "Corresponding R2 file not found for $r1_file"
        fi
    fi
done

cd ${output_dir}

for file in *_R1.fastq*; do
    # Check if the file exists
    if [ -e "$file" ]; then
        # Extract the filename without the "_R1.fastq" string
        new_name=$(echo "$file" | sed 's/_R1.fastq//')

        # Rename the file
        mv "$file" "$new_name"
        echo "Renamed: $file to $new_name"
    fi
done
