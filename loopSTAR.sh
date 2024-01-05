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
            dataset_name=$(basename "$r1_file")

            # Run STAR aligner for the current pair
            STAR \
                --runMode alignReads \
		--runThreadN 16 \
                --outSAMtype BAM SortedByCoordinate \
                --genomeDir ${genome_index} \
                --readFilesIn ${r1_file} ${r2_file} \
		--outFileNamePrefix ${output_dir}/${dataset_name}_ \
		--outReadsUnmapped None \
		--twopassMode Basic \
		--outSAMstrandField intronMotif \
		--outSAMunmapped Within \
		--chimSegmentMin 12 \
          	--quantMode GeneCounts \
		--chimOutJunctionFormat 1 \
          	--alignMatesGapMax 100000 \
         	--alignSJstitchMismatchNmax 5 -1 5 5 \
          	--chimMultimapScoreRange 3 \
          	--chimScoreJunctionNonGTAG -4 \
          	--chimMultimapNmax 20 \
          	--chimNonchimScoreDropMin 10 \
          	--peOverlapNbasesMin 12 \
          	--peOverlapMMp 0.1 \
          	--alignInsertionFlush Right \
          	--alignSplicedMateMapLminOverLmate 0 \
          	--alignSplicedMateMapLmin 30
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

# Começando alterações 18/12

diretorio_f="./Fusions"

# Verifica se o diretório existe
if [ -d "$diretorio_f" ]; then
    echo "O diretório $diretorio já existe."
else
    # Cria o diretório se não existir
    mkdir -p "$diretorio_f"
    echo "Diretório $diretorio criado com sucesso."
fi

# Iterate over *Chimeric.out.junction files in the input directory
for junction_file in "${output_dir}"/*Chimeric.out.junction; do
    if [ -e "$junction_file" ]; then
        # Extract sample name from the file
        sample_name=$(basename "$junction_file" "_Chimeric.out.junction")

        # Run STAR-Fusion for the current file
        /home/hcb/tools/STAR-Fusion-v1.12.0/STAR-Fusion \
 --genome_lib_dir /home/hcb/tools/STAR-Fusion-v1.12.0/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ \
 -J ${junction_file} --output_dir "${diretorio_f}/${sample_name}_STAR-Fusion_output/"
    else
        echo "Chimeric.out.junction file not found for $junction_file"
    fi
done
