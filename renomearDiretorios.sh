#!/bin/bash

## Oobjetivo desse script é renomear os diretórios extraídos do basespace
## Para cada corrida deve ser atualizado o directory_path
## No final da execução as pastas terão nome mais limpo

# Directory where the target directories are located
directory_path="/home/hcb/rawData/HCBagilent1/"

# List all directories in the specified path
directories=($(ls -d "$directory_path"/*/))

# Loop through each directory
for old_directory_path in "${directories[@]}"; do
    # Extract the directory name from the path
    old_directory_name=$(basename "$old_directory_path")

    # Extract the desired directory name using awk
    new_directory_name=$(echo "$old_directory_name" | awk -F '_' '{print $1}')

    # Check if the old directory exists
    if [ -d "$old_directory_path" ]; then
        # Rename the directory
        mv "$old_directory_path" "$directory_path/$new_directory_name"
        echo "Directory renamed successfully from $old_directory_name to $new_directory_name"
    else
        echo "Error: The directory $old_directory_name does not exist."
    fi
done




