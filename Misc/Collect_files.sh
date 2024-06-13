#!/bin/bash
#SBATCH --account almass_genes_ana
#SBATCH -c 8
#SBATCH --mem 32g
#SBATCH --time 00:10:00

# Source directory containing subfolders
source_dir="output"
# Destination directory
destination_dir="collected_filesBIGRUN"

# File name to search for
file_name="PopulationInQuadrants.txt"

# Iterate over each file found in subdirectories
find "$source_dir" -type f -name "$file_name" | while read -r file_path; do
    # Extract the parent directory names
    inner_parent_dir=$(dirname "$file_path")
    inner_parent_dir_name=$(basename "$inner_parent_dir")
    outer_parent_dir=$(dirname "$inner_parent_dir")
    outer_parent_dir_name=$(basename "$outer_parent_dir")
    echo "Found file: $file_path"
    head -n1 "$file_path"

    # New file name with addon
    new_file_name="${outer_parent_dir_name}_${inner_parent_dir_name}_${file_name}" # Use both inner and outer folder names

    # Move the file to the destination directory with the new name
    cp "$file_path" "$destination_dir/$new_file_name"
done

# File name to search for
file_name="FST_output_dataset_style.txt"

# Iterate over each file found in subdirectories
find "$source_dir" -type f -name "$file_name" | while read -r file_path; do
    # Extract the parent directory names
    inner_parent_dir=$(dirname "$file_path")
    inner_parent_dir_name=$(basename "$inner_parent_dir")
    outer_parent_dir=$(dirname "$inner_parent_dir")
    outer_parent_dir_name=$(basename "$outer_parent_dir")
    echo "Found file: $file_path"
    head -n1 "$file_path"

    # New file name with addon
    new_file_name="${outer_parent_dir_name}_${inner_parent_dir_name}_${file_name}" # Use both inner and outer folder names

    # Move the file to the destination directory with the new name
    cp "$file_path" "$destination_dir/$new_file_name"
done

# File name to search for
file_name="FIS_output_dataset_style.txt"

# Iterate over each file found in subdirectories
find "$source_dir" -type f -name "$file_name" | while read -r file_path; do
    # Extract the parent directory names
    inner_parent_dir=$(dirname "$file_path")
    inner_parent_dir_name=$(basename "$inner_parent_dir")
    outer_parent_dir=$(dirname "$inner_parent_dir")
    outer_parent_dir_name=$(basename "$outer_parent_dir")
    echo "Found file: $file_path"
    head -n1 "$file_path"

    # New file name with addon
    new_file_name="${outer_parent_dir_name}_${inner_parent_dir_name}_${file_name}" # Use both inner and outer folder names

    # Move the file to the destination directory with the new name
    cp "$file_path" "$destination_dir/$new_file_name"
done

# File name to search for
file_name="DistanceAndGeneticDifference.txt"

# Iterate over each file found in subdirectories
find "$source_dir" -type f -name "$file_name" | while read -r file_path; do
    # Extract the parent directory names
    inner_parent_dir=$(dirname "$file_path")
    inner_parent_dir_name=$(basename "$inner_parent_dir")
    outer_parent_dir=$(dirname "$inner_parent_dir")
    outer_parent_dir_name=$(basename "$outer_parent_dir")
    echo "Found file: $file_path"
    head -n1 "$file_path"

    # New file name with addon
    new_file_name="${outer_parent_dir_name}_${inner_parent_dir_name}_${file_name}" # Use both inner and outer folder names

    # Move the file to the destination directory with the new name
    cp "$file_path" "$destination_dir/$new_file_name"
done

# File name to search for
file_name="Heterozyg_output.txt"

# Iterate over each file found in subdirectories
find "$source_dir" -type f -name "$file_name" | while read -r file_path; do
    # Extract the parent directory names
    inner_parent_dir=$(dirname "$file_path")
    inner_parent_dir_name=$(basename "$inner_parent_dir")
    outer_parent_dir=$(dirname "$inner_parent_dir")
    outer_parent_dir_name=$(basename "$outer_parent_dir")
    echo "Found file: $file_path"
    head -n1 "$file_path"

    # New file name with addon
    new_file_name="${outer_parent_dir_name}_${inner_parent_dir_name}_${file_name}" # Use both inner and outer folder names

    # Move the file to the destination directory with the new name
    cp "$file_path" "$destination_dir/$new_file_name"
done
