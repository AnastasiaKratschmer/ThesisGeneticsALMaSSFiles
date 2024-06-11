#!/bin/bash
#SBATCH --account almass_genes_ana
#SBATCH -c 8
#SBATCH --mem 32g
#SBATCH --time 4:00:00
# Directory containing the CSV files
directory="collected_filesBIGRUN"

# Define the keyword you're looking for in the file names
keyword="FST_output_dataset_style"

# Create a temporary directory for storing intermediate files
temp_dir=$(mktemp -d)

# Find files containing the keyword in their titles
find "$directory" -type f -name "*$keyword*" | while read -r file; do
    # Copy files to temporary directory
    cp "$file" "$temp_dir"
done

# Merge the files in the temporary directory
cat "$temp_dir"/* > FST_mergedBIGRUN1.txt

# Optionally, clean up the temporary directory
rm -r "$temp_dir"
####
# Define the keyword you're looking for in the file names
keyword="FIS_output_dataset_style"

# Create a temporary directory for storing intermediate files
temp_dir=$(mktemp -d)

# Find files containing the keyword in their titles
find "$directory" -type f -name "*$keyword*" | while read -r file; do
    # Copy files to temporary directory
    cp "$file" "$temp_dir"
done

# Merge the files in the temporary directory
cat "$temp_dir"/* > FIS_mergedBIGRUN1.txt

# Optionally, clean up the temporary directory
rm -r "$temp_dir"
###
# Define the keyword you're looking for in the file names
keyword="PopulationInQuadrants"

# Create a temporary directory for storing intermediate files
temp_dir=$(mktemp -d)

# Find files containing the keyword in their titles
find "$directory" -type f -name "*$keyword*" | while read -r file; do
    # Copy files to temporary directory
    cp "$file" "$temp_dir"
done

# Merge the files in the temporary directory
cat "$temp_dir"/* > PopInQuad_mergedBIGRUN1.txt

# Optionally, clean up the temporary directory
rm -r "$temp_dir"
###
# Define the keyword you're looking for in the file names
keyword="DistanceAndGeneticDifference"

# Create a temporary directory for storing intermediate files
temp_dir=$(mktemp -d)

# Find files containing the keyword in their titles
find "$directory" -type f -name "*$keyword*" | while read -r file; do
    # Copy files to temporary directory
    cp "$file" "$temp_dir"
done

# Merge the files in the temporary directory
cat "$temp_dir"/* > DifferenceSimil_mergedBIGRUN1.txt

# Optionally, clean up the temporary directory
rm -r "$temp_dir"
###
# Define the keyword you're looking for in the file names
keyword="Heterozyg_output"

# Create a temporary directory for storing intermediate files
temp_dir=$(mktemp -d)

# Find files containing the keyword in their titles
find "$directory" -type f -name "*$keyword*" | while read -r file; do
    # Copy files to temporary directory
    cp "$file" "$temp_dir"
done

# Merge the files in the temporary directory
cat "$temp_dir"/* > Heterozyg_mergedBIGRUN1.txt

# Optionally, clean up the temporary directory
rm -r "$temp_dir"
