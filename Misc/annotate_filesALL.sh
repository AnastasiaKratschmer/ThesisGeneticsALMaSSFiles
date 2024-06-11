#!/bin/bash
#SBATCH --account almass_genes_ana
#SBATCH -c 8
#SBATCH --mem 32g
#SBATCH --time 10:00:00
# Directory containing the CSV files
directory="collected_filesBIGRUN"

# Iterate through each CSV file
for file in "$directory"/*.txt; do
    # Extract filename without extension
    filename=$(basename "$file" .txt)

    # Extract relevant parts from filename
    parts=($(echo "$filename" | awk -F'[_.]' '{print $2 "_" $3, $4 "_" $5}'))

    # Add "landscape, rep" to the very beginning of the first line
    sed -i '1s/^/landscape, rep, /' "$file"

    # Iterate through each line of the CSV file starting from the second line
    awk 'NR>1' "$file" | while IFS= read -r line; do
        # Append filename parts to the end of each line
        echo "$line,${parts[0]},${parts[1]}"
    done > "$file.tmp"  # Write output to a temporary file

    # Replace original file with modified one
    mv "$file.tmp" "$file"
done

