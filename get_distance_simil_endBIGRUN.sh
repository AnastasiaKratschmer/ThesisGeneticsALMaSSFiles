#!/bin/bash
#SBATCH --account almass_genes_ana
#SBATCH -c 8
#SBATCH --mem 32g
#SBATCH --time 4:00:00

# Check if a directory path is provided as an argument
if [ $# -ne 1 ]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

# Store the directory path provided as an argument
directory="$1"

# Define the word to search for in file titles
search_word="Distance"

# Loop through each file
for file in "${directory}"/*${search_word}*; do
    # Reverse the file and search for lines starting with "31"
    # Once a line starting with "30" is encountered, terminate the search
    tac "$file" | awk '/^30/ { exit } /^31/' | tac > "$file.tmp" && mv "$file.tmp" "$file"
done
