# Chromosome 5 synteny analysis

# Aggregate genomes



## Step 1: format_genomes.sh

This script performs various formatting operations on the input  genome assembly files required for the GENESPACE program. The input genomes need to be stored in the "genomes_original" directory and indexed using a file named "genome_index.txt". It only takes one second to run, so it can be performed interactively.

````
# Set up the folders
mkdir -p genespace genespace/genomes_original genespace/genomes_processed

# Navigate to directory
cd genespace

#Make the index file
ls genomes_original > genome_index.txt

# Set index as a variable
INDEXFILE="/scicore/home/ebertd/dexter0000/aphid/genespace/genome_index.txt"

# Loop through each line of the genome index file
while IFS= read -r GENOME; do
    # Step 1: Convert all nucleotides to uppercase
    # This AWK command reads the original genome file and converts all lowercase nucleotides to uppercase.
    # It skips header lines (those starting with '>') and outputs them unchanged.
    awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' genomes_original/"$GENOME" |
    
    # Step 2: Remove non-numeric characters from fasta headers and leading zeros
    # This sed command processes the output from the previous step.
    # It modifies fasta header lines (those starting with '>') to remove any non-numeric characters and leading zeros.
    sed -E '/^>/ { s/[^0-9>]//g; s/>0+/>/; }' |
    
    # Step 3: Ensure all contigs are numbered
    # This sed command ensures that contig headers are numbered properly.
    # If a header line becomes empty after the previous step, it adds a zero to it.
    sed -e '/^>$/s/^>$/>0/' > genomes_processed/"$GENOME"

    # Check if the final file was created successfully
    if [ -s genomes_processed/"$GENOME" ]; then 
        printf "Finished formatting %s\n" "$GENOME"
    else 
        printf "There was a problem formatting %s\n" "$GENOME"
    fi
done < "$INDEXFILE"
````



# Aggregate genomes

The next steps are easier if we place all the genomes into one folder

````
# Navigate to project folder
cd aphids

# Move all the filtered assemblies into one folder
find assemblies -type f -name "*filtered.fa" -exec cp {} assemblies_filtered/ \;

# Make an index of all the genomes in the folder
ls assemblies_filtered/ > genome_index.txt
````



# Prokka

This runs in less than 10 minutes

````
#!/bin/bash

#SBATCH --job-name=PROKKA                     # Job name
#SBATCH --cpus-per-task=8                     # Number of cores reserved
#SBATCH --mem-per-cpu=4G                      # Memory reserved per core
                                              # Total memory reserved: 32GB
#SBATCH --time=1:00:00                       # Maximum time the job will run
#SBATCH --qos=6hours                            # The job queue (time based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/prokka_out_%A_%a.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/prokka_err_%A_%a.log
#SBATCH --array=1-18%18                      # Array job specifications

# Define the number of threads to use
THREADS=8

# Load required modules
module load prokka

# Define variables
# Directory containing genome files
INPUT_DIR="/scicore/home/ebertd/dexter0000/aphid/assemblies_filtered"    

# Output directory for Prokka annotations
OUTPUT_DIR="/scicore/home/ebertd/dexter0000/aphid/prokka" 

# File index
INDEXFILE="/scicore/home/ebertd/dexter0000/aphid/genome_index.txt"

# Navigate to the input directory
cd "$INPUT_DIR"

# Retrieve the sample name and genome file for this task ID
SAMP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$INDEXFILE")

# Run Prokka on the genome file
prokka --outdir "$OUTPUT_DIR/$SAMP" --prefix "$SAMP" --cpus "$THREADS" "$INPUT_DIR/$SAMP" --addgenes --force

# Check if Prokka ran successfully
if [ $? -eq 0 ]; then
    echo "Prokka annotation for $SAMP completed successfully."
else
    echo "Prokka annotation for $SAMP failed."
fi
````

Need to convert GFF to bed format for Orthovenn3. Stupid script throws an error you can ignore, but you have to remove the _gene suffix from the entries in the gff file.

````
#Navigate to prokka directory
cd prokka

# Define the directory containing the Prokka output folders
PROKKA_OUTPUT_DIR="/scicore/home/ebertd/dexter0000/aphid/prokka"

# Define the path to the Python script
PYTHON_SCRIPT="/scicore/home/ebertd/dexter0000/aphid/scripts/gff_to_bed.py"

# Loop through each subdirectory in the Prokka output directory
for dir in "$PROKKA_OUTPUT_DIR"/*; do
    if [ -d "$dir" ]; then
        # Define the GFF and BED file paths
        GFF_FILE="$dir/$(basename "$dir").gff"
        BED_FILE="$dir/$(basename "$dir").bed.gff"
        
        # Check if the GFF file exists
        if [ -f "$GFF_FILE" ]; then
            echo "Processing $GFF_FILE..."
            
            # Run the Python script to convert GFF to BED
            python "$PYTHON_SCRIPT" "$GFF_FILE" "$BED_FILE"
            
            # Check if the BED file was created successfully
            if [ -f "$BED_FILE" ]; then
                # Use sed to remove the "_gene" suffix from the second field
                sed -i 's/\(.*\t.*\)_gene/\1/' "$BED_FILE"
                echo "Modified $BED_FILE successfully."
            else
                echo "Failed to create BED file for $GFF_FILE."
            fi
        else
            echo "GFF file not found in $dir."
        fi
    fi
done
````



# Move files for orthovenn3 into one folder

````
cd prokka

mkdir -p orthovenn

# Move all the filtered assemblies into one folder
find . -type f -name "*faa" -exec cp {} orthovenn/ \;
find . -type f -name "*.bed.gff" -exec cp {} orthovenn/ \;

````







