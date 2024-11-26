

# Create reduced read set

````
#!/bin/bash

#SBATCH --job-name=ASSEMBLE2	             	# Job name
#SBATCH --cpus-per-task=8                       # Number of cores reserved
#SBATCH --mem-per-cpu=4G                        # Memory reserved per core
                                                # Total memory reserved: 32GB
#SBATCH --time=24:00:00                         # Maximum time the job will run
#SBATCH --qos=1day                              # The job queue (time based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/HDEF_out_%A_%a.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/HDEF_err_%A_%a.log
#SBATCH --array=1-14%14                         # Array job specifications

# Define the number of threads to use
THREADS=8

# Navigate to the project folder
cd /scicore/home/ebertd/dexter0000/aphid

# File index
INDEXFILE=scripts/fastq_index.txt

# Name index
INDEXNAME=scripts/sampleID.txt

# File i from index
SAMP=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXFILE)

# Name i from index
NAME=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXNAME)

# Load required module
module load SAMtools

# Assign input arguments to variables

fasta_file="genespace/genomes_original/$NAME.p_ctg.filtered.fa"

bam_file="bamsSelfmap/$NAME.bam"

fastq_file="$SAMP"

output_fastq="reads_filtered/$NAME.fq"

# Step 1: Extract contig names from the .fasta file
grep "^>" "$fasta_file" | sed 's/^>//' > "$NAME"_contig_list

# Step 2: Identify read names from the .bam file that map to the contigs
samtools view "$bam_file" | awk -v OFS='\t' 'NR==FNR {contigs[$1]; next} $3 in contigs {print $1}' "$NAME"_contig_list - | sort | uniq > "$NAME"_mapped_read_names

# Step 3: Extract reads from the .fastq file that match the read names
# Note: seqtk must be installed and available in your PATH

module purge
module load seqtk
seqtk subseq "$fastq_file" "$NAME"_mapped_read_names > "$output_fastq"

# Clean up intermediate files
#rm "$NAME"_mapped_read_names "$mapped_read_names"

echo "Filtering complete. Output written to $output_fastq."
````



# Assemble metagenomes

````
#!/bin/bash

#SBATCH --job-name=FLYE2			             # Job name
#SBATCH --cpus-per-task=8                       # Number of cores reserved
#SBATCH --mem-per-cpu=4G                        # Memory reserved per core
                                                # Total memory reserved: 32GB
#SBATCH --time=24:00:00                         # Maximum time the job will run
#SBATCH --qos=1day                              # The job queue (time based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/FLYE_out_%A_%a.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/FLYE_err_%A_%a.log
#SBATCH --array=1-14%14                         # Array job specifications

# Define the number of threads to use
THREADS=8

# Initialize conda for bash
eval "$(conda shell.bash hook)"

#Load conda environment
module load Flye

# Navigate to the project folder
cd /scicore/home/ebertd/dexter0000/aphid

# Name index
INDEXNAME=scripts/sampleID.txt

# Name i from index
NAME=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXNAME)

# Ensure the output directory exists
mkdir -p FLYE2/"$NAME"

# Assemble genome
flye --pacbio-hifi reads_filtered/"$NAME".fq --out-dir FLYE2/"$NAME" --threads "$THREADS"
````

# Backmap reads

````bash
#!/bin/bash

#SBATCH --job-name=SAMPLE_BACKMAP	            # Job name
#SBATCH --cpus-per-task=8                       # Number of cores reserved
#SBATCH --mem-per-cpu=4G                        # Memory reserved per core
                                                # Total memory reserved: 16GB
#SBATCH --time=24:00:00                         # Maximum time the job will run
#SBATCH --qos=1day                              # The job queue (time based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/backmap_out_%A_%a.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/backmap_err_%A_%a.log
#SBATCH --array=1-14%14                        # Array job specifications

# Define the number of threads to use
THREADS=8

# Load required modules
module load minimap2
module load SAMtools

# Navigate to the project folder
cd /scicore/home/ebertd/dexter0000/aphid

#File index
INDEXFILE=scripts/fastq_index.txt

#Name index
INDEXNAME=scripts/sampleID.txt

#File i from index
SAMP=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXFILE)

#Name i from index
NAME=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXNAME)

OUTPREFIX=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXFILE)

# Define the reference genome
REF="FLYE2/$NAME/assembly.fasta"

# Ensure output directory exists
mkdir -p FLYE2/bamsSelfmap

# Map reads to the reference genome
minimap2 -ax map-pb -t "$THREADS" "$REF" "$SAMP" | samtools view -bS -@ $THREADS - | samtools sort -@ $THREADS -o FLYE2/bamsSelfmap/"$NAME".bam

# Index the reference genome (must be .csi format index)
samtools index -c FLYE2/bamsSelfmap/"$NAME".bam
````

# Find busco

````
#!/bin/bash

#SBATCH --job-name=BUSCO			            # Job name
#SBATCH --cpus-per-task=8                       # Number of cores reserved
#SBATCH --mem-per-cpu=4G                        # Memory reserved per core
                                                # Total memory reserved: 16GB
#SBATCH --time=24:00:00                         # Maximum time the job will run
#SBATCH --qos=1day                              # The job queue (time based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/busco_out_%A_%a.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/busco_err_%A_%a.log
#SBATCH --array=1-14%14                        # Array job specifications

# Define the number of threads to use
THREADS=8

# Load required modules
module load BUSCO

# Navigate to the project folder
cd /scicore/home/ebertd/dexter0000/aphid

#Name index
INDEXNAME=scripts/sampleID.txt

#Name i from index
NAME=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXNAME)

# Define the reference genome
REF="/scicore/home/ebertd/dexter0000/aphid/FLYE2/$NAME/assembly.fasta"

# Move into output directory
mkdir -p FLYE2/busco
cd FLYE2/busco

# Calculate BUSCO scores
busco -i "$REF" \
-m genome \
-l /scicore/home/ebertd/dexter0000/aphid/blobtools/busco/enterobacterales_odb10 \
-c 8 \
-o $NAME \
-f --offline
````

# BLAST

````
#!/bin/bash

#SBATCH --job-name=FL_BLAST			            # Job name
#SBATCH --cpus-per-task=16                      # Number of cores reserved
#SBATCH --mem-per-cpu=8G                        # Memory reserved per core
                                                # Total memory reserved: 256GB
#SBATCH --time=24:00:00                        # Maximum time the job will run
#SBATCH --qos=1day                             # The job queue (time based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/blast_out_%A_%a.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/blast_err_%A_%a.log
#SBATCH --array=1-14%14                        # Array job specifications

# Define the number of threads to use
THREADS=16

# Load required modules
module load BLAST+

# Navigate to the project folder
cd /scicore/home/ebertd/dexter0000/aphid

#Name index
INDEXNAME=scripts/sampleID.txt

#Name i from index
NAME=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXNAME)

# Define the reference genome
REF="/scicore/home/ebertd/dexter0000/aphid/FLYE2/$NAME/assembly.fasta"

# Navigate to the results directory
cd blobtools

# Set up results directory
mkdir -p FLYE2_"$NAME"

# Perform BLAST against the NT database
blastn -db nt/nt \
       -query "$REF" \
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 10 \
       -max_hsps 1 \
       -evalue 1e-25 \
       -num_threads 16 \
       -out FLYE2_"$NAME"/"$NAME".ncbi.blastn.out
````

# Run blobtools

````
# Activate environment
conda activate btk

# Navigate to project folder
cd aphid

mkdir -p FLYE2/blobtools

while IFS= read -r SAMP; do
  echo "Processing sample ID: $SAMP"

  # Run blobtools create with the current sample ID
  blobtools create --replace \
    --fasta FLYE2/$SAMP/assembly.fasta \
    --cov FLYE2/bamsSelfmap/"$SAMP".bam \
    --hits blobtools/FLYE2_"$SAMP"/"$SAMP".ncbi.blastn.out \
    --taxdump blobtools/taxdump \
    --busco FLYE2/busco/"$SAMP"/run_enterobacterales_odb10/full_table.tsv \
    FLYE2/blobtools/"$SAMP"

done < scripts/sampleID.txt
````



