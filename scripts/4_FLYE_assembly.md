

# Create reduced read set

In the previous BLOBTOOLS module, we back-mapped the raw sequence reads to the filtered meta-genome assemblies. Now we use this mapping data to produce a filtered read set that will be used as input for a 2nd round of genome assembly using the FLYE assembler.

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



# Second round genome assembly using FLYE

FLYE produced less assembly errors with plasmids and circular elements than HIFIASM, but does not handle non-target contaminants as well as HIFISAM. Therefore HIFIASM is used for the first round of assembly to remove contaminants and FLYE is used for the second round to produce the final assemblies.

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



# Separate parameters for  S07 to complete

````
flye --pacbio-hifi reads_filtered/S07.fq --out-dir FLYE2/S07 --threads 16 --min-overlap 10000 --meta

flye --pacbio-hifi reads_filtered/S07.fq --out-dir FLYE2/temp --threads 16 --min-overlap 5000 --meta
# Notes on 7
4 appears to be part of - APSE CAN BE DELETED
3 appears to be part of - APSE CAN BE DELETED
6 (piece of APSE) is entirely contained inside 5 and mostly inside 4 - CAN BE DELETED
1 (Piece of APSE) is contained inside 5 and and mostly inside 3
````

