# Meta-genome assembly

This script assembles a meta-genome from all of the reads in the sample, and calculate a few summary statistics about the assemblies. There are many species represented here, so this is just a rough first step that will require filtering of individuals contigs based on BLOBTOOLS results. This script requires an index file of the sequencer file names (fastq_index.txt) and and an index file of the sample IDs (sampleID.txt).

````bash
#!/bin/bash

#SBATCH --job-name=HDEF_assemble             	# Job name
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

# Initialize conda for bash
eval "$(conda shell.bash hook)"

#Load conda environment
conda activate eukaryotic_genome_assembly

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

# Ensure the output directory exists
mkdir -p assemblies/"$NAME"

# Assemble genome
hifiasm --primary -t "$THREADS" -o assemblies/"$NAME"/"$NAME" $SAMP

#Convert HIFIASM output to fasta format
scripts/gfa2fa assemblies/"$NAME"/"$NAME".p_ctg.gfa > assemblies/"$NAME"/"$NAME".p_ctg.fa

#Get some assembly stats
scripts/asmstats assemblies/"$NAME"/"$NAME".p_ctg.fa > assemblies/"$NAME"/"$NAME".p_ctg.fa.stats
````

