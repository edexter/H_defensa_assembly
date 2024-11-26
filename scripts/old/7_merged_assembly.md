# Assemble

````
#!/bin/bash

#SBATCH --job-name=HDEF_assemble             	# Job name
#SBATCH --cpus-per-task=16                       # Number of cores reserved
#SBATCH --mem-per-cpu=8G                        # Memory reserved per core
                                                # Total memory reserved: 32GB
#SBATCH --time=168:00:00                         # Maximum time the job will run
#SBATCH --qos=1week                              # The job queue (time based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/meta_out_%A_%a.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/meta_err_%A_%a.log

# Define the number of threads to use
THREADS=16

# Initialize conda for bash
eval "$(conda shell.bash hook)"

#Load conda environment
conda activate eukaryotic_genome_assembly

# Navigate to the project folder
cd /scicore/home/ebertd/dexter0000/aphid

# Ensure the output directory exists
mkdir -p assembly_merged

# Assemble genome
hifiasm --primary -t "$THREADS" -o assembly_merged/merged \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2177--bc2177.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2178--bc2178.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2179--bc2179.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2180--bc2180.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2181--bc2181.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2182--bc2182.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2183--bc2183.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2184--bc2184.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2185--bc2185.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2186--bc2186.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2187--bc2187.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2188--bc2188.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2189--bc2189.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2190--bc2190.hifi_reads.fastq.gz

#Convert HIFIASM output to fasta format
scripts/gfa2fa assembly_merged/merged.p_ctg.gfa > assembly_merged/merged.p_ctg.fa

#Get some assembly stats
scripts/asmstats assembly_merged/merged.p_ctg.fa > assembly_merged/merged.p_ctg.fa.stats
````

# Backmap reads to the assembled genomes

Blobtools needs read coverage for each genome assembly, so we align the raw reads from each sample against each metagenome assembly. This runs in less than 10 minutes.

````
#!/bin/bash

#SBATCH --job-name=SAMPLE_BACKMAP	            # Job name
#SBATCH --cpus-per-task=16                       # Number of cores reserved
#SBATCH --mem-per-cpu=8G                        # Memory reserved per core
                                                # Total memory reserved: 16GB
#SBATCH --time=24:00:00                         # Maximum time the job will run
#SBATCH --qos=1day                              # The job queue (time based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/backmap2_out_%A_%a.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/backmap2_err_%A_%a.log

# Define the number of threads to use
THREADS=16

# Load required modules
module load minimap2
module load SAMtools

# Navigate to the project folder
cd /scicore/home/ebertd/dexter0000/aphid

# Define the reference genome
REF=assembly_merged/merged.p_ctg.fa

# Map reads to the reference genome
minimap2 -ax map-hifi -t "$THREADS" "$REF" data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2177--bc2177.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2178--bc2178.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2179--bc2179.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2180--bc2180.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2181--bc2181.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2182--bc2182.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2183--bc2183.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2184--bc2184.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2185--bc2185.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2186--bc2186.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2187--bc2187.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2188--bc2188.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2189--bc2189.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2190--bc2190.hifi_reads.fastq.gz \
| samtools view -bS -@ "$THREADS" - | samtools sort -@ "$THREADS" -o assembly_merged/backmap.bam

# Index the reference genome (must be .csi format index)
samtools index -c assembly_merged/backmap.bam
````

# BUSCO

````
# Load required module
module load BUSCO

# Navigate to directory
cd assembly_merged

# Calculate BUSCO scores
busco -i merged.p_ctg.fa \
-m genome \
-l /scicore/home/ebertd/dexter0000/aphid/blobtools/busco/hemiptera_odb10/ \
-c 16 \
-o merged.busco.out \
-f --offline


````

# BLAST

````
#!/bin/bash

#SBATCH --job-name=BLAST			            # Job name
#SBATCH --cpus-per-task=16                      # Number of cores reserved
#SBATCH --mem-per-cpu=8G                        # Memory reserved per core
                                                # Total memory reserved: 256GB
#SBATCH --time=168:00:00                        # Maximum time the job will run
#SBATCH --qos=1week                             # The job queue (time based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/blast_out_merge.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/blast_err_merge.log

# Define the number of threads to use
THREADS=16

# Load required modules
module load BLAST+

# Navigate to the project folder
cd /scicore/home/ebertd/dexter0000/aphid

#Name
NAME=merged

# Define the reference genome
REF="/scicore/home/ebertd/dexter0000/aphid/assembly_merged/merged.p_ctg.fa"

# Navigate to the results directory
cd blobtools

# Set up results directory
mkdir -p "$NAME"

# Perform BLAST against the NT database
blastn -db nt/nt \
       -query "$REF" \
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 10 \
       -max_hsps 1 \
       -evalue 1e-25 \
       -num_threads 16 \
       -out "$NAME"/"$NAME".ncbi.blastn.out
````

# Populate the blobtools directory

````
# Activate environment
conda activate /scicore/home/ebertd/dexter0000/anaconda3/envs/btk

# Navigate to project folder
cd aphid

# Run blobtools create with the current sample ID
blobtools create --replace \
	--fasta assembly_merged/merged.p_ctg.fa \
    --cov assembly_merged/backmap.bam \
    --hits blobtools/merged/merged.ncbi.blastn.out \
    --taxdump blobtools/taxdump \
    --busco assembly_merged/merged.busco.out/run/full_table.tsv \
    viewer/merged
````



# Interactively filter the assembly

````
blobtools view --remote .
````

Open local viewer in another terminal

````
ssh -L 8001:127.0.0.1:8001 -L 8000:127.0.0.1:8000 dexter0000@login.scicore.unibas.ch
````

navigate to page in a web browser and use the interactive filtering tool. Filter as needed and save the filter parameters using the lists Menu

````
http://localhost:8001/view/all
````



# Create filtered assemblies

````
# Request resources for interactive node
srun --nodes=1 --cpus-per-task=8 --mem=32G --pty bash

# Activate environment
conda activate btk

# Navigate to project folder
cd aphid

# Run blobtools filter to extract Buchnera contigs
blobtools filter \
     --json blobtools/JSON/buchnera_blobtools_filter.json \
     --fasta assembly_merged/merged.p_ctg.fa \
     --output blobtools/filtered/buchnera \
     viewer/merged

# Run blobtools filter to extract Aphid contigs
blobtools filter \
     --json blobtools/JSON/aphid_blobtools_filter.json \
     --fasta assembly_merged/merged.p_ctg.fa \
     --output blobtools/filtered/aphid \
     viewer/merged

````

# Calculate ASM stats on filtered assemblies

````
# Activate environment
conda activate eukaryotic_genome_assembly

# Navigate to project folder
cd aphids

# Calculate ASM stats across all filtered assemblies
while IFS= read -r SAMP; do
  echo "Processing sample ID: $SAMP"

  # Run blobtools filter with the current sample ID
  scripts/asmstats curated/"$SAMP".curated.fasta > curated/"$SAMP".stats

done < scripts/sampleID.txt
````

# Recalculate BUSCO for Buchnera

````

````

# Multi-thread blast

````
/scicore/home/ebertd/dexter0000/aphid/scripts/fasta-splitter.pl --part-size 500000 --measure seq /scicore/home/ebertd/dexter0000/aphid/assembly_merged/merged.p_ctg.fa

ls merged.p_ctg.part* > partIndex.txt
````

````
#!/bin/bash

#SBATCH --job-name=BLASTmulti			            # Job name
#SBATCH --cpus-per-task=1                       # Number of cores reserved
#SBATCH --mem-per-cpu=8G                        # Memory reserved per core
                                                # Total memory reserved: 256GB
#SBATCH --time=168:00:00                        # Maximum time the job will run
#SBATCH --qos=1week                             # The job queue (time based)
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
REF="/scicore/home/ebertd/dexter0000/aphid/assemblies/"$NAME"/"$NAME".p_ctg.fa"

# Navigate to the results directory
cd blobtools

# Set up results directory
mkdir -p "$NAME"

# Perform BLAST against the NT database
blastn -db nt/nt \
       -query "$REF" \
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 10 \
       -max_hsps 1 \
       -evalue 1e-25 \
       -num_threads 16 \
       -out "$NAME"/"$NAME".ncbi.blastn.out

````

# Append species names to assembly

````

sed '/^>/ s/$/ [organism=Aphae fabae [strain=407]/' aphid.filtered.fa > aphid.filtered_with_info.fa

````

