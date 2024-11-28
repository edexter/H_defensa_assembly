# Prepare working directory and data

The sequence data is compressed inside a tarball object. In needs to be untarred and decompressed into the project folder data directory. This only requires a minutes of processing time and can be performed interactively.

```
# Request computing resources
srun --nodes=1 --cpus-per-task=8 --mem=32G --pty bash

#Make project directory and a few subfolders
mkdir -p aphid/scripts aphid data_round1 aphid/data_round2 aphid/results aphid/scratch aphid/assemblies aphid/logs

# Navidate to data directory
cd data_round2

#Decompress the data tarball from the sequencing center
tar -xvf p30633_o35228_multiplexed_1_A01.dmxData.tar.gz
```

# Sample QC

Perform a basic set of sample QC on the reads to ensure that sequencing was successful. This only requires 10-15 minutes of processing time and can be performed interactively.

````bash
#Request lots of resources for interactive jobs
srun --nodes=1 --cpus-per-task=8 --mem=32G --pty bash

# Navigate to the project folder
cd aphid

#Load JELLYFISH
module load Jellyfish

#Count Kmers with JELLYFISH
jellyfish count -C -m 21 -s 2000M -t 8 -o results/kmer_count.jf <(zcat data_round2/fastx/m64141e_240627_145720.hifi_reads*.gz)

#Make a summary histogram file from JELLYFISH output
jellyfish histo results/kmer_count.jf > results/kmer_count.histo

# Use the R script included with genomescope to plot the histogram data 

# Load required module
module purge
module load R

Rscript /scicore/home/ebertd/dexter0000/software/genomescope/genomescope.R results/kmer_count.histo 21 10000 results/genomescope

#Run a custom script to extract some stats

#Load required module
module purge
module load seqtk

#Make sure the script is executable (if it isn't already)
chmod u+x scripts/asmstats

#Run ASM script
zcat data_round2/fastx/m64141e_240627_145720.hifi_reads*.gz | scripts/asmstats > results/asm_read_stats

#Plot the read length distribution
python scripts/plot_fasta_length.py data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2177--bc2177.hifi_reads.fastq.gz t1.length.png
````



# Assemble genome

Assemble the genome using HIFASM and, convert the primary assembly to fasta, and calculate some metrics.

````bash
#!/bin/bash

#SBATCH --job-name=HIFIASM_APHID				#Job name
#SBATCH --cpus-per-task=8	                  	#Number of cores reserved
#SBATCH --mem-per-cpu=32G              			#Memory reserved per core.
												#Total memory reserved: 256GB
#SBATCH --time=168:00:00	        			#Maximum time the job will run
#SBATCH --qos=1week          					#The job queue (time based)

#This is the stdout file
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/HIFASM_out

#This is the stderr file
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/HIFASM_err

#Things remember when running on sciCore:
#This job runs from the current working directory
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.

################################################################################
# Load modules and set variables
################################################################################

# Initialize conda for bash
eval "$(conda shell.bash hook)"

#Load conda environment
conda activate eukaryotic_genome_assembly

# Prefix for assembly output files
ASSEMBLY="APHID_HIFIASM"

# Directory containing the reads
READS_DIR="/scicore/home/ebertd/dexter0000/aphid/data_round2/fastx"

# Wildcard pattern to match all .fasta.gz files
READS_FILES="$READS_DIR/m64141e_240627_145720.hifi_reads*.gz"

################################################################################
# Assemble genome
################################################################################

# Navigate to project directory
cd /scicore/home/ebertd/dexter0000/aphid

# Create output directory (if it doesn't exist)
mkdir -p assemblies/HIFASM_round1

# This runs amazingly fast (~1 day)
hifiasm -f0 --primary -t 8 -o assemblies/HIFASM_round1/"$ASSEMBLY" $READS_FILES

#Convert HIFIASM output to fasta format
scripts/gfa2fa assemblies/HIFASM_round1/"$ASSEMBLY".p_ctg.gfa > assemblies/HIFASM_round1/"$ASSEMBLY".p_ctg.fa

#Get some assembly stats
scripts/asmstats assemblies/HIFASM_round1/"$ASSEMBLY".p_ctg.fa > assemblies/HIFASM_round1/"$ASSEMBLY".p_ctg.fa.stats
````



# Separate QC passes

First create a list of files

```
ls /scicore/home/ebertd/dexter0000/aphid/data_round2/fastx/m64141e_240627_145720.hifi_reads*.gz > scripts/fastq_index.txt
```

Calculate QC metrics over list of files

````
#!/bin/bash

#SBATCH --job-name=APHID_KMER_ANALYSIS          # Job name
#SBATCH --cpus-per-task=8                       # Number of cores reserved
#SBATCH --mem-per-cpu=4G                        # Memory reserved per core
                                                # Total memory reserved: 16GB
#SBATCH --time=24:00:00                         # Maximum time the job will run
#SBATCH --qos=1day                              # The job queue (time based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/kmer_analysis_out_%A_%a.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/kmer_analysis_err_%A_%a.log
#SBATCH --array=1-14%14                        # Array job specifications

# Initialize conda for bash (if needed)
eval "$(conda shell.bash hook)"

#Load conda environment
conda activate eukaryotic_genome_assembly

# Navigate to the project folder
cd /scicore/home/ebertd/dexter0000/aphid

# Create results directory if it doesn't exist
mkdir -p results

# Get the file to process for this task
READ_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scripts/fastq_index.txt)
BASENAME=$(basename $READ_FILE.gz)

# Count Kmers with JELLYFISH
module purge
module load Jellyfish
jellyfish count -C -m 21 -s 2000M -t 8 -o results/${BASENAME}_kmer_count.jf <(zcat $READ_FILE)

# Make a summary histogram file from JELLYFISH output
jellyfish histo results/${BASENAME}_kmer_count.jf > results/${BASENAME}_kmer_count.histo

# Use the R script included with genomescope to plot the histogram data
module purge
module load R
Rscript /scicore/home/ebertd/dexter0000/software/genomescope/genomescope.R results/${BASENAME}_kmer_count.histo 21 10000 results/${BASENAME}_genomescope

# Run ASM script
module purge
module load seqtk
zcat $READ_FILE | scripts/asmstats > results/${BASENAME}_asm_read_stats

# Plot the read length distribution
python scripts/plot_fasta_length.py $READ_FILE results/${BASENAME}_length.png

````

Merge results together from each sample

````
RESULTS_DIR="results"
OUTPUT_FILE="results/merged_asm_read_stats.csv"

# Find all files ending in "asm_read_stats" and save them to a temporary file
find "$RESULTS_DIR" -name "*reads.fastq.gz.gz_asm_read_stats" -print > temp_results.txt

# Add a header to the output file
echo "Filename,sum,n,mean,largest,smallest,N50,L50" > "$OUTPUT_FILE"

# Process each file to extract the required fields and append to the output file
while read -r file; do
  sum=$(awk -F '[=,]' '/sum =/ {gsub(/ /, "", $2); print $2}' "$file")
  n=$(awk -F '[=,]' '/sum =/ {gsub(/ /, "", $4); print $4}' "$file")
  mean=$(awk -F '[=,]' '/sum =/ {gsub(/ /, "", $6); print $6}' "$file")
  largest=$(awk -F '[=,]' '/sum =/ {gsub(/ /, "", $8); print $8}' "$file")
  smallest=$(awk -F '[=,]' '/sum =/ {gsub(/ /, "", $10); print $10}' "$file")
  N50=$(awk -F '[=,]' '/N50 =/ {gsub(/ /, "", $2); print $2}' "$file")
  L50=$(awk -F '[=,]' '/N50 =/ {gsub(/ /, "", $4); print $4}' "$file")
  
  echo "$file,$sum,$n,$mean,$largest,$smallest,$N50,$L50" >> "$OUTPUT_FILE"
done < temp_results.txt

# Clean up temporary file
rm temp_results.txt

echo "File processing complete. Output saved to $OUTPUT_FILE"
````



````
BASE_DIR="results"
OUTPUT_FILE="results/aggregated_genoscope.csv"

# Add header to the output file
echo "Sample_ID,Heterozygosity_Min,Heterozygosity_Max,Genome_Haploid_Length_Min,Genome_Haploid_Length_Max,Genome_Repeat_Length_Min,Genome_Repeat_Length_Max,Genome_Unique_Length_Min,Genome_Unique_Length_Max,Model_Fit_Min,Model_Fit_Max,Read_Error_Rate" > "$OUTPUT_FILE"

# Function to extract values from summary.txt
extract_values() {
  local file=$1
  local sample_id=$(basename $(dirname "$file"))

  local heterozygosity_min=$(awk '/Heterozygosity/ {print $2}' "$file")
  local heterozygosity_max=$(awk '/Heterozygosity/ {print $3}' "$file")
  local genome_haploid_length_min=$(awk '/Genome Haploid Length/ {print $4}' "$file")
  local genome_haploid_length_max=$(awk '/Genome Haploid Length/ {print $5}' "$file")
  local genome_repeat_length_min=$(awk '/Genome Repeat Length/ {print $4}' "$file")
  local genome_repeat_length_max=$(awk '/Genome Repeat Length/ {print $5}' "$file")
  local genome_unique_length_min=$(awk '/Genome Unique Length/ {print $4}' "$file")
  local genome_unique_length_max=$(awk '/Genome Unique Length/ {print $5}' "$file")
  local model_fit_min=$(awk '/Model Fit/ {print $3}' "$file")
  local model_fit_max=$(awk '/Model Fit/ {print $4}' "$file")
  local read_error_rate=$(awk '/Read Error Rate/ {print $3}' "$file")

  echo "$sample_id,$heterozygosity_min,$heterozygosity_max,$genome_haploid_length_min,$genome_haploid_length_max,$genome_repeat_length_min,$genome_repeat_length_max,$genome_unique_length_min,$genome_unique_length_max,$model_fit_min,$model_fit_max,$read_error_rate" >> "$OUTPUT_FILE"
}

# Find all summary.txt files and process them
find "$BASE_DIR" -type f -name "summary.txt" | while read file; do
  extract_values "$file"
done

echo "Aggregation complete. Output saved to $OUTPUT_FILE"
````



### Look for HDEF contigs in genome assembly

````
# Request resources
srun --nodes=1 --cpus-per-task=8 --mem=32G --pty bash

# Load required module
module load BLAST

# Build a blastDB from assembly
makeblastdb -in APHID_HIFIASM.p_ctg.fa -dbtype nucl -out blast/whole_db

BLASTQUERY=/scicore/home/ebertd/dexter0000/aphid/reference_genomes/GCA_000021705.1_ASM2170v1_genomic.fna

blastn -query "$BLASTQUERY" -db blast/whole_db -out blast/HDEF_blast_result.txt
````

