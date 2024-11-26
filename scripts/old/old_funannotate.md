Script to run PASA only

````
#!/bin/bash

#SBATCH --job-name=fun_train2                  	# Job name
#SBATCH --cpus-per-task=64                     # Number of cores reserved
#SBATCH --mem-per-cpu=2G                       # Memory per core (total: 128GB)
#SBATCH --time=168:00:00                       # Max runtime
#SBATCH --qos=1week                            # Queue (time-based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/fun_train2_out.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/fun_train2_err.log

# Initialize conda for bash
eval "$(conda shell.bash hook)"

# Activate conda environment
conda activate funannotate

# Set working directory
cd /scicore/home/ebertd/dexter0000/aphid/annotation

# Set path variables
ASSEMBLY="/scicore/home/ebertd/dexter0000/aphid/annotation/output/repeatmasker_output/aphid.filtered.cleaned.renamed.fa.masked"
OUTDIR="/scicore/home/ebertd/dexter0000/aphid/annotation/output/funannotate_train_full"
TRINITY_FILE="/scicore/home/ebertd/dexter0000/aphid/annotation/output/funannotate_train_full/training/trinity.fasta.clean"
LEFT_NORM="/scicore/home/ebertd/dexter0000/aphid/annotation/output/funannotate_train_full/training/normalize/trimmed_left.fastq.gz.normalized_K25_maxC50_minC5_maxCV10000.fq"
RIGHT_NORM="/scicore/home/ebertd/dexter0000/aphid/annotation/output/funannotate_train_full/training/normalize/trimmed_right.fastq.gz.normalized_K25_maxC50_minC5_maxCV10000.fq"

# Define the number of threads to use
THREADS=64

# Create a custom temporary directory in your project folder
SCRATCH_DIR="/scicore/home/ebertd/dexter0000/aphid/annotation/tmp"
mkdir -p "$SCRATCH_DIR"
export TMPDIR="$SCRATCH_DIR"

# Run only PASA step in funannotate train
funannotate train -i "$ASSEMBLY" -o "$OUTDIR" \
    --trinity "$TRINITY_FILE" \
    --left_norm "$LEFT_NORM" \
    --right_norm "$RIGHT_NORM" \
    --species "Aphis fabae" \
    --stranded no \
    --cpus "$THREADS" \
    --no_trimmomatic \
    --no_normalize_reads
   

# Clean up the temporary directory after the job completes
rm -rf "$SCRATCH_DIR"
````



Script to just train using trinity

````
srun --nodes=1 --cpus-per-task=16 --mem=64G --pty bash

# Initialize conda for bash
eval "$(conda shell.bash hook)"

# Activate conda environment
conda activate funannotate

# Set file paths
ASSEMBLY="/scicore/home/ebertd/dexter0000/aphid/annotation/output/repeatmasker_output/aphid.filtered.cleaned.renamed.fa.masked"
OUTDIR="/scicore/home/ebertd/dexter0000/aphid/annotation/output/funannotate_train_short"
THREADS=16

# Set working directory
cd /scicore/home/ebertd/dexter0000/aphid/annotation

if [ -d "$OUTDIR" ]; then
    rm -r "$OUTDIR"
fi

funannotate train -i "$ASSEMBLY" -o "$OUTDIR" \
	--trinity RNAseq/transcripts/Transcripts_Aphis_fabae.fasta \
	--species "Aphis fabae" \
    --cpus "$THREADS" \
    --pasa_db sqlite
````



````
#!/bin/bash

#SBATCH --job-name=fun_predict                    # Job name
#SBATCH --cpus-per-task=32                        # Number of cores
#SBATCH --mem-per-cpu=2G                          # Memory per CPU (total 64GB)
#SBATCH --time=168:00:00                          # Max runtime
#SBATCH --qos=1week                               # Queue (time-based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/fun_pred_out.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/fun_pred_err.log

# Load environment with Funannotate
eval "$(conda shell.bash hook)"
conda activate funannotate

# Set paths to your data files
ASSEMBLY="/scicore/home/ebertd/dexter0000/aphid/annotation/output/repeatmasker_output/aphid.filtered.cleaned.renamed.fa.masked"
TRAIN_DIR="/scicore/home/ebertd/dexter0000/aphid/annotation/output/funannotate_train_short/training"
OUTDIR="/scicore/home/ebertd/dexter0000/aphid/annotation/output/funannotate_predict"
SCRATCH_DIR="/scicore/home/ebertd/dexter0000/aphid/annotation/scratch"

# Define the number of threads to use
THREADS=32

# Navigate to project folder
cd /scicore/home/ebertd/dexter0000/aphid/annotation

# Create scratch directory if it doesn't exist
if [ ! -d "$SCRATCH_DIR" ]; then
    mkdir -p "$SCRATCH_DIR"
fi

# Make sure the output directory doesn't already exist. Funannotate will throw an error if a previous run exists.
if [ -d "$OUTDIR" ]; then
    rm -r "$OUTDIR"
fi

# Run funannotate predict with updated parameters, using scratch directory as tmpdir
funannotate predict -i "$ASSEMBLY" \
    -o "$OUTDIR" \
    -s "Aphis fabae" \
    --cpus "$THREADS" \
    --augustus_species "insect" \
    --transcript_evidence "$TRAIN_DIR/trinity.fasta.clean"  \
    --genemark_mode ET \
    --other_gff "$TRAIN_DIR/pasa.step1.gff3:10" \
    --max_intronlen 3000 \
    --organism other \
    --repeats2evm \
    --busco_db insecta \
    --tmpdir "$SCRATCH_DIR"

# Optional: clean up scratch directory if desired after run completion
rm -r "$SCRATCH_DIR"
````

This skips the RNAseq transcripts



````
#!/bin/bash

#SBATCH --job-name=fun_predict                    # Job name
#SBATCH --cpus-per-task=32                        # Number of cores
#SBATCH --mem-per-cpu=2G                          # Memory per CPU (total 64GB)
#SBATCH --time=168:00:00                          # Max runtime
#SBATCH --qos=1week                               # Queue (time-based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/fun_pred_out.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/fun_pred_err.log

# Load environment with Funannotate
eval "$(conda shell.bash hook)"
conda activate funannotate

# Set paths to your data files
ASSEMBLY="/scicore/home/ebertd/dexter0000/aphid/annotation/output/repeatmasker_output/aphid.filtered.cleaned.renamed.fa.masked"
TRAIN_DIR="/scicore/home/ebertd/dexter0000/aphid/annotation/output/funannotate_train_short/training"
OUTDIR="/scicore/home/ebertd/dexter0000/aphid/annotation/output/funannotate_predict"

# Define the number of threads to use
THREADS=32

# Navigate to project folder
cd /scicore/home/ebertd/dexter0000/aphid/annotation

# Make sure the output directory doesn't already exist. Funnanotate train will throw an error if a previous run exists.
if [ -d "$OUTDIR" ]; then
    rm -r "$OUTDIR"
fi

# Run funannotate predict with updated parameters
funannotate predict -i "$ASSEMBLY" \
    -o "$OUTDIR" \
    -s "Aphis fabae" \
    --cpus "$THREADS" \
    --augustus_species "insect" \
    --transcript_evidence "$TRAIN_DIR/trinity.fasta.clean"  \
    --genemark_mode ET \
    --other_gff "$TRAIN_DIR/pasa.step1.gff3:10" \
    --max_intronlen 3000 \
    --organism other \
    --repeats2evm \
    --busco_db insecta
````

