# Aphid genome annotation pipeline

This pipeline performs gene prediction and annotation using the FUNANNOTATE pipeline tool. This pipeline contains many individual submodules (e.g. RepeatModeler, RepeatMasker, Augustus, InterproScan, etc.) and therefore is an all-in-one solution to eukaryotic genome annotation.

# Installation and set up

There is a lot of installation and configuration required to get the pipeline in place:

### Set up the new directory structure

````bash
# Navigate to project folder
cd aphid

# Make required directories
mkdir annotation
mkdir annotation/funannotate_db
mkdir annotation/inputProcessed
mkdir annotation/inputRaw
mkdir annotation/output
mkdir annotation/repeatDB
mkdir annotation/eggnog_db
````

### Install and configure tools

````bash
# Add the necessary Conda channels to ensure all dependencies are sourced correctly
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Create conda environment for the project
conda create -n funannotate "python>=3.6,<3.9" funannotate

# Install additional dependencies. THese are the easy ones. Some additional dependencies require more hands-on installation, and are detailed in their own code blocks below.
conda install -c conda-forge distro
conda install -c bioconda eggnog-mapper
conda install -c bioconda repeatmodeler
conda install -c bioconda repeatmasker
conda install -c bioconda sra-tools #If fetching RNAseq data from SRA repository

# Provide path to directory where databases will be stored
echo "export FUNANNOTATE_DB=/scicore/home/ebertd/dexter0000/aphid/annotation/funannotate_db" > /scicore/home/ebertd/dexter0000/miniconda3/envs/funannotate/etc/conda/activate.d/funannotate.sh

# Provide empty path for clearning databases (not really used)
echo "unset FUNANNOTATE_DB" > /scicore/home/ebertd/dexter0000/miniconda3/envs/funannotate/etc/conda/deactivate.d/funannotate.sh
    
# Activate and test the environment
conda activate funannotate
funannotate --help

#check that all modules are installed
funannotate check --show-versions

#Install the required databases
funannotate setup -i all

# Add the insecta busco database
funannotate setup -b insecta
````

### Install more complicated dependencies (Eggnog)

This one isn't available through conda or pip and has to be manually downloaded.

````bash
# Download the eggnog database file
download_eggnog_data.py --data_dir /scicore/home/ebertd/dexter0000/aphid/annotation/eggnog_db

# Tell funannotate where the eggnog database file is located
echo "export EGGNOG_DATA_DIR=/scicore/home/ebertd/dexter0000/aphid/annotation/eggnog_db" > /scicore/home/ebertd/dexter0000/miniconda3/envs/funannotate/etc/conda/activate.d/eggnog.sh

# Restart conda environment for changes to take effect
conda deactivate
conda activate funannotate

#check that eggnog is no longer on the missing dependancy list
funannotate check --show-versions
````

### Install more complicated dependencies (signalP)

This one is difficult to install because of licensing issues. The instructions are described here: https://anaconda.org/predector/signalp5

````
# Set up placeholder for signal installation
conda install predector::signalp5

# Download signalp-5.0b.Linux.tar.gz and move to directory where it is stored

# Complete the installation with the downloaded file
signalp5-register signalp-5.0b.Linux.tar.gz

# Confirm installation
funannotate check --show-versions

# Delete installation file that is no longer needed
rm signalp-5.0b.Linux.tar.gz
````

### Install more complicated dependencies (Genemark)

This is one is also difficult to install because of licensing issues. The instructions are described here: http://topaz.gatech.edu/GeneMark/license_download.cgi. Don't forget to separately extract and place the license key in the home directory (not the installation directory)!

````bash
# Decompress the downloaded installation files
tar -xzf gmes_linux_64_4.tar.gz

# Add the path to the conda environment
echo "export GENEMARK_PATH=/scicore/home/ebertd/dexter0000/aphid/annotation/gmes_linux_64_4" > /scicore/home/ebertd/dexter0000/miniconda3/envs/funannotate/etc/conda/activate.d/genemark.sh

# Restart conda environment for changes to take effect
conda deactivate
conda activate funannotate

# Confirm installation
funannotate check --show-versions

# The path still isn't always recognized, so further steps are needed
Steps to Add gmes_petap.pl to PATH in the Conda Environment
Locate the Conda Activation Directory: In your Funannotate environment, there should be an activate.d directory where you can add custom environment variables.

bash
Copy code
mkdir -p /scicore/home/ebertd/dexter0000/miniconda3/envs/funannotate/etc/conda/activate.d
Create a New Script to Add GeneMark to PATH: Create a new script file in this directory (let’s name it genemark_path.sh):

bash
Copy code
nano /scicore/home/ebertd/dexter0000/miniconda3/envs/funannotate/etc/conda/activate.d/genemark_path.sh
Add the PATH Export Command: In this file, add the following line to add gmes_petap.pl to the PATH every time the Funannotate environment is activated:

bash
Copy code
export PATH=$PATH:/scicore/home/ebertd/dexter0000/aphid/annotation/gmes_linux_64_4
Save and Exit (in nano, press CTRL + X, then Y, and Enter).

Deactivate and Reactivate the Environment: Deactivate and then reactivate the Funannotate environment to apply the change:

bash
Copy code
conda deactivate
conda activate funannotate
Verify the PATH: Check that gmes_petap.pl is now in your PATH by running:

bash
Copy code
which gmes_petap.pl
This approach will only apply the gmes_petap.pl path when the Funannotate environment is active, keeping your global environment clean. Let me know if this works as expected!
````



# Data pre-processing

The input files will require a lot of pre-processing before we can run the annotation pipeline

### Remove duplicate contigs

We want to identify and remove repetitive contigs that are contained in other scaffolds of the assembly, as those are likely assembly errors. If the repeats are indeed unique, then we want to keep them in the assembly. Funannotate has a module to clean up repetitive contigs  This is done using a “leave one out” approach using minimap2 or mummer (nucmer), where the the shortest contigs/scaffolds are aligned to the rest of the assembly to determine if it is repetitive. The script loops through the contigs starting with the shortest and workings its way to the N50 of the assembly, dropping contigs/scaffolds that are greater than the percent coverage of overlap (`--cov`) and the percent identity of overlap (`--pident`). This scripts sorts contigs by size, starting with shortest contigs it uses minimap2 to find contigs duplicated elsewhere, and then removes duplicated contigs. This step can be performed interactively because it only takes a few minutes to run.

```
# Request compute resources
srun --nodes=1 --cpus-per-task=16 --mem=64G --pty bash

# Activate conda environment
conda activate funannotate

# Navigate to the annotation directory
cd aphid/annotation

# Search for and clean out small duplicate contigs
funannotate clean --input inputRaw/aphid.filtered.fa --out inputProcessed/aphid.filtered.cleaned
```

## Sort/Rename FASTA Headers

Augustus also has problems with long contig/scaffold names. Funannotate provides a script to sort by size and rename the headers, we don't need to sort and rename them for this assembly. Instead we'll just remove the metadata from the fasta headers.

````
# Remove all metadata from contig headers
sed '/^>/ s/ .*//' inputProcessed/aphid.filtered.cleaned > inputProcessed/aphid.filtered.cleaned.renamed.fa

````

### Repeat Modeling/Masking

Annotation requires that repetitive regions are soft-masked (lower-case letters) in the genome assembly. Otherwise it will produce a lot of spurious gene predictions. The funannotate pipeline includes the option to pre-process the genome using repeatModeler and repeatMasker, but it doesn't use the new and improved versions of those tools, so we're going to call them separately. This step requires at least several days to one week of processing and so it submitted as a script on a cluster.

### RepeatMod.sh

Repeat modeling and repeat masking are closely related tasks, but we'll run them separately to make the pipeline modular. First we generate a de novo repeat library for the species. If we start with a good genome assembly, then we can reuse this repeat library for other genome assemblies for the same species.

````bash
#!/bin/bash

#SBATCH --job-name=repeatModel                 # Job name
#SBATCH --cpus-per-task=16                     # Number of cores reserved
#SBATCH --mem-per-cpu=4G                       # Memory per core (total: 64GB)
#SBATCH --time=168:00:00                       # Max runtime
#SBATCH --qos=1week                            # Queue (time-based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/repeatMod_out.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/repeatMod_err.log

# Define the number of threads to use
THREADS=16

# Define input/output paths
INFILE="/scicore/home/ebertd/dexter0000/aphid/annotation/inputProcessed/aphid.filtered.cleaned.renamed.fa"
OUTDIR="/scicore/home/ebertd/dexter0000/aphid/annotation/output/repeatmodeler_output"
REPEAT_LIB="$OUTDIR/consensi.fa.classified"

# Initialize conda for bash
eval "$(conda shell.bash hook)"

# Load conda environment
conda activate funannotate

# Confirm environment activation
if [[ "$CONDA_DEFAULT_ENV" != "funannotate" ]]; then
    echo "Conda environment failed to activate" >&2
    exit 1
fi

# Disable usage reporting in BLAST
export BLAST_USAGE_REPORT=false

# Set working directory
cd /scicore/home/ebertd/dexter0000/aphid/annotation

# Create the output directory if it does not exist
mkdir -p "$OUTDIR"

# Step 1: Build the database
BuildDatabase -name "$OUTDIR/aphid_repeat_db" "$INFILE"

# Step 2: Run RepeatModeler with LTRStruct, directing output to a fixed location
RepeatModeler -LTRStruct -database "$OUTDIR/aphid_repeat_db" -threads "$THREADS" -dir "$OUTDIR"

# Verify RepeatModeler output
if [[ -f "$REPEAT_LIB" ]]; then
    echo "RepeatModeler completed successfully. Repeat library saved at $REPEAT_LIB"
else
    echo "Error: RepeatModeler did not generate consensi.fa.classified" >&2
    exit 1
fi
````

### RepeatMask.sh

We can use the repeat model that we developed in the previous step to soft-mask (lowercase letters) repeat elements. This is mandatory for gene prediction and annotation, and it take a very long time to run (up to a week).

````bash
#!/bin/bash

#SBATCH --job-name=repeatMask                  # Job name
#SBATCH --cpus-per-task=32                     # Number of cores reserved
#SBATCH --mem-per-cpu=2G                       # Memory per core (total: 64GB)
#SBATCH --time=168:00:00                       # Max runtime
#SBATCH --qos=1week                            # Queue (time-based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/repeatMask_quick_out.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/repeatMask_quick_err.log

# Define the number of threads to use
THREADS=32

# Define input/output paths
INFILE="/scicore/home/ebertd/dexter0000/aphid/annotation/inputProcessed/aphid.filtered.cleaned.renamed.fa"
OUTDIR="/scicore/home/ebertd/dexter0000/aphid/annotation/output/repeatmasker_output_fast"
REPEAT_LIB="/scicore/home/ebertd/dexter0000/aphid/annotation/output/repeatmodeler_output/consensi.fa.classified"

# Initialize conda for bash
eval "$(conda shell.bash hook)"

# Load conda environment
conda activate funannotate

# Confirm environment activation
if [[ "$CONDA_DEFAULT_ENV" != "funannotate" ]]; then
    echo "Conda environment failed to activate" >&2
    exit 1
fi

# Disable usage reporting in BLAST
export BLAST_USAGE_REPORT=false

# Set working directory
cd /scicore/home/ebertd/dexter0000/aphid/annotation

# Create the output directory if it does not exist
mkdir -p "$OUTDIR"

# Run repeatMasker
RepeatMasker -qq -lib "$REPEAT_LIB" -pa "$THREADS" -gff -xsmall -dir "$OUTDIR" "$INFILE" -noisy
````



# Now run the freakin script

### Train2

````bash
#!/bin/bash

#SBATCH --job-name=fun_train                  	# Job name
#SBATCH --cpus-per-task=64                     # Number of cores reserved
#SBATCH --mem-per-cpu=2G                       # Memory per core (total: 128GB)
#SBATCH --time=168:00:00                       # Max runtime
#SBATCH --qos=1week                            # Queue (time-based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/fun_train_out.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/fun_train_err.log

# Initialize conda for bash
eval "$(conda shell.bash hook)"

# Activate conda environment
conda activate funannotate

# Set working directory
cd /scicore/home/ebertd/dexter0000/aphid/annotation

# Set path variables
ASSEMBLY="/scicore/home/ebertd/dexter0000/aphid/annotation/output/repeatmasker_output/aphid.filtered.cleaned.renamed.fa.masked"
OUTDIR="/scicore/home/ebertd/dexter0000/aphid/annotation/output/funannotate_train_full"
FASTQ_DIR="RNAseq/fastq"

# Define the number of threads to use
THREADS=64

# Create a custom temporary directory in your project folder
SCRATCH_DIR="/scicore/home/ebertd/dexter0000/aphid/annotation/tmp"
mkdir -p "$SCRATCH_DIR"
export TMPDIR="$SCRATCH_DIR"

# Make sure the output directory doesn't already exist. Funnanotate train will throw an error if a previous run exists.
if [ -d "$OUTDIR" ]; then
    rm -r "$OUTDIR"
fi

# Gather left and right reads as individual variables
LEFT_READS=( "$FASTQ_DIR"/*_1.renamed.fastq.gz )
RIGHT_READS=( "$FASTQ_DIR"/*_2.renamed.fastq.gz )

funannotate train -i "$ASSEMBLY" -o "$OUTDIR" \
    --left "${LEFT_READS[@]}" \
    --right "${RIGHT_READS[@]}" \
    --species "Aphis fabae" \
    --stranded no \
    --cpus "$THREADS"
    
# Clean up the temporary directory after the job completes
rm -rf "$SCRATCH_DIR"
````



# Predict2

Note that the output folder has to be the same as in the training step.

````
#!/bin/bash

#SBATCH --job-name=fun_predict2                    # Job name
#SBATCH --cpus-per-task=32                        # Number of cores
#SBATCH --mem-per-cpu=2G                          # Memory per CPU (total 64GB)
#SBATCH --time=168:00:00                          # Max runtime
#SBATCH --qos=1week                               # Queue (time-based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/fun_pred2_out.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/fun_pred2_err.log

# Load environment with Funannotate
eval "$(conda shell.bash hook)"
conda activate funannotate

# Set paths to your data files
ASSEMBLY="/scicore/home/ebertd/dexter0000/aphid/annotation/output/repeatmasker_output/aphid.filtered.cleaned.renamed.fa.masked"
OUTDIR="/scicore/home/ebertd/dexter0000/aphid/annotation/output/funannotate_train_full"
SCRATCH_DIR="/scicore/home/ebertd/dexter0000/aphid/annotation/scratch"

# Define the number of threads to use
THREADS=32

# Navigate to project folder
cd /scicore/home/ebertd/dexter0000/aphid/annotation

# Create scratch directory if it doesn't exist
if [ ! -d "$SCRATCH_DIR" ]; then
    mkdir -p "$SCRATCH_DIR"
fi

# Run funannotate predict with updated parameters, using scratch directory as tmpdir
funannotate predict -i "$ASSEMBLY" \
    -o "$OUTDIR" \
    -s "Aphis fabae" \
    --cpus "$THREADS" \
    --organism other \
    --repeats2evm \
    --busco_db insecta \
    --ploidy 2 \
    --tmpdir "$SCRATCH_DIR" \
    --optimize_augustus \
    --max_intronlen 50000
    
# Optional: clean up scratch directory if desired after run completion
rm -r "$SCRATCH_DIR"
````



















