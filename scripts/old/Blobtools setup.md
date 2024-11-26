# Set up environment and install blobtools

Create and activate the conda environment

````bash
conda create -y -n btk -c conda-forge python=3.9
conda activate btk
````

Install blobtoolkit

````
pip install "blobtoolkit[full]"
````

Install personal copy of BLAST+ because the server installation is broken

````
# Download BLAST+
wget ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz

# Decompress BLAST+
tar -zxvf ncbi-blast-2.16.0+-x64-linux.tar.gz
````



# Prepare NCBI databases

First install the NT Database. This uses FTP so sciCore won't allow the transfer on an interactive node. Need to use to use the transfer node.

````
# Need to load Perl separately
module load Perl

# Need to use a local instal of Blast+ because the cluster install is broken.
/scicore/home/ebertd/dexter0000/aphid/blobtools/blast/ncbi-blast-2.16.0+/bin/update_blastdb.pl --decompress nt [*]
````

Next install and prepare the UniProt database.

````
mkdir -p uniprot
wget -q -O uniprot/reference_proteomes.tar.gz \
 ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/$(curl \
     -vs ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/ 2>&1 | \
     awk '/tar.gz/ {print $9}')
cd uniprot
tar xf reference_proteomes.tar.gz

touch reference_proteomes.fasta.gz
find . -mindepth 2 | grep "fasta.gz" | grep -v 'DNA' | grep -v 'additional' | xargs cat >> reference_proteomes.fasta.gz

echo -e "accession\taccession.version\ttaxid\tgi" > reference_proteomes.taxid_map
zcat */*/*.idmapping.gz | grep "NCBI_TaxID" | awk '{print $1 "\t" $1 "\t" $3 "\t" 0}' >> reference_proteomes.taxid_map

diamond makedb -p 16 --in reference_proteomes.fasta.gz --taxonmap reference_proteomes.taxid_map --taxonnodes ../taxdump/nodes.dmp -d reference_proteomes.dmnd
cd -
````

Finally, fetch the NCBI Taxdump.

````

mkdir -p taxdump;
cd taxdump;
curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -;
cd ..;
````



# Add BUSCO databases

````
mkdir -p busco

# Most general aphid database
wget -q -O eukaryota_odb10.gz "https://busco-data.ezlab.org/v4/data/lineages/eukaryota_odb10.2020-09-10.tar.gz" \
        && tar xf eukaryota_odb10.gz -C busco

# Insecta database
wget -q -O insecta_odb10.gz "https://busco-data.ezlab.org/v4/data/lineages/insecta_odb10.2020-09-10.tar.gz" \
        && tar xf insecta_odb10.gz -C busco

# Hemiptera database
wget -q -O hemiptera_odb10.gz "https://busco-data.ezlab.org/v4/data/lineages/hemiptera_odb10.2020-08-05.tar.gz" \
        && tar xf hemiptera_odb10.gz -C busco

# Most specific Hamiltonela and Buchneria database
wget -q -O enterobacterales_odb10.gz "https://busco-data.ezlab.org/v4/data/lineages/enterobacterales_odb10.2020-03-06.tar.gz" \
        && tar xf enterobacterales_odb10.gz -C busco

# More general Hamiltonela and Buchneria database
wget -q -O gammaproteobacteria_odb10.gz "https://busco-data.ezlab.org/v4/data/lineages/gammaproteobacteria_odb10.2020-03-06.tar.gz" \
        && tar xf gammaproteobacteria_odb10.gz -C busco
````











# Backmap reads to the assembled genomes

Blobtools needs read coverage for each genome assembly, so we align the raw reads from each sample against each metagenome assembly. This runs in less than 10 minutes.

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
REF="assemblies/"$NAME"/"$NAME".p_ctg.fa"

# Map reads to the reference genome
minimap2 -ax map-pb -t "$THREADS" "$REF" "$SAMP" | samtools view -bS -@ $THREADS - | samtools sort -@ $THREADS -o bamsSelfmap/"$NAME".bam

# Index the reference genome (must be .csi format index)
samtools index -c bamsSelfmap/"$NAME".bam
````

# Calculate BUSCO scores for each assembly

Note that the BUSCO software is really inconsistent about relative file paths, hence the inconsistencies in the following script.

````bash
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
REF="/scicore/home/ebertd/dexter0000/aphid/assemblies/"$NAME"/"$NAME".p_ctg.fa"

# Move into output directory
cd busco

# Calculate BUSCO scores
busco -i "$REF" \
-m genome \
-l /scicore/home/ebertd/dexter0000/aphid/blobtools/busco/enterobacterales_odb10 \
-c 8 \
-o $NAME \
-f --offline
````

# Create a minimal blob directory

The minimum requirement to create a new BlobDir is an assembly FASTA file. This command creates a new directory in the location specified by the last argument (in this case “AssemblyName”) that contains a set of files containing values for GC-content (gc.json), length (length.json), number of Ns (ncount.json) and sequence names (identifiers.json) for each sequence in the assembly. A final file (meta.json) contains metadata for the dataset describing the datatypes of the available fields and the ranges of values for each of these fields:

````
blobtools create \
    --fasta /scicore/home/ebertd/dexter0000/aphid/assemblies/S01/S01.p_ctg.fa \
    S01
````



# Add blast hits

````
blastn -db nt \
       -query /scicore/home/ebertd/dexter0000/aphid/assemblies/S01/S01.p_ctg.fa \
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 10 \
       -max_hsps 1 \
       -evalue 1e-25 \
       -num_threads 16 \
       -out S01.ncbi.blastn.out"
````



# Add BUSCO hits

Can add as many BUSCO files as are available

````
blobtools add \
    --busco ASSEMBLY_NAME.busco.nematoda_odb9.full_summary.tsv \
    --busco ASSEMBLY_NAME.busco.metazoa_odb9.full_summary.tsv \
    --busco ASSEMBLY_NAME.busco.eukaryota_odb9.full_summary.tsv \

blobtools add \
    --busco testS01/run_enterobacterales_odb10/full_table.tsv \
    S01
````



# Add the coverage data

````
blobtools add \
    --cov /scicore/home/ebertd/dexter0000/aphid/bamsSelfmap/S01.bam \
    S01
````

# Open the interactive viewer

https://blobtoolkit.genomehubs.org/btk-viewer/viewer-tutorials/hosting-a-local-instance/

````
conda activate btk
blobtools host /path/to/dir/containing/blobdirs

````

````
ssh -L 8000:127.0.0.1:8000 -L 8080:127.0.0.1:8080 dexter0000@login.scicore.unibas.ch
````





# Test data

Initialize the blob

````
blobtools create \
    --fasta testData/assembly.fasta \
    testBlob
````

Add coverage data

````
blobtools add \
    --cov testData/assembly.reads.bam \
    testBlob
````

Add hits data

````
blobtools add \
    --hits testData/blast.out \
    --taxdump taxdump \
    testBlob
````

Open viewer

````
blobtools host .

````

Open local viewer

````
ssh -L 8000:127.0.0.1:8000 -L 8080:127.0.0.1:8080 dexter0000@login.scicore.unibas.ch
````

navigate to page

````
http://localhost:8080
````



````
ssh -L 8001:127.0.0.1:8001 -L 8080:127.0.0.1:8080 dexter0000@login.scicore.unibas.ch
````

````
blobtools view --view snail --host https://blobtoolkit.genomehubs.org mSciVul1_1


````

