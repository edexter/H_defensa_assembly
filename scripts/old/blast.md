Build the BLAST database


````
# Request compute resources
srun --nodes=1 --cpus-per-task=8 --mem=32G --pty bash

# Make sure that a combined fasta does not exist
rm blast/combined_assemblies.fasta

# The genome name has to be added to each contig and then all of the individual fasta concatenated together
for file in assemblies_filtered/*.fa
do
    # Extract the filename without the .fa extension
    name=$(basename "$file" .fa)
    
    # Decompress, modify FASTA headers, and append to combined file
    cat "$file" | awk -v name="$name" '/^>/{print ">" name "_" $0; next} {print}' >> blast/combined_assemblies.fasta
done

# Make a custom database from the concatenated fasta
module load BLAST+
makeblastdb -in blast/combined_assemblies.fasta -dbtype nucl -out blast/hdef_db
````



### Run local BLAST

````
# For APSE
blastn -query blast/APSE.fa -db blast/hdef_db -out blast/APSE_results.txt

grep '^>' blast/APSE_results.txt | sort | uniq

# For PHAT
blastn -query blast/PHD5AT.fa -db blast/hdef_db -out blast/PHD5AT_results.txt

grep '^>' blast/PHD5AT_results.txt | sort | uniq
````



````
srun --nodes=1 --cpus-per-task=8 --mem=32G --pty bash
module load BLAST+
cd aphid

blastn -query blast/APSE.fa -db blast/hdef_db -out blast/APSE_results.txt  -outfmt "6 qseqid staxids bitscore std"

blastn -query blast/PHD5AT.fa -db blast/hdef_db -out blast/PHD5AT_results.txt  -outfmt "6 qseqid staxids bitscore std"

blastn -query blast/M147_plasmid.fasta -db blast/hdef_db -out blast/M147_results.txt  -outfmt "6 qseqid staxids bitscore std"

blastn -query blast/plasmid_P4M47.fa -db blast/hdef_db -out blast/P4M47_results.txt  -outfmt "6 qseqid staxids bitscore std"
````

