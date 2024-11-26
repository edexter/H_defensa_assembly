# Set up Amazon EC2 instance

First, set up the EC2 instance and download the access key file into local home directory. Note that all of these commands are for "Amazon 2023 Linux", and are quite different from the syntax that one finds on the internet for the older "Amazon 2 Linux".

````
# Set permission level for key
chmod 400 pgap-key.pem

# Login with keypair
ssh -i pgap-key.pem ec2-user@3.80.203.237

# Start with the usual update 
sudo yum update -y

# Install docker
sudo yum install -y docker

# Start docker
sudo service docker start

# Check docker version
docker --version

# Make sure docker always starts when you start this instance
sudo systemctl enable docker

# Give the ec2-user permission to use docker (requires log-out to take effect)
sudo usermod -a -G docker ec2-user
exit

# Log back in
ssh -i pgap-key.pem ec2-user@3.80.203.237

# Check docker version
docker --version

# Check docker access
docker ps

# Download the docker image
docker pull ncbi/pgap-utils:2024-07-18.build7555

# Confirm image is available
docker images

# Install curl (Amazon 2023 OS specific)
sudo dnf swap curl-minimal curl-full

# Get the rest of the PGAP install that is not included with the docker image
curl -OL https://github.com/ncbi/pgap/raw/prod/scripts/pgap.py

# Make the script executable
chmod +x pgap.py

# Run the update script to get the rest of the files (takes ~5 minutes)
./pgap.py --update

# Run the test data
./pgap.py -r -o mg37_results -g $HOME/.pgap/test_genomes/MG37/ASM2732v1.annotation.nucleotide.1.fasta -s 'Mycoplasmoides genitalium'
````

### Transfer files to server

````bash
scp -i ~/pgap-key.pem \
-r /mnt/c/Users/ericd/Downloads/genomes \
ec2-user@3.80.203.237:/home/ec2-user/
````

# Transfer files from server

````
rsync -av -e "ssh -i ~/pgap-key.pem" ec2-user@3.80.203.237:/home/ec2-user/S* /mnt/c/Users/ericd/Downloads/results/

rsync -av -e "ssh -i ~/pgap-key.pem" ec2-user@3.80.203.237:/home/ec2-user/buchnera* /mnt/c/Users/ericd/Downloads/results/
````

# Run using genomes

````bash
# Buchnera
./pgap.py -r \
-o buchnera_results \
-g $HOME/genomes/buchnera.curated.fa \
-s 'Buchnera aphidicola'

#S01
./pgap.py -r \
-o S01_results \
-g $HOME/genomes/S01.curated.fasta \
-s 'Candidatus Hamiltonella defensa'

#S05
./pgap.py -r \
-o S05_results \
-g $HOME/genomes/S05.curated.fasta \
-s 'Candidatus Hamiltonella defensa'

#S06
./pgap.py -r \
-o S06_results \
-g $HOME/genomes/S06.curated.fasta \
-s 'Candidatus Hamiltonella defensa'

#S06
./pgap.py -r \
-o S06_results \
-g $HOME/genomes/S06.curated.fasta \
-s 'Candidatus Hamiltonella defensa'

#S08
./pgap.py -r \
-o S08_results \
-g $HOME/genomes/S08.curated.fasta \
-s 'Candidatus Hamiltonella defensa'

# Run across all samples with nohup
nohup bash -c '
while IFS= read -r genome_file; do
    # Extract the base name without the extension for the output directory name
    base_name=$(basename "$genome_file" .curated.fasta)
    
    # Set the output directory
    output_dir="${base_name}_results"
    
    # Run PGAP
    ./pgap.py -r \
        -o "$output_dir" \
        -g "$HOME/genomes/$genome_file" \
        -s "Candidatus Hamiltonella defensa"
    
    echo "Completed processing for $genome_file. Results saved in $output_dir."
done < index.txt' > pgap_batch.log 2>&1 &

#
nohup bash -c '
while IFS= read -r genome_file; do
    # Extract the base name without the extension for the output directory name
    base_name=$(basename "$genome_file" .curated.fasta)
    
    # Set the output directory
    output_dir="${base_name}_results"
    
    # Run PGAP and log output for each sample separately
    ./pgap.py -r \
        -o "$output_dir" \
        -g "$HOME/genomes/$genome_file" \
        -s "Candidatus Hamiltonella defensa" > "${output_dir}.log" 2>&1

    echo "Completed processing for $genome_file. Results saved in $output_dir."
done < index.txt' > pgap_batch.log 2>&1 &

````

# Job scipt

````
#!/bin/bash

#REDO
./pgap.py -r \
-o S08_results \
-g $HOME/genomes/S08.curated.fasta \
-s 'Candidatus Hamiltonella defensa'

./pgap.py -r \
-o S09_results \
-g $HOME/genomes/S09.curated.fasta \
-s 'Candidatus Hamiltonella defensa'

./pgap.py -r \
-o S10_results \
-g $HOME/genomes/S10.curated.fasta \
-s 'Candidatus Hamiltonella defensa'

./pgap.py -r \
-o S11_results \
-g $HOME/genomes/S11.curated.fasta \
-s 'Candidatus Hamiltonella defensa'

./pgap.py -r \
-o S12_results \
-g $HOME/genomes/S12.curated.fasta \
-s 'Candidatus Hamiltonella defensa'

./pgap.py -r \
-o S13_results \
-g $HOME/genomes/S13.curated.fasta \
-s 'Candidatus Hamiltonella defensa'

./pgap.py -r \
-o S14_results \
-g $HOME/genomes/S14.curated.fasta \
-s 'Candidatus Hamiltonella defensa'
````

````
nohup bash batch.sh > pgap_batch_main.log 2>&1 &
````

