### Sample S01

````
module load SAMtools
module load seqtk
cd /scicore/home/ebertd/dexter0000/aphid/FLYE2/S01

samtools 
mkdir -p curated

#Contig 1
samtools faidx assembly.fasta contig_1:1-538739 > curated/curated_C1.fa

#Contig 4 (APSE)
samtools faidx assembly.fasta contig_4 > curated/curated_C4.fa

#Contig 2
samtools faidx -i assembly.fasta contig_2 > curated/rev_C2.fa
samtools faidx curated/rev_C2.fa contig_2/rc:10701-1400900 > curated/curated_C2.fa

cd curated

cat curated_C1.fa curated_C4.fa curated_C2.fa > chrom.fa

C1: 1-538739
C4: All
C2_RC: 10701-1400900
````

### Sample S04

````
cd /scicore/home/ebertd/dexter0000/aphid/FLYE2/S04

mkdir -p curated

#Contig 2
samtools faidx -i assembly.fasta contig_2 > curated/rev_C2.fa
samtools faidx curated/rev_C2.fa contig_2/rc:1-1918799 > curated/curated_C2.fa

#Contig 3 (APSE)
samtools faidx assembly.fasta contig_3 > curated/curated_C3.fa

#Contig 1
samtools faidx assembly.fasta contig_1:9782-166856 > curated/curated_C1.fa

cd curated

cat curated_C2.fa curated_C3.fa curated_C1.fa > chrom.fa

C2_RC: 1-1918799
C3: All
C1: 9782-166856
````

### S06

````
cd /scicore/home/ebertd/dexter0000/aphid/FLYE2/S04

mkdir -p curated

#Contig 1
samtools faidx -i assembly.fasta contig_1:11838-1406076 > curated/curated_C1.fa

#Contig 6
samtools faidx -i assembly.fasta contig_6 > curated/curated_C6.fa

#Contig 2
samtools faidx assembly.fasta contig_2:10890-549672 > curated/curated_C2.fa

cd curated

cat curated_C1.fa curated_C6.fa curated_C2.fa > chrom.fa


````

### S12

````
cd /scicore/home/ebertd/dexter0000/aphid/FLYE2/S04

mkdir -p curated

#Contig 1
samtools faidx -i assembly.fasta contig_1:8451-1100801 > curated/curated_C1.fa

#Contig 3
samtools faidx -i assembly.fasta contig_3 > curated/curated_C3.fa

#Contig 2
samtools faidx -i assembly.fasta contig_2:1-976513 > curated/curated_C2.fa

cd curated

cat curated_C1.fa curated_C3.fa curated_C2.fa > chrom.fa
````

### S13

````
cd /scicore/home/ebertd/dexter0000/aphid/FLYE2/S04

mkdir -p curated

#Contig 1
samtools faidx -i assembly.fasta contig_1:8545-1927204 > curated/curated_C1.fa

samtools faidx -i assembly.fasta contig_4 > curated/curated_C4.fa

samtools faidx assembly.fasta contig_2:8600-166361 > curated/curated_C2.fa

cd curated

cat curated_C1.fa curated_C4.fa curated_C2.fa > chrom.fa
````

