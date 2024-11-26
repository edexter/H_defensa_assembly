```
srun --nodes=1 --cpus-per-task=8 --mem=32G --pty bash

module load circlator

cd circlator
circlator all /scicore/home/ebertd/dexter0000/aphid/FLYE2/S01/assembly.fasta /scicore/home/ebertd/dexter0000/aphid/reads_filtered/S01.fq S01
```

