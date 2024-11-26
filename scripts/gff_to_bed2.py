import argparse
import os

parser = argparse.ArgumentParser(description='Convert gff file to bed5 file')
parser.add_argument('gff_file', type=str, help='input gff file')
parser.add_argument('bed_file', type=str, help='output bed5 file')
args = parser.parse_args()

if not os.path.isfile(args.gff_file):
    raise ValueError(f'Invalid file path, file not found: {args.gff_file}')

# file format should be gff or gff3
if not args.gff_file.endswith('.gff') and not args.gff_file.endswith('.gff3'):
    raise ValueError(f'Invalid file format, should be gff or gff3: {args.gff_file}')


def has_gene_annotation(gff_file):
    with open(gff_file, 'r') as f:
        for line in f:
            # Skip comments and empty lines
            if line.startswith('#') or line.strip() == '':
                continue
            fields = line.strip().split('\t')
            # Check if feature is gene
            if fields[2] == 'gene':
                return True

    # No gene features found
    return False


def parse_gff_gene(gff_file, bed_file):
    with open(gff_file, 'r') as f_in, open(bed_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[2].lower() != 'gene':
                continue
            chr_name = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8].split(';')
            gene_id = None
            for attr in attributes:
                if attr.startswith('ID=gene:'):
                    gene_id = attr[8:]
                    break
                elif attr.startswith('ID=gene-'):
                    gene_id = attr[8:]
                    break
                elif attr.startswith('ID='):
                    gene_id = attr[3:]
                    break
                else:
                    continue

            if gene_id is None:
                print(f'No ID found: {line}')
                continue
            f_out.write(f'{chr_name}\t{gene_id}\t{start}\t{end}\t{strand}\n')


def parse_gff_mrna(gff_file, bed_file):
    with open(gff_file, 'r') as f_in, open(bed_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[2].lower() != 'mrna':
                continue
            chr_name = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8].split(';')
            gene_id = None
            for attr in attributes:
                if attr.startswith('ID=gene:'):
                    gene_id = attr[8:]
                    break
                elif attr.startswith('ID=gene-'):
                    gene_id = attr[8:]
                    break
                elif attr.startswith('ID='):
                    gene_id = attr[3:]
                    break
                else:
                    continue

            if gene_id is None:
                print(f'No ID found: {line}')
                continue
            f_out.write(f'{chr_name}\t{gene_id}\t{start}\t{end}\t{strand}\n')
            
            # Remove the "_gene" suffix if it exists
            if gene_id.endswith('_gene'):
                gene_id = gene_id[:-5]


# main function
gff_file = args.gff_file
bed_file = args.bed_file
if has_gene_annotation(gff_file):
    parse_gff_gene(gff_file, bed_file)
else:
    parse_gff_mrna(gff_file, bed_file)


# Usage method：
# ```bash
# python gff3_to_bed.py input.gff3 output.bed
# ```
# Note: 'input.gff3' is the input gff3 file and 'output.bed' is the output bed5 file.
