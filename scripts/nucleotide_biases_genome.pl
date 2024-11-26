#!/usr/bin/perl
## Pombert Lab, 2022
my $name = 'nucleotide_biases.pl';
my $version = '0.4';
my $updated = '2023-09-04';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use File::Path qw(make_path);

my $usage = <<"OPTIONS";
NAME        $name
VERSION     $version
UPDATED     $updated
SYNOPSIS    Generates tab-delimited sliding windows of GC, AT, purine, and pyrimidine
             distributions for easy plotting with MS Excel or other tools.

COMMAND     $name \\
          -f *.fasta \\
          -o output_file.tsv \\
          -w 1000 \\
          -s 500

-f (--fasta)    Fasta file(s) to process
-o (--outfile)  Output file name [Default: combined_gc_content.tsv]
-w (--winsize)  Sliding window size [Default: 1000]
-s (--step)     Sliding window step [Default: 500]
OPTIONS
die "\n$usage\n" unless @ARGV;

my @fasta;
my $outfile = 'combined_gc_content.tsv'; # Default output file name
my $winsize = 1000;
my $step = 500;
GetOptions(
    'f|fasta=s@{1,}' => \@fasta,
    'o|outfile=s' => \$outfile,
    'w|winsize=i' => \$winsize,
    's|step=i' => \$step
);

### Open the specified output file
open my $combined_fh, '>', $outfile or die "Can't create $outfile: $!\n";

### Print header to the combined output file
print $combined_fh "# Contig\tLocation\t% GC\t% AT\t% Purines\t% Pyrimidines\t% GT\t% AC\n";

### Iterating through FASTA file(s)
while (my $fasta = shift @fasta){

    my ($basename) = fileparse($fasta);
    my ($fileprefix) = $basename =~ /(\S+)\.\w+$/;

    open FASTA, "<", $fasta or die "Can't open $fasta: $!\n";

    ### Creating database of sequences (could be multifasta)
    my %sequences;
    my $seqname;

    while (my $line = <FASTA>){
        chomp $line;
        if ($line =~ /^>(\S+)/){
            $seqname = $1;
        }
        else {
            $sequences{$seqname} .= $line;
        }
    }

    ### Iterating through each sequence in the FASTA file
    foreach my $sequence (sort (keys %sequences)){

        ### Sliding windows
        my $seq = $sequences{$sequence};
        my $csize = length $seq;
        my $x;
        for ($x = 0; $x <= ($csize - $step); $x += $step){
            my $subseq = substr($seq, $x, $winsize);
            my $gc = $subseq =~ tr/GgCc//;
            my $at = $subseq =~ tr/AaTt//;
            my $pur = $subseq =~ tr/GgAa//;
            my $pyr = $subseq =~ tr/CcTt//;
            my $gt = $subseq =~ tr/GgTt//;
            my $ac = $subseq =~ tr/AaCc//;
            $gc = ($gc)/$winsize * 100;
            $at = ($at)/$winsize * 100;
            $pur = ($pur)/$winsize * 100;
            $pyr = ($pyr)/$winsize * 100;
            $gt = ($gt)/$winsize * 100;
            $ac = ($ac)/$winsize * 100;
            print $combined_fh "$sequence\t$x\t$gc\t$at\t$pur\t$pyr\t$gt\t$ac\n";
        }

        ### Working on leftover string < $winsize
        my $modulo = $csize % $winsize;
        my $subseqleft = substr($seq, -$modulo, $modulo);
        my $leftover_size = length $subseqleft;
        my $gc = $subseqleft =~ tr/GgCc//;
        my $at = $subseqleft =~ tr/AaTt//;
        my $pur = $subseqleft =~ tr/GgAa//;
        my $pyr = $subseqleft =~ tr/CcTt//;
        my $gt = $subseqleft =~ tr/GgTt//;
        my $ac = $subseqleft =~ tr/AaCc//;
        $gc = ($gc)/$leftover_size * 100;
        $at = ($at)/$leftover_size * 100;  # Corrected from $winsize
        $pur = ($pur)/$leftover_size * 100;  # Corrected from $winsize
        $pyr = ($pyr)/$leftover_size * 100;  # Corrected from $winsize
        $gt = ($gt)/$leftover_size * 100;  # Corrected from $winsize
        $ac = ($ac)/$leftover_size * 100;  # Corrected from $winsize
        print $combined_fh "$sequence\t$x\t$gc\t$at\t$pur\t$pyr\t$gt\t$ac\n";

    }

    close FASTA;

}

### Close the concatenated output file
close $combined_fh;
