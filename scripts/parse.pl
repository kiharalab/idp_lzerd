#!/usr/bin/perl
# Call rosetta fragment picker functions to parse PSI-BLAST binary checkpoint file

use strict;
use warnings;

unless (@ARGV == 3) {
    die <<EO_USAGE;
Convert PSI-BLAST binary chk file to fmt expected by Rosetta Fragment Picker
Usage: script.pl path/to/rosetta/make_fragments.pl myfile.fasta myfile.chk
EO_USAGE
}

my ($make_fragments_path, $fasta_file, $chk_file) = @ARGV;

unless (-s $make_fragments_path) { die("make_fragments.pl not found at $make_fragments_path"); }
unless (-s $fasta_file ) { die("FASTA file $fasta_file not found"); }
unless (-s $chk_file ) { die("PSI-BLAST chk file $chk_file not found"); }

do $make_fragments_path;

my @parts = split(/\./, $chk_file);
pop @parts;
my $file_no_ext = join('.', @parts);
my $outfile = "$file_no_ext.checkpoint";

my $sequence = read_fasta($fasta_file);
my @checkpoint_matrix;
@checkpoint_matrix = &parse_checkpoint_file($chk_file);
@checkpoint_matrix =
  &finish_checkpoint_matrix( $sequence, @checkpoint_matrix );
&write_checkpoint_file($outfile, $sequence,
    @checkpoint_matrix );
