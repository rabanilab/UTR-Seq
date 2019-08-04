#!/usr/bin/perl

# =============================================================================
# Include
# =============================================================================
use strict;

# =============================================================================
# Main part
# =============================================================================

# reading arguments
if ($ARGV[0] eq "--help") {
  print STDOUT <DATA>;
  exit;
}


my $file_ref;
my $file_name = $ARGV[0];
if (length($file_name) < 1 or $file_name =~ /^-/) {
  $file_ref = \*STDIN;
}
else {
  shift(@ARGV);
  open(FILE, $file_name) or die("Could not open $file_name.\n");
  $file_ref = \*FILE;
}


my $letter_size = 1;
while (scalar(@ARGV) > 0) {
  my $carg = shift(@ARGV);
  if ($carg =~ m/^-k$/g) {
      $letter_size = shift(@ARGV);
  }
}
print STDERR ">> k = $letter_size\n";


# count k-mers
my %kmers_cnt;
my %kmers_seq;
my %kmers_weight;
while (<$file_ref>) {
    chomp $_;
    my ($id, @seq) = split("\t", "$_");

    my %kmers;
    foreach my $s (@seq) {
	my $length = length($s) - $letter_size + 1;
	for (my $i = 0; $i < $length; $i++) {
	    my $kmr = substr($s, $i, $letter_size);
	    $kmers{$kmr}++;
	}
    }

    foreach my $kmr (keys %kmers) {
	$kmers_cnt{$kmr}++;
	$kmers_seq{$kmr}.="$id,";
	$kmers_weight{$kmr}.="$kmers{$kmr},";
    }
}

# print output
foreach my $kmr (keys %kmers_cnt) {
    chop($kmers_weight{$kmr});
    chop($kmers_seq{$kmr});
    print "$kmr\t$kmers_cnt{$kmr}\t$kmers_seq{$kmr}\t$kmers_weight{$kmr}\n";
}


my $T = 4**$letter_size;
my $C = scalar(keys %kmers_cnt);
my $P = int(1000*$C/$T)/10;
print STDERR ">> Expected = $T, observed = $C ($P%)\n";


# =============================================================================
# Subroutines
# =============================================================================


# ------------------------------------------------------------------------
# Help message
# ------------------------------------------------------------------------
__DATA__

extract_kmer_counts.pl <file_name> [options]

Extract all k-mers of a given size in the input sequences.
Input:  <id> <seq>
Output format is a list:
  <kmer> <seq count> [<seq>,<seq>,...] [<weight>,<weight>,...]
  (weight = number of appearances of kmer in seq)

OPTIONS
  -k <num>      K-mer size (default=1)
