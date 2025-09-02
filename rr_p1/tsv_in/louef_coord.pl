#!/usr/bin/perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;

# Load the Ensembl registry from the database
Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'useastdb.ensembl.org',
    -user => 'anonymous'
);

# Open the input and output files
open(my $in,  "<", "/Users/jon/Documents/mats/rsrcs/gnomad.v4.1.constraint_metrics.tsv") or die "Cannot open input file: $!";
open(my $out, ">", "/Users/jon/Documents/mats/rsrcs/gnomad.v4.1.constraint_metrics_with_coordinates.tsv") or die "Cannot open output file: $!";

# Print the header line from the input file to the output file
my $header = <$in>;
chomp $header;
print $out $header;

# Get adaptors
my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Human', 'Core', 'Gene');
my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Human', 'Core', 'Slice');

# Process each line in the input file
while (my $line = <$in>) {
    chomp $line;
    my @fields = split(/\t/, $line);
    my $gene_symbol = $fields[0];  # Gene symbol is in the first column

    # Fetch all genes with this symbol
    my $genes = $gene_adaptor->fetch_all_by_external_name($gene_symbol);

    if (@$genes) {
        # Use the first gene if multiple exist (or handle as needed)
        my $gene = $genes->[0];
        push @fields, $gene->start, $gene->end;
        print $out join("\t", @fields) . "\n";
    } else {
        # Log missing genes
        warn "Gene symbol $gene_symbol not found in Ensembl\n";
    }
}

close($in);
close($out);

