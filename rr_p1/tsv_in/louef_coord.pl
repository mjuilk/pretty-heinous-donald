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
open(my $in,  "<", "gnomad.v4.1.constraint_metrics.tsv") or die "Cannot open input file: $!";
open(my $out, ">", "gnomad.v4.1.constraint_metrics_with_coordinates.tsv") or die "Cannot open output file: $!";

# Print the header line from the input file to the output file
my $header = <$in>;
print $out $header;

# Process each line in the input file
while (my $line = <$in>) {
    chomp $line;
    my @fields = split(/\t/, $line);
    my $gene_id = $fields[0];  # Assuming the gene ID is in the first column

    # Fetch the gene object from Ensembl
    my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Human', 'Core', 'Gene');
    my $gene = $gene_adaptor->fetch_by_stable_id($gene_id);

    if ($gene) {
        # Append the start and end positions to the line
        push @fields, $gene->start, $gene->end;
        print $out join("\t", @fields) . "\n";
    } else {
        warn "Gene $gene_id not found in Ensembl\n";
    }
}

close($in);
close($out);

