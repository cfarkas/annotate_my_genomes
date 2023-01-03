#!/usr/bin/perl

use strict;
use warnings;
use LWP::UserAgent;
use HTTP::Date;
use HTTP::Status;  # Add "use" statement for HTTP::Status module

# Print usage statement
if (@ARGV != 1) {
  print STDERR "Usage: $0 TAXON_ID\n";
  print STDERR "Downloads FASTA records for the specified taxonomy ID from UniProt.\n";
  exit 1;  # Exit with non-zero status code to indicate error
}

my $taxon = $ARGV[0]; # Taxonomy identifier of organism.

my $query = "https://rest.uniprot.org/uniprotkb/search?query=organism_id:$taxon&format=fasta";

my $file = $taxon . '.fasta';

my $contact = ''; # Please set a contact email address here to help us debug in case of problems (see https://www.uniprot.org/help/privacy).
my $agent = LWP::UserAgent->new(agent => "libwww-perl $contact");
my $response = $agent->mirror($query, $file);

if ($response->is_success) {
  # Print number of results and release info
  my $results = $response->header('X-Total-Results');
  my $release = $response->header('X-UniProt-Release');
  my $date = sprintf("%4d-%02d-%02d", HTTP::Date::parse_date($response->header('X-UniProt-Release-Date')));
  print "Downloaded $results entries of UniProt release $release ($date) to file $file\n";
}
elsif ($response->code == HTTP::Status::RC_NOT_MODIFIED) {
  # Print message if data is already up-to-date
  print "Data for taxon $taxon is up-to-date.\n";
}
else {
  # Print error message and exit with non-zero status code
  warn 'Failed, got ' . $response->status_line .
    ' for ' . $response->request->uri . "\n";
  exit 1;
}
