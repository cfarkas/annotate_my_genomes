#!/usr/bin/perl

use warnings;
use strict;
use LWP::UserAgent;
use HTTP::Date;

# Check that a taxonomy identifier was passed as a command line argument
if (!$ARGV[0]) {
  die "Error: No taxonomy identifier specified.\nUsage: perl download_proteome_uniprot.pl <taxonomy_id>\n";
}

# Taxonomy identifier of top node for query, e.g. 2 for Bacteria, 2157 for Archaea, etc.
# (see https://www.uniprot.org/taxonomy)
my $top_node = $ARGV[0];

# Create a user agent for making HTTP requests
my $agent = LWP::UserAgent->new;

# Get a list of all reference proteomes of organisms below the given taxonomy node.
my $query_list = "https://rest.uniprot.org/proteomes/stream?query=reference:true+taxonomy_id:$top_node&format=list";

my $response_list = $agent->get($query_list);

# Check for HTTP errors
if (!$response_list->is_success) {
  die 'Failed to get proteome list: ' . $response_list->status_line .
    ' for ' . $response_list->request->uri . "\n";
}

# For each proteome, mirror its set of UniProt entries in compressed FASTA format.
for my $proteome (split(/\n/, $response_list->content)) {
  my $file = $proteome . '.fasta.gz';
  my $query_proteome = "https://rest.uniprot.org/uniprotkb/stream?query=proteome:$proteome&format=fasta&compressed=true";
  my $response_proteome = $agent->mirror($query_proteome, $file);

  # Check for HTTP errors
  if ($response_proteome->is_success) {
    my $release = $response_proteome->header('x-uniprot-release');
    my $date = $response_proteome->header('x-uniprot-release-date');
    print "File $file: downloaded entries of UniProt release $release ($date)\n";
  }
  elsif ($response_proteome->code == HTTP::Status::RC_NOT_MODIFIED) {
    print "File $file: up-to-date\n";
  }
  else {
    die 'Failed to download proteome: ' . $response_proteome->status_line .
      ' for ' . $response_proteome->request->uri . "\n";
  }
}
