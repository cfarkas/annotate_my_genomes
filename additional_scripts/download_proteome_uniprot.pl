use strict;
use warnings;
use LWP::UserAgent;
use HTTP::Date;

my $taxon = $ARGV[0]; # Taxonomy identifier of organism.

my $query = "https://rest.uniprot.org/uniprotkb/search?query=organism_id:$taxon&format=fasta";
    
my $file = $taxon . '.fasta';

my $contact = ''; # Please set a contact email address here to help us debug in case of problems (see https://www.uniprot.org/help/privacy).
my $agent = LWP::UserAgent->new(agent => "libwww-perl $contact");
my $response = $agent->mirror($query, $file);

if ($response->is_success) {
  my $results = $response->header('X-Total-Results');
  my $release = $response->header('X-UniProt-Release');
  my $date = sprintf("%4d-%02d-%02d", HTTP::Date::parse_date($response->header('X-UniProt-Release-Date')));
  print "Downloaded $results entries of UniProt release $release ($date) to file $file\n";
  }
  elsif ($response->code == HTTP::Status::RC_NOT_MODIFIED) {
    print "Data for taxon $taxon is up-to-date.\n";
  }
  else {
    die 'Failed, got ' . $response->status_line .
      ' for ' . $response->request->uri . "\n";
  }
