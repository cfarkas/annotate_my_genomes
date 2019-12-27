#!/usr/bin/perl -w

use strict;
use Getopt::Std;
use lib '.';
use GTF;

use vars qw($opt_c $opt_e $opt_f $opt_p $opt_s $opt_t $opt_l $opt_m $opt_o);
getopts('ce:fpst:mlo');
my $usage = "usage: $0 [-fsmbB] [-t tx output filename] <gtf file> [sequence file]
Options:
 -t <file>: output transcript file
 -f: create a fixed gtf file (This may not be possible.  
     Always check the \"fixed\" file) 
 -o: create one gene with multiple transcript from overlapping genes. 
     Useful for cleaning UCSC RefSeq or MGC downloads
     only works with -f flag!
 -e <count>: sets the maximum number of detailed error messages to return per 
     error to <count> (default is 5).
 -s: output list of inframe stop genes.
 -c: suppress warnings about missing start/stop
 -p: suppress warnings about bad splice site sequence
 -l: output a list of bad genes for training applications
 -m: output a list of bad genes for evaluation purposes
";

die $usage unless ((@ARGV == 1)||(@ARGV == 2));
my ($filename,$seqname) = @ARGV;
if($opt_t){
    unless(defined($seqname)){
	die "Must give a sequence file when you give the -t option.\n";
    }
    open(TX,">$opt_t");
}
else{
    open(TX,">/dev/null");
}
my $fix_file = $opt_f;
my $fix_filename;
my $suppress = [];
if($opt_c){
    $$suppress[17] = 1;
    $$suppress[19] = 1;
}
if($opt_p){
    $$suppress[26] = 1;
    $$suppress[27] = 1;
}
my $bad_genes = [];
if($opt_m){
	$$bad_genes[24] = 1;
	$$bad_genes[25] = 1;
	$$bad_genes[26] = 1;
	$$bad_genes[27] = 1;
	$$bad_genes[28] = 1;
	$$bad_genes[29] = 1;
	$$bad_genes[37] = 1;
	$$bad_genes[38] = 1;
	$$bad_genes[41] = 1;
	$$bad_genes[48] = 1;
}
if($opt_l){
	$$bad_genes[16] = 1;
	$$bad_genes[18] = 1;
	$$bad_genes[21] = 1;
	$$bad_genes[22] = 1;
	$$bad_genes[24] = 1;
	$$bad_genes[25] = 1;
	$$bad_genes[26] = 1;
	$$bad_genes[27] = 1;
	$$bad_genes[28] = 1;
	$$bad_genes[29] = 1;
	$$bad_genes[33] = 1;
	$$bad_genes[34] = 1;
	$$bad_genes[35] = 1;
	$$bad_genes[36] = 1;
	$$bad_genes[37] = 1;
	$$bad_genes[38] = 1;
	$$bad_genes[41] = 1;
	$$bad_genes[48] = 1;
}
my $info = {gtf_filename  => $filename,
	    warning_fh    => \*STDOUT,
	    fix_gtf       => $fix_file,
	    seq_filename  => $seqname,
	    inframe_stops => $opt_s,
	    tx_out_fh     => \*TX,
	    warning_skips => $suppress,
	    bad_list      => $bad_genes};
if($opt_e){
    if($opt_e =~ /^(\d+)$/){
	$info->{detailed_error_count} = $opt_e;
    }
    else{
	print STDERR "Bad value, $opt_e, for -e flag.  Should be an integer.\n";
    }
}

my $gtf = GTF::new($info);

close(TX);
if($fix_file){
    if($filename =~ /(\S*).g[t,f]f/){
	$fix_filename = "$1.fixed.gtf";
    }
    else{
	$fix_filename = "$filename.fixed.gtf";
    }
    if(open(FIX, ">$fix_filename")){
	if($opt_o){
		$gtf = cleanup_overlaps($gtf);
		$gtf->output_gtf_file(\*FIX);
          
        }else{
		$gtf->output_gtf_file(\*FIX);
	}
    }
    else{
	$fix_file = 0;
	print STDERR "Could not open $fix_filename for output.\n";
	print STDERR "Will not create fixed gtf file.\n";
    }
}

exit;

  
sub cleanup_overlaps{
  my($gtf) = (@_);
  
  my%hs = ();
  my%nonexon_overlap = ();
  my%allgenes = ();
  my$previous_end = "0";
  my@removals = ();

  my $genes = $gtf->genes;
  foreach my$gene(@$genes){
    my$gene_id = $gene->id;
    my$gene_start = $gene->start;
    my$gene_stop = $gene->stop;
    my$seqname = $gene->seqname;
#    if($seqname =~ /random/){
#      print "$gene_id is located on random chromosome\n";
#      $gtf -> remove_gene($gene_id);
#      push(@removals, $gene_id);
#    }
  
    if(exists($allgenes{$gene_id})){
      die "$gene_id already exists in file";
    }
    my $transcripts = $gene->transcripts;
    foreach my$transcript(@$transcripts){
      my$nr_of_exons = "0";
      my$transcript_id = $transcript->id;
      my$orientation = $transcript->strand;
  
# we're only interested in the CDS of the RefSeqs
      my$exons = $transcript->cds;
      my$exonnr = scalar @$exons;
                                                                                                                               
# GTF knows to take the start of the first or the end of the last
# exon if there is no start or stop codon.
      my$start = $transcript->start;
      my$stop = $transcript->stop;
      foreach my$exon(@$exons){
        my$return = "";
        $nr_of_exons++;
        my$exonstart = $exon->start;
        my$exonend = $exon->stop;
# check if exonstart already exists. If so, change transcriptname and exit loop
# unless of course the transcript is already part of the correct gene
        if(exists ($hs{$exonstart.$exonend})){
          next if ($hs{$exonstart.$exonend} -> gene_id eq $gene_id);
# see if gene already has a transcript that is identical to the current transcript
          my$parent_transcripts = $hs{$exonstart.$exonend}-> transcripts;
          my$sameflag = "0";
          foreach my$tx(@$parent_transcripts){
            if($transcript ->equals($tx)){
              $sameflag = "1";
              print "$gene_id is identical to ".$tx->gene_id." of gene ".$hs{$exonstart.$exonend}->gene_id."\n";
              $gtf -> remove_gene($gene_id);
              push(@removals, $gene_id);
              last;
            }
          }
          if ($sameflag){
            if($gene_stop > $previous_end){
              $previous_end = $gene_stop;
            }
            last;
          }
    
# if the transcript is not identical to an existing one,
# add transcript to parent gene
          my$parent_gene = $hs{$exonstart.$exonend} -> gene_id;
          $hs{$exonstart.$exonend} -> add_transcript($transcript);
          print "added $gene_id as transcript of $parent_gene\n";
# and remove current gene
          $gtf -> remove_gene($gene_id);
          push(@removals, $gene_id);
          last;
# Exon not found in any other transcript, create entry for exonstart
        }else{
          $hs{$exonstart.$exonend} = $gene;
          if($gene_start <= $previous_end){
            $nonexon_overlap{$gene_id} = $exonnr;
          }
        }
      }
    }
    if($gene_stop > $previous_end){
      $previous_end = $gene_stop;
    }
  }
  
# of all genes that have exons in introns of other genes, remove the ones that are single exon genes
# so multiexon genes that overlap other genes are kept in the set! THIS DOESN'T WORK PROPERLY
#  my@nonexon_overlap = sort keys %nonexon_overlap; 
#  my@overlap_nonexon = subtract_array(\@nonexon_overlap,\@removals); 
#  foreach my$id(@overlap_nonexon){
#    my$exonnr = $nonexon_overlap{$id};
#    if($exonnr == 1){
#      print "$id is single exon gene in intron of multiexon gene:removed\n";
#      $gtf -> remove_gene($id);
#      push(@removals, $id);
#    }else{
#      print "$id multiexon gene in intron(s) of other genes: kept in set\n";
#    }
#  }
  
#  $gtf->output_gtf_file;
  
  foreach(@removals){
    print "removed\t$_\n";
  }
  return $gtf;
}  
  
  
# subtract_array takes in 2 arrays and removes all entries from the first array
# that appear in the second array
  
sub subtract_array{
  my($total, $subset) = (@_);
  my%hs = ();
  my@rest = ();
   
  foreach my$entry(@$subset){
    $hs{$entry}="1";
  }
  foreach my$line(@$total){
    unless(exists$hs{$line}){
      push(@rest, $line);
    }
  }
  return @rest;
}

