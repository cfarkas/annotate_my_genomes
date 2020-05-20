#!/usr/bin/perl
use warnings;
use strict;

my $filename = 'perl.script';

open(FH, '<', $filename) or die $!;

while(<FH>){
   do perl -pe perl.script replace.gtf > merged.gtf;
}

close(FH);