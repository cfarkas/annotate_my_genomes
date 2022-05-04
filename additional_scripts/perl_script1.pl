#!/usr/bin/env perl
perl -lne 'print "@m" if @m=(/((?:transcript_id|gene_id)\s+\S+)/g);' coding_transcripts.gtf > final_annotated.tab

