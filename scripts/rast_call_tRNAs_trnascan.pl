#!/usr/bin/env perl
# This is a SAS Component
########################################################################
# Copyright (c) 2003-2013 University of Chicago and Fellowship
# for Interpretations of Genomes. All Rights Reserved.
#
# This file is part of the SEED Toolkit.
# 
# The SEED Toolkit is free software. You can redistribute
# it and/or modify it under the terms of the SEED Toolkit
# Public License. 
#
# You should have received a copy of the SEED Toolkit Public License
# along with this program; if not write to the University of Chicago
# at info@ci.uchicago.edu or the Fellowship for Interpretation of
# Genomes at veronika@thefig.info or download a copy from
# http://www.theseed.org/LICENSE.TXT.
########################################################################

=head1 NAME

rast_call_RNAs

=head1 SYNOPSIS

rast_call_RNAs [--input genome_file] [--output genome_file] [--url service-url] [< genome-file] [> genome-file]   

=head1 DESCRIPTION

Find instances of tRNAs in a genome-type object. Uses tRNAscan-SE directly

Example:

    rast_call_tRNAs < input_genome > output_genome_with_tRNAs_called

=head1 COMMAND-LINE OPTIONS

Usage: rast_call_tRNAs  < input_genome_object  > output_genome_object
Usage: rast_call_tRNAs  --input input_genome_object --output output_genome_object

    --input      --- Read input genome-typed object from file instead of STDIN

    --output     --- Read output genome-typed object from file outstead of STDOUT

    --tmpdir     --- Use named temporary-file directory instead of the default temporary directory

    --id_prefix  --- Use the specified feature prefix instead of the default of 'rast|0'

=head1 AUTHORS

L<The SEED Project|http://www.theseed.org>

=cut


use strict;
use warnings;
use Data::Dumper;
use File::Copy;
use Time::HiRes 'gettimeofday';

use gjoseqlib;
use Bio::KBase::GenomeAnnotation::Client;
eval {
    require  Bio::KBase::IDServer::Client;
};
use JSON::XS;

use IDclient;
use Find_RNAs;
use GenomeTypeObject;

my $help;
my $input_file;
my $output_file;
my $temp_dir;
my $id_prefix = 'rast|0';
my $id_server;

use Getopt::Long;
my $rc = GetOptions('help'        => \$help,
		    'input=s'     => \$input_file,
		    'output=s'    => \$output_file,
		    'tmpdir=s'    => \$temp_dir,
		    'id-prefix=s' => \$id_prefix,
		    'id-server=s' => \$id_server,
		    );

if (!$rc || $help || @ARGV != 0) {
    seek(DATA, 0, 0);
    while (<DATA>) {
	last if /^=head1 COMMAND-LINE /;
    }
    while (<DATA>) {
	last if (/^=/);
	print $_;
    }
    exit($help ? 0 : 1);
}

my $in_fh;
if ($input_file) {
    open($in_fh, "<", $input_file) or die "Cannot open $input_file: $!";
} else { $in_fh = \*STDIN; }

my $out_fh;
if ($output_file) {
    open($out_fh, ">", $output_file) or die "Cannot open $output_file: $!";
} else { $out_fh = \*STDOUT; }

my $json = JSON::XS->new;

my $genomeTO = GenomeTypeObject->create_from_file($in_fh);

my $id_client;
if ($id_server)
{
    $id_client = Bio::KBase::IDServer::Client->new($id_server);
}
else
{	
    $id_client = IDclient->new($genomeTO);
}

my($domain) = ($genomeTO->{domain} =~ m/^([AB])/io);
my $dflag = $domain ? "-$domain" : "-G";

my $contigs_tmp = $genomeTO->extract_contig_sequences_to_temp_file();
close($contigs_tmp);

my $scan_out = File::Temp->new();
close($scan_out);
my $scan_log = File::Temp->new();
close($scan_log);
my $scan_stats = File::Temp->new();
close($scan_stats);

my @cmd = ("tRNAscan-SE", $dflag,
	   "-o", "$scan_out",
	   "-m", "$scan_stats",
	   "-l", "$scan_log",
	   "-Q",
	   "-b",
	   "$contigs_tmp");

my $hostname = `hostname`;
chomp $hostname;

my $event = {
    tool_name => 'tRNAscan-SE',
    execute_time => scalar gettimeofday,
    parameters => \@cmd,
    hostname => $hostname,
};
my $event_id = $genomeTO->add_analysis_event($event);

print STDERR "@cmd\n";
$rc = system(@cmd);
print STDERR "==============================\ntRNAscan-SE stats:\n==============================\n";
copy("$scan_stats", \*STDERR);
print STDERR "==============================\ntRNAscan-SE log:\n==============================\n";
copy("$scan_log", \*STDERR);

if ($rc != 0)
{
    die "tRNAscan-SE failed with $rc\n";
}

open(RES, "<", "$scan_out") or die "Cannot open output file $scan_out: $!";

while (<RES>)
{
    chomp;
    my($contig, $idx, $beg, $end, $type, $anti, $intro_beg, $intron_end, $cove_score) = split(/\t/);
    $contig =~ s/^\s*//;
    $contig =~ s/\s*$//;

    $beg += 0;
    $end += 0;

    my $func = join("-", "tRNA", $type, $anti);

    my $length = 1 + abs($end - $beg);
    my $strand = ($beg < $end) ? q(+) : q(-);
    
    $genomeTO->add_feature({ -id_client  => $id_client,
						-id_prefix  => $id_prefix,
						-type       => 'rna',
						-location   => [[ $contig, $beg, $strand, $length ]],
						-function   => $func,
						-annotator  => 'search_for_rnas',
						-annotation => "Add feature called by tRNAscan-SE with score $cove_score",
						-analysis_event_id => $event_id,
					    }
				   );
}

$genomeTO->destroy_to_file($out_fh);
close($out_fh);

__DATA__
