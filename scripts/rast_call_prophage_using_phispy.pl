#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use SeedUtils;

use File::Slurp;
use Time::HiRes 'gettimeofday';
use File::Temp;
use File::Basename;
use Bio::KBase::GenomeAnnotation::Client;
use JSON::XS;

use IDclient;
use Prodigal;
use GenomeTypeObject;
eval {
    require  Bio::KBase::IDServer::Client;
};

my $help;
my $input_file;
my $output_file;
my $temp_dir;
my $id_prefix;
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
    die "Bad ARGV";
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

if ($genomeTO->{domain} !~ m/^([ABV])/o) {
    die "Invalid domain: \"$genomeTO->{domain}\"";
}

my $gto_file = File::Temp->new(UNLINK => 1);
close($gto_file);
$genomeTO->destroy_to_file("$gto_file");

$genomeTO = GenomeTypeObject->create_from_file("$gto_file");


my $id_client;
if ($id_server)
{
    $id_client = Bio::KBase::IDServer::Client->new($id_server);
}
else
{	
    $id_client = IDclient->new($genomeTO);
}

my $gb_file = File::Temp->new(UNLINK => 1);
close($gb_file);

$rc = system("rast_export_genome", "-i", "$gto_file", "-o", "$gb_file", "genbank");
die "Failed exporting genome: rc=$rc\n" unless $rc == 0;

my $phispy = "PhiSpy.py";

my $sci = $genomeTO->{scientific_name} // "";

my $training_set;
my $default_training_set;

# PhiSpy.py -l long | cat
# Training Set	# of genomes used
# data/testSet_genericAll.txt	48
# - Generic Test Set
# data/trainSet_122586.26.txt	1
# - Neisseria_meningitidis_MC58.gb.gz
# data/trainSet_122587.18.txt	1
# - Neisseria_meningitidis_Z2491.gb.gz
# data/trainSet_1280.10152.txt	1


open(T, "$phispy -l long 2>&1 |") or die "Cannot list training datasets: $!";
my $cur;
while (<T>)
{
    chomp;
    if (/^(data\S+)/)
    {
	$cur = $1;
	next;
    }
    if (/^-\s+(.*)(\.gb\.gz)?$/)
    {
	my $n = $1;
	$n =~ s/_/ /g;

	# print "$cur: $n\n";
	if ($sci =~ /^$n\b/)
	{
	    $training_set = $cur;
	    last;
	}

    }
}
close(T);

my $out = File::Temp->newdir(undef, CLEANUP => 1);

my @cmd = ($phispy,
	   "-o", "" . $out,
	   (defined($training_set) ? ("-t", $training_set) : ()),
	   "$gb_file");
my $cmd = "@cmd > $out/phispy.stdout 2> $out/phispy.stderr";
print STDERR "Run $cmd\n";
$rc = system($cmd);
if ($rc != 0)
{
    my $err = read_file("$out/phispy.stderr");
    die "Error $rc running @cmd\n$err";
}

$id_prefix = $genomeTO->{id} unless $id_prefix;

my $type = 'prophage';
my $typed_prefix = join(".", $id_prefix, $type);

my $hostname = `hostname`;
chomp $hostname;

my $event = {
    tool_name => "phispy",
    execute_time => scalar gettimeofday,
    parameters => \@cmd,
    hostname => $hostname,
};
my $event_id = &GenomeTypeObject::add_analysis_event($genomeTO, $event);

# The columns of the file are:
#   - 1. Prophage number
#   - 2. The contig upon which the prophage resides
#   - 3. The start location of the prophage
#   - 4. The stop location of the prophage
# If we can detect the _att_ sites, the additional columns are:
#   - 11. start of _attL_;
#   - 12. end of _attL_;
#   - 13. start of _attR_;
#   - 14. end of _attR_;
#   - 15. sequence of _attL_;
#   - 16. sequence of _attR_;
#   - 17. The explanation of why this _att_ site was chosen for this prophage.

if (open(O, "<", "$out/prophage_coordinates.tsv"))
{
    while (<O>)
    {
	chomp;
	my($xid, $contig, $beg, $end) = split(/\t/);
	
	my $len = $end - $beg + 1;
	if ($contig)
	{
	    &GenomeTypeObject::add_feature($genomeTO, {
		-id_client => $id_client,
		-id_prefix => $id_prefix,
		-type       => $type,
		-location   => [[ $contig, $beg, '+', $len ]],
		-annotator  => 'phispy',
		-annotation => 'Add feature called by phispy',
		-analysis_event_id   => $event_id,
		-function => 'phiSpy-predicted prophage',
	    });
	}
    }
    close(O);
}
else
{
    my $err = read_file("$out/phispy.stderr");
    my $out = read_file("$out/phispy.stdout");
    print STDERR "Prophage run did not find any output. Stdout:\n$out\nStderr:\n$err\n";
}

$genomeTO->destroy_to_file($out_fh);
close($out_fh);
