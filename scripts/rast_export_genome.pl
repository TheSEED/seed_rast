#!/usr/bin/env perl

use strict;
use warnings;
no warnings 'uninitialized';
use Data::Dumper;
use File::Copy;
use SeedUtils;
use gjoseqlib;
use Date::Parse;
use POSIX;
use DB_File;

use Bio::KBase::GenomeAnnotation::Client;
use JSON::XS;
use Digest::MD5 'md5_hex';

use GenomeTypeObject;
use URI::Escape;
use Bio::SeqIO;
use Bio::Seq;
use Bio::Seq::RichSeq;
use Bio::Species;
use Bio::Location::Split;
use Bio::Location::Simple;
use Bio::SeqFeature::Generic;
use Getopt::Long::Descriptive;
use Spreadsheet::Write;
use File::Temp;

my $temp_dir;

#
# Master list of supported formats. When changes happen here please propagate them to
# the rast2-export-genome script and to the documentation in the API spec.
#   

my @formats = ([genbank => "Genbank format"],
	       [genbank_merged => "Genbank format as single merged locus, suitable for Artemis"],
	       [spreadsheet_txt => "Spreadsheet (tab-separated text format)"],
	       [spreadsheet_xls => "Spreadsheet (Excel XLS format)"],
	       [feature_data => "Tabular form of feature data"],
	       [protein_fasta => "Protein translations in fasta format"],
	       [contig_fasta => "Contig DNA in fasta format"],
	       [feature_dna => "Feature DNA sequences in fasta format"],
	       [seed_dir => "SEED directory"],
	       [patric_features => "PATRIC features.tab format"],
	       [patric_specialty_genes => "PATRIC spgenes.tab format"],
	       [patric_genome_metadata => "PATRIC genome metadata format"],
	       [gff => "GFF format"],
	       [embl => "EMBL format"]);

my($opt, $usage) = describe_options("%c %o format",
				    ["feature-type=s@", 'Select a feature type to dump'],
				    ["input|i=s", "Input file"],
				    ["output|o=s", "Output file"],
				    ["genbank-roundtrip", "For features that have genbank types, output those"],
				    ["with-headings", "For downloads with optional headings (feature_data) include headings"],
				    ["specialty-gene-lookup-db=s", "Use this reference database for exporting spgene data"],
				    ["emit-contig-ids", "In generated data, use contig ids instead of accessions"],
				    ["help|h", "Show this help message"],
				    [],
				    ["Supported formats:\n"],
				    map { [ "$_->[0]  : $_->[1]" ] } @formats,
				   );

print($usage->text), exit if $opt->help;
print($usage->text), exit 1 if @ARGV != 1;

my $feature_type_ok;
if (ref($opt->feature_type))
{
    my $feature_type = { map { $_ => 1 } @{$opt->feature_type} };
    $feature_type_ok = sub {
	my($feat) = @_;
	return $feature_type->{$feat->{type}} ? 1 : 0;
    };
}
else
{
    $feature_type_ok = sub { 1 };
}

my $format = shift;

if (lc($format) eq 'gff')
{
    $format = 'GTF';
}

my $in_fh;
if ($opt->input) {
    open($in_fh, "<", $opt->input) or die "Cannot open " . $opt->input . ": $!";
} else { $in_fh = \*STDIN; }

my $out_fh;
if ($opt->output) {
    open($out_fh, ">", $opt->output) or die "Cannot open " . $opt->output . ": $!";
} else { $out_fh = \*STDOUT; }

my $genomeTO = GenomeTypeObject->create_from_file($in_fh);

#
# For each protein, if the function is blank, make it hypothetical protein.
for my $f (@{$genomeTO->{features}})
{
    next unless &$feature_type_ok($f);
    if (!defined($f->{function}) || $f->{function} eq '')
    {
	$f->{function} = "hypothetical protein";
    }
}

#
# Compute sequence ID to accession map.
#
my %seq_to_accession;
for my $ctg (@{$genomeTO->{contigs}})
{
    $seq_to_accession{$ctg->{id}} = id_for_contig($ctg);
}


#
# Simple exports.
#

if ($format eq 'protein_fasta')
{
    export_protein_fasta($genomeTO, $out_fh);
    exit;
}
elsif ($format eq 'contig_fasta')
{
    export_contig_fasta($genomeTO, $out_fh);
    exit;
}
elsif ($format eq 'patric_features')
{
    export_patric_features($genomeTO, $out_fh);
    exit;
}
elsif ($format eq 'patric_specialty_genes')
{
    export_patric_specialty_genes($genomeTO, $out_fh);
    exit;
}
elsif ($format eq 'patric_genome_metadata')
{
    export_patric_genome_metadata($genomeTO, $out_fh);
    exit;
}
elsif ($format eq 'feature_data')
{
    export_feature_data($genomeTO, $opt->with_headings, $out_fh);
    exit;
}
elsif ($format =~ /spreadsheet_(xls|txt)/)
{
    export_spreadsheet($genomeTO, $1, $out_fh);
    exit;
}
elsif ($format eq 'seed_dir')
{
    export_seed_dir($genomeTO, $out_fh);
    exit;
}
elsif ($format eq 'feature_dna')
{
    export_feature_dna($genomeTO, $out_fh);
    exit;
}

my $bio = {};
my $bio_list = [];

my $gs = $genomeTO->{scientific_name};
my $genome = $genomeTO->{id};

my $offset = 0;
my %contig_offset;

my $species;
my @tax;

# Noncompliant GTO
if (ref($genomeTO->{taxonomy}))
{
    @tax = @{$genomeTO->{taxonomy}};
}
else
{
    @tax = split(/;\s+/, $genomeTO->{taxonomy});
}

if (@tax == 0)
{
    #
    # If we have a NCBI lineage structure, use that for taxonomy.
    #
    my $lin = $genomeTO->{ncbi_lineage} // [];
    @tax = map { $_->[0] } @$lin;
}

#
# To avoid issues with EMBL exports on long taxonomies, truncate
# any field in the taxonomy to 74 chars.
#
if (lc($format) eq 'embl')
{
    for my $t (@tax)
    {
	if ((my $l = length($t)) > 73)
	{
	    
	    $t = substr($t, 0, 73);
	    print STDERR "EMBL export: truncating long taxonomy string (length=$l) to $t\n";
	}
    }
}

if (@tax)
{
    $species = Bio::Species->new(-classification => [reverse @tax]);
}
my $gc = $genomeTO->{genetic_code} || 11;

my $tax_id = $genomeTO->{ncbi_taxonomy_id};
if (!$tax_id)
{
    ($tax_id) = $genome =~ /^(\d+)\.\d+/;
}
my @taxon = defined($tax_id) ? (db_xref => "taxon:$tax_id") : ();

#
# code is similar but subtly different for the
# merged/non-merged cases.
#

if ($format eq 'genbank_merged')
{
    my $dna = '';

    my @feats;
    for my $c (@{$genomeTO->{contigs}})
    {
	my $contig = id_for_contig($c);
	$dna .= $c->{dna};
	
	my $contig_start = $offset + 1;
	my $contig_end   = $offset + length($c->{dna});
	
	$contig_offset{$contig} = $offset;
	
	my $feature = Bio::SeqFeature::Generic->new(-start => $contig_start,
						    -end => $contig_end,
						    -tag => {
							organism => $gs,
							@taxon,
							mol_type => "genomic DNA",
							note => $genome,
							note => $contig,
						    },
						    -primary => 'source');
	my $fa_record = Bio::SeqFeature::Generic->new(-start => $contig_start,
						      -end => $contig_end,
						      -tag => {
							  @taxon,
							  label => $contig,
							  note => $contig,
						      },
						      -primary => 'fasta_record');
	push(@feats, $feature, $fa_record);
	
	$offset += length($c->{dna});
    }
    my $bseq = Bio::Seq::RichSeq->new(-id => $genome,
				      -display_name => $gs,
				      -accession_number => id_for_contig($genomeTO->{contigs}->[0]),
				      (defined($species) ? (-species => $species) : ()),
				      -seq => $dna);

    for my $c (@{$genomeTO->{contigs}})
    {
	my $contig = id_for_contig($c);
	$bio->{$contig} = $bseq;
    }

    $bseq->add_SeqFeature($_) foreach @feats;
    @$bio_list = $bseq;
}
else
{
    for my $c (@{$genomeTO->{contigs}})
    {
	my $contig = id_for_contig($c);
	$bio->{$contig} = Bio::Seq::RichSeq->new(-id => $contig,
						 -display_name => $gs,
						 -accession_number => $contig,
						 (defined($species) ? (-species => $species) : ()),
						 -seq => $c->{dna});
	$bio->{$contig}->desc($gs);
	push(@$bio_list, $bio->{$contig});

	# print "Set div for $bio->{$contig}\n";
	# $bio->{$contig}->division("PHG");
	
	my $contig_start = $offset + 1;
	my $contig_end   = $offset + length($c->{dna});
	
	$contig_offset{$contig} = $offset;
	
	my $feature = Bio::SeqFeature::Generic->new(-start => $contig_start,
						    -end => $contig_end,
						    -tag => {
							  @taxon,
							organism => $gs,
							mol_type => "genomic DNA",
							note => $genome,
						    },
						    -primary => 'source');
	$bio->{$contig}->add_SeqFeature($feature);
	
	my $fa_record = Bio::SeqFeature::Generic->new(-start => $contig_start,
						      -end => $contig_end,
						      -tag => {
							  @taxon,
							  label => $contig,
							  note => $contig,
						      },
						      -primary => 'fasta_record');
	$bio->{$contig}->add_SeqFeature($fa_record);
	
	if ($format eq 'genbank_merged')
	{
	    $offset += length($c->{dna});
	}
    }
}
#
# Reset format to genbank since we've computed offsets.
#
$format = 'genbank' if $format eq 'genbank_merged';

my %protein_type = (CDS => 1, peg => 1, mat_peptide => 1);
my $strip_ec;
my $gff_export = [];

my $anno_source = $genomeTO->{source} // "RAST2";

my @features;
for my $f (@{$genomeTO->{features}})
{
    next unless &$feature_type_ok($f);
    my $peg = $f->{id};
    my $note = {};
    my $contig;

    my $func = "";
    if (defined($f->{function}) && $f->{function} ne '')
    {
	$func = $f->{function};
    }
    elsif ($protein_type{$f->{type}})
    {
	$func = "hypothetical protein";
    }

    push @{$note->{db_xref}}, "$anno_source:$peg";

    my %ecs;
    if ($func)
    {
	foreach my $role (SeedUtils::roles_of_function($func))
	{
	    my ($ec) = ($role =~ /\(EC (\d+\.\d+\.\d+\.\d+)\)/);
	    $ecs{$ec} = 1 if ($ec);
	}
	
	# add ECs
	push @{$note->{"EC_number"}}, keys(%ecs);
    }

    my($creation_event, $annotation_event, $annotation_tool) = $genomeTO->get_creation_info($f);

    $annotation_tool ||= $annotation_event->{tool_name};
    push(@{$note->{note}}, "rasttk_feature_creation_tool=$creation_event->{tool_name}") if $creation_event;
    push(@{$note->{note}}, "rasttk_feature_annotation_tool=$annotation_tool") if $annotation_tool;

    my $loc;

    my @loc_obj;
    my @loc_info;
    for my $l (@{$f->{location}})
    {
	my($ctg, $start, $strand, $len) = @$l;
	$contig = $ctg;
	my $offset = $contig_offset{$contig};
	my $end = $strand eq '+' ? ($start + $len - 1) : ($start - $len + 1);

	my $bstrand = 0;
	if ($protein_type{$f->{type}})
	{
	    $bstrand = ($strand eq '+') ? 1 : -1;
	}

	$start += $offset;
	$end += $offset;

	($start, $end) = ($end, $start) if $strand eq '-';

	my $sloc = new Bio::Location::Simple(-start => $start, -end => $end, -strand => $strand);
	push(@loc_obj, $sloc);

	#
	# Compute loc_info for GFF stuff.
	#
	my $frame = $start % 3;

	#
	# If we are mapping contig ids to accessions, we need to look up the accession here.
	#
	my $nctg = maybe_accession_of($ctg, $genomeTO);
	push(@loc_info, [$nctg, $start, $end, (($len == 0) ? "." : $strand), $frame]);
    }

    if (@loc_obj == 1)
    {
	$loc = $loc_obj[0];
    }
    elsif (@loc_obj > 1)
    {
	$loc = new Bio::Location::Split();
	$loc->add_sub_Location($_) foreach @loc_obj;
    }

    my $feature;
    my $source = "rast2_export";
    my $type = $f->{type};
	
    # strip EC from function
    my $func_ok = $func;
    
    if ($strip_ec) {
	$func_ok =~ s/\s+\(EC \d+\.(\d+|-)\.(\d+|-)\.(\d+|-)\)//g;
	$func_ok =~ s/\s+\(TC \d+\.(\d+|-)\.(\d+|-)\.(\d+|-)\)//g;
    }

    if ($protein_type{$f->{type}})
    {
	my $type = $opt->genbank_roundtrip ? ($f->{genbank_type} // $f->{type}) : 'CDS';
	$feature = Bio::SeqFeature::Generic->new(-location => $loc,
						 -primary  => $type,
						 -tag      => {
						     product     => $func_ok,
						     translation => $f->{protein_translation},
						 },
						);
	push @{$note->{transl_table}}, $gc;
	
	foreach my $tagtype (keys %$note) {
	    $feature->add_tag_value($tagtype, @{$note->{$tagtype}});
	}
	
	# work around to get annotations into gff
	# this is probably still wrong for split locations.
	$func_ok =~ s/ #.+//;
	$func_ok =~ s/;/%3B/g;
	$func_ok =~ s/,/%2C/g;
	$func_ok =~ s/=//g;
	for my $l (@loc_info)
	{
	    my $ec = "";
	    my @ecs = ($func =~ /[\(\[]*EC[\s:]?(\d+\.[\d-]+\.[\d-]+\.[\d-]+)[\)\]]*/ig);
	    if (scalar(@ecs)) {
		$ec = ";Ontology_term=".join(',', map { "KEGG_ENZYME:" . $_ } @ecs);
	    }
	    my($contig, $start, $stop, $strand, $frame) = @$l;
	    # for bacterial genomes phase=0
	    my $phase = 0;
	    push @$gff_export, "$contig\t$source\tCDS\t$start\t$stop\t.\t$strand\t$phase\tID=".$peg.";Name=".$func_ok.$ec."\n";
	}
    } elsif ($type eq "rna") {
	my $primary;
	if ( $func =~ /tRNA/ ) {
	    $primary = 'tRNA';
	} elsif ( $func =~ /(Ribosomal RNA|5S RNA|rRNA)/ ) {
	    $primary = 'rRNA';
	} else {
	    $primary = 'misc_RNA';
	}
	
	$feature = Bio::SeqFeature::Generic->new(-location => $loc,
						 -primary  => $primary,
						 -tag      => {
						     product => $func,
						 },
						 
						);
	$func_ok =~ s/ #.+//;
	$func_ok =~ s/;/%3B/g;
	$func_ok =~ s/,/%2C/g;
	$func_ok =~ s/=//g;
	foreach my $tagtype (keys %$note) {
	    $feature->add_tag_value($tagtype, @{$note->{$tagtype}});
	}
	# work around to get annotations into gff
	for my $l (@loc_info)
	{
	    my($contig, $start, $stop, $strand, $frame) = @$l;
	    push @$gff_export, "$contig\t$source\t$primary\t$start\t$stop\t.\t$strand\t.\tID=$peg;Name=$func_ok\n";
	}
	
    } elsif ($type eq "crispr_repeat" || $type eq 'repeat') {
	my $primary = "repeat_region";
	$feature = Bio::SeqFeature::Generic->new(-location => $loc,
						 -primary  => $primary,
						 -tag      => {
						     product => $func,
						 },
						 
						);
	$feature->add_tag_value("rpt_type",
				($type eq 'crispr_repeat' ? "direct" : "other"));
	$func_ok =~ s/ #.+//;
	$func_ok =~ s/;/%3B/g;
	$func_ok =~ s/,/%2C/g;
	$func_ok =~ s/=//g;
	foreach my $tagtype (keys %$note) {
	    $feature->add_tag_value($tagtype, @{$note->{$tagtype}});
	}
	# work around to get annotations into gff
	for my $l (@loc_info)
	{
	    my($contig, $start, $stop, $strand, $frame) = @$l;
	    push @$gff_export, "$contig\t$source\t$primary\t$start\t$stop\t.\t$strand\t.\tID=$peg;Name=$func_ok\n";
	}
	
    } else {
	my $primary = "misc_feature";
	$feature = Bio::SeqFeature::Generic->new(-location => $loc,
						 -primary  => $primary,
						 -tag      => {
						     product => $func,
						     note => $type,
						 },
						 
						);
	$func_ok =~ s/ #.+//;
	$func_ok =~ s/;/%3B/g;
	$func_ok =~ s/,/%2C/g;
	$func_ok =~ s/=//g;
	foreach my $tagtype (keys %$note) {
	    $feature->add_tag_value($tagtype, @{$note->{$tagtype}});
	}
	# work around to get annotations into gff
	for my $l (@loc_info)
	{
	    my($contig, $start, $stop, $strand, $frame) = @$l;
	    push @$gff_export, "$contig\t$source\t$primary\t$start\t$stop\t.\t$strand\t.\tID=$peg;Name=$func_ok\n";
	}
	

    }

    $contig = maybe_accession_of($contig);
    my $bc = $bio->{$contig};
    if (ref($bc))
    {
	$bc->add_SeqFeature($feature) if $feature;
    }
    else
    {
	print STDERR "No contig found for $contig on $feature\n";
    }
    # hack for gtf via bioperl push(@features, $feature) if $feature;
}

# check for FeatureIO or SeqIO
if ($format eq "GTF") {
    if (0)
    {
	my $fio = Bio::FeatureIO->newFh(-fh => $out_fh, -format => "GTF");
	foreach my $feature (@features) {
	    $fio->write_feature($feature);
	}
    }
    else
    {
	print $out_fh "##gff-version 3\n";
	foreach (@$gff_export) {
	    print $out_fh $_;
	}
    }
} else {
#	my $sio = Bio::SeqIO->new(-file => ">$filename", -format => $format);

    #
    # bioperl writes lowercase dna. We want uppercase for biophython happiness.
    #
#    my $sio = Bio::SeqIO->new(-file => ">tmpout", -format => $format);
    my $sio = Bio::SeqIO->new(-fh => $out_fh, -format => $format);

    foreach my $seq (@$bio_list)
    {
	$sio->write_seq($seq);
    }
}

sub export_protein_fasta
{
    my($genomeTO, $out_fh) = @_;

    for my $f (@{$genomeTO->{features}})
    {
	next unless &$feature_type_ok($f);
	my $peg = $f->{id};
	if ($f->{protein_translation})
	{
	    print_alignment_as_fasta($out_fh, [$peg,
					       "$f->{function} [$genomeTO->{scientific_name} | $genomeTO->{id}]",
					       $f->{protein_translation}]);
	}
    }
}

sub export_feature_dna
{
    my($genomeTO, $out_fh) = @_;

    for my $f (@{$genomeTO->{features}})
    {
	next unless &$feature_type_ok($f);
	my $id = $f->{id};
	print_alignment_as_fasta($out_fh, [$id,
					   "$f->{function} [$genomeTO->{scientific_name} | $genomeTO->{id}]",
					   $genomeTO->get_feature_dna($id)]);
    }
}

sub export_contig_fasta
{
    my($genomeTO, $out_fh) = @_;

    for my $c (@{$genomeTO->{contigs}})
    {
	my $contig = id_for_contig($c);
	print_alignment_as_fasta($out_fh, [$contig, undef, $c->{dna}]);
    }
}

sub export_feature_data
{
    my($genomeTO, $with_headings, $out_fh) = @_;

    if ($with_headings)
    {
	print $out_fh join("\t", qw(feature_id location type function aliases protein_md5)), "\n";
    }

    my $features = $genomeTO->sorted_features();
    foreach my $feature (@$features)
    {
	next unless &$feature_type_ok($feature);
	my $fid = $feature->{id};
	my $loc = join(",",map { my($contig,$beg,$strand,$len) = @$_; 
				 "$contig\_$beg$strand$len" 
			       } @{$feature->{location}});
	my $type = $feature->{type};
	my $func = $feature->{function};
	my $md5 = "";
	$md5 = md5_hex(uc($feature->{protein_translation})) if $feature->{protein_translation};

	my @aliases = $genomeTO->flattened_feature_aliases($feature);
	my $aliases = join(",", @aliases);

	print $out_fh join("\t", $fid,$loc,$type,$func,$aliases,$md5), "\n";
    }
}

sub export_patric_features
{
    my($genomeTO, $out_fh) = @_;

    my @headings = qw(genome_id genome_name accession
		      annotation feature_type patric_id
		      refseq_locus_tag alt_locus_tag uniprotkb_accession
		      start end strand na_length gene product
		      figfam_id plfam_id pgfam_id go ec pathway);
		      
    print $out_fh join("\t", @headings), "\n";

    my $features = $genomeTO->sorted_features();
    my %common_dat = (genome_id => $genomeTO->{id},
		      genome_name => $genomeTO->{scientific_name});
		      
    foreach my $feature (@$features)
    {
	my %dat = %common_dat;

	my $fid = $feature->{id};

	my($contig, $min, $max, $dir) = SeedUtils::boundaries_of(map { my($c,$s,$d,$l) = @$_; "${c}_$s$d$l" } @{$feature->{location}});

	$dat{accession} = $contig;
	$dat{annotation} = 'PATRIC';
	$dat{feature_type} = $feature->{type};
	$dat{patric_id} = $fid;
	$dat{product} = $feature->{function};
	$dat{na_length} = $max - $min + 1;
	$dat{strand} = $dir;
	if ($dir eq '+')
	{
	    $dat{start} = $min;
	    $dat{end} = $max;
	}
	else
	{
	    $dat{start} = $max;
	    $dat{end} = $min;
	}

	for my $fam (@{$feature->{family_assignments}})
	{
	    my($ftype, $id, $val, $version) = @$fam;
	    if ($ftype eq 'PGFAM')
	    {
		$dat{pgfam_id} = $id;
	    }
	    elsif ($ftype eq 'PLFAM')
	    {
		$dat{plfam_id} = $id;
	    }
	    if ($ftype eq 'FIGFAM')
	    {
		$dat{figfam_id} = $id;
	    }
	}

	for my $ap (@{$feature->{alias_pairs}})
	{
	    my($ak, $v) = @$ap;
	    if ($ak eq 'locus_tag')
	    {
		$dat{refseq_locus_tag} = $v;
	    }
	    elsif ($ak eq 'gene')
	    {
		$dat{gene} = $v;
	    }
	}
	
	print $out_fh join("\t", @dat{@headings}), "\n";
    }
}

sub export_patric_specialty_genes
{
    my($genomeTO, $out_fh) = @_;

    my %lookup;

    if ($opt->specialty_gene_lookup_db)
    {
	tie %lookup, "DB_File", $opt->specialty_gene_lookup_db, O_RDONLY, 0, $DB_HASH
	    or warn "Could not tie specialty gene lookup database " . $opt->specialty_gene_lookup_db . ": $!";
    }

    my @headings = qw(genome_name patric_id refseq_locus_tag alt_locus_tag gene product property source
		      evidence source_id organism function classification pmid query_coverage subject_coverage
		      identity e_value);
		      
    print $out_fh join("\t", @headings), "\n";

    my $features = $genomeTO->sorted_features();
    my %common_dat = (genome_id => $genomeTO->{id},
		      genome_name => $genomeTO->{scientific_name},
		      organism => $genomeTO->{scientific_name});
		      
    # typedef tuple <string source, string source_id,
    #   float query_coverage, float subject_coverage, float identity, float e_value>
    # similarity_association;


    foreach my $feature (@$features)
    {
	my $assoc = $feature->{similarity_associations};
	next unless ref($assoc) && @$assoc;

	for my $val (@$assoc)
	{
	    my %dat = %common_dat;
	    my($db, $dbid, $qry, $subj, $iden, $eval) = @$val;

	    $dbid =~ s/^$db\|//;

	    my $ltxt = $lookup{$db, $dbid};
	    my $ldat = {};
	    if ($ltxt)
	    {
		$ldat = decode_json($ltxt);
	    }

	    $dat{source} = $db;
	    $dat{evidence} = 'BLASTP';
	    $dat{source_id} = $dbid;
	    $dat{query_coverage} = $qry;
	    $dat{subject_coverage} = $subj;
	    $dat{identity} = $iden;
	    $dat{e_value} = $eval;

	    my $fid = $feature->{id};

	    $dat{patric_id} = $fid;
	    $dat{function} = $feature->{function};

	    for my $ap (@{$feature->{alias_pairs}})
	    {
		my($ak, $v) = @$ap;
		if ($ak eq 'locus_tag')
		{
		    $dat{refseq_locus_tag} = $v;
		}
		elsif ($ak eq 'gene')
		{
		    $dat{gene} = $v;
		}
	    }

	    $dat{gene} = $ldat->{gene_name} if $ldat->{gene_name};
	    $dat{property} = $ldat->{property};
	    $dat{classification} = $ldat->{classification};
	    $dat{pmid} = $ldat->{pmid};

	    print $out_fh join("\t", @dat{@headings}), "\n";
	}
    }
}

sub export_patric_genome_metadata
{
    my($genomeTO, $out_fh) = @_;

    my @headings = qw(genome_id genome_name organism_name taxon_id genome_status strain serovar biovar 
		      pathovar mlst other_typing culture_collection type_strain 
		      completion_date publication bioproject_accession biosample_accession 
		      assembly_accession genbank_accessions refseq_accessions 
		      sequencing_centers sequencing_status sequencing_platform sequencing_depth assembly_method 
		      chromosomes plasmids contigs sequences genome_length gc_content patric_cds brc1_cds refseq_cds 
		      isolation_site isolation_source isolation_comments collection_date isolation_country 
		      geographic_location latitude longitude altitude depth other_environmental 
		      host_name host_gender host_age host_health body_sample_site body_sample_subsite 
		      other_clinical antimicrobial_resistance antimicrobial_resistance_evidence gram_stain 
		      cell_shape motility sporulation temperature_range optimal_temperature salinity 
		      oxygen_requirement habitat disease comments additional_metadata source source_id);

    print $out_fh join("\t", @headings), "\n";

    my $dat = getGenomeInfo($genomeTO);

    $dat->{source} = $genomeTO->{source};
    $dat->{source_id} = $genomeTO->{source_id};

    for my $f (@{$genomeTO->{features}})
    {
	if ($f->{type} eq 'CDS')
	{
	    $dat->{patric_cds}++;
	}
    }
    
    print $out_fh join("\t", @$dat{@headings}), "\n";
}

sub export_seed_dir
{
    my($genomeTO, $out_fh) = @_;

    my $td = File::Temp::tempdir(CLEANUP => 1);
    my $dir = "$td/$genomeTO->{id}";
    mkdir($dir) or die "Cannot mkdir $dir: $@";
    $genomeTO->write_seed_dir($dir);

    my $fh;
    open($fh, "cd $td; tar czf - '$genomeTO->{id}' |") or die "cannot open tar: $!";
    copy($fh, $out_fh);
    close($fh);
}

sub export_spreadsheet
{
    my($genomeTO, $suffix, $out_fh) = @_;
    
    my $features = $genomeTO->{features};

    my @cols = qw(contig_id feature_id type location start stop strand
		  function aliases plfam pgfam figfam evidence_codes nucleotide_sequence aa_sequence);

    my $tmp;
    my $ss;
    
    if ($suffix eq 'xls')
    {	
	$tmp = File::Temp->new(SUFFIX => ".$suffix");
	close($tmp);
	
	my $sheetname = substr("Features in $genomeTO->{scientific_name}", 0, 31);

	$ss = Spreadsheet::Write->new(file => "$tmp",
				      format => 'xls',
				      sheet => $sheetname,
				      styles  => {
					  header => { font_weight => 'bold' },
				      });

	$ss->addrow(map { { content => $_, style => 'header' } } @cols);
    }
    else
    {
	print $out_fh join("\t", @cols), "\n";
    }

    foreach my $feature ($genomeTO->sorted_features())
    {
	next unless &$feature_type_ok($feature);
	my %dat;

	my $fid = $feature->{id};

	$dat{feature_id} = $fid;

	my($contig, $min, $max, $dir) = SeedUtils::boundaries_of(map { my($c,$s,$d,$l) = @$_; "${c}_$s$d$l" } @{$feature->{location}});
	if (!$contig)
	{
	    die Dumper($feature);
	}
	$dat{contig_id} = $contig;
	($dat{start}, $dat{stop}) = ($dir eq '+') ? ($min, $max) : ($max, $min);
	$dat{strand} = $dir;
	
	$dat{location} = join(",",map { my($contig,$beg,$strand,$len) = @$_; 
				 "$contig\_$beg$strand$len" 
			       } @{$feature->{location}});

	$dat{type} = $feature->{type};
	$dat{function} = $feature->{function};

	$dat{aa_sequence} = $feature->{protein_translation} ? $feature->{protein_translation} : '';
	$dat{nucleotide_sequence} = $genomeTO->get_feature_dna($fid);

	$dat{evidence_codes} = '';

	for my $fam (@{$feature->{family_assignments}})
	{
	    my($ftype, $id, $val, $version) = @$fam;
	    if ($ftype eq 'PGFAM')
	    {
		$dat{pgfam} = $id;
	    }
	    elsif ($ftype eq 'PLFAM')
	    {
		$dat{plfam} = $id;
	    }
	    if ($ftype eq 'FIGFAM')
	    {
		$dat{figfam} = $id;
	    }
	}


	$dat{aliases} = join(",", $genomeTO->flattened_feature_aliases($feature));

	if ($ss)
	{
	    $ss->addrow(@dat{@cols});
	}
	else
	{
	    print $out_fh join("\t", @dat{@cols}), "\n";
	}
    }

    if ($ss)
    {
	undef $ss;
	copy("$tmp", $out_fh);
    }
}

#
# Hoisted from rast2solr; compute aggregate genome metadata.
#
sub getGenomeInfo
{
    my($genomeObj) = @_;
    my $genome = {};
    
    my ($chromosomes, $plasmids, $contigs, $sequences, $cds, $genome_length, $gc_count, $taxon_lineage_ranks);
    
    $genome->{owner} = $genomeObj->{owner}? $genomeObj->{owner} : "PATRIC";
    
    $genome->{genome_id} = $genomeObj->{id};
    $genome->{genome_name} = $genomeObj->{scientific_name};
    $genome->{common_name} = $genomeObj->{scientific_name};
    $genome->{common_name}=~s/\W+/_/g;
    $genome->{organism_name} = $genome->{common_name};
    
    $genome->{taxon_id}    =  $genomeObj->{ncbi_taxonomy_id};
    
    foreach my $type (@{$genomeObj->{typing}}){
	$genome->{mlst} .= "," if $genome->{mlst};
	$genome->{mlst} .= $type->{typing_method}.".".$type->{database}.".".$type->{tag};
	#$genome->{mlst} .= "," if $genome->{mlst};
	#$genome->{mlst} .= @{$type}[0].".".@{$type}[1].".".@{$type}[2];
    }
    
    foreach my $seqObj (@{$genomeObj->{contigs}})
    {
	if ($seqObj->{genbank_locus}->{definition}=~/chromosome|complete genome/i)
	{
	    $chromosomes++;
	} elsif ($seqObj->{genbank_locus}->{definition}=~/plasmid/i)
	{
	    $plasmids++;
	} else
	{
	    $contigs++;
	}
	
	$sequences++;
	$genome_length += length($seqObj->{dna});
	$gc_count += $seqObj->{dna}=~tr/GCgc//;
	
	if ($sequences == 1){
	    foreach my $dblink (@{$seqObj->{genbank_locus}->{dblink}}){
		$genome->{bioproject_accession} = $1 if $dblink=~/BioProject:\s*(.*)/;
		$genome->{biosample_accession} = $1 if $dblink=~/BioSample:\s*(.*)/;
	    }
	    foreach my $reference (@{$seqObj->{genbank_locus}->{references}}){
		$genome->{publication} .= $reference->{PUBMED}."," unless $genome->{publication}=~/$reference->{PUBMED}/;
	    }
	    $genome->{publication}=~s/,*$//g;
	    $genome->{completion_date} = strftime "%Y-%m-%dT%H:%M:%SZ", localtime str2time($seqObj->{genbank_locus}->{date}) if $seqObj->{genbank_locus}->{date};
	}
	
	if ($seqObj->{genbank_locus}->{accession}[0]=~/^([A-Z]{4})\d{8}$/){ # wgs, capture only master accession
	    $genome->{genbank_accessions} .= $1."00000000" unless $genome->{genbank_accessions}=~/$1.0000000/;
	}else{
	    $genome->{genbank_accessions} .= $seqObj->{genbank_locus}->{accession}[0].",";
	}
	
    }
    $genome->{genbank_accessions}=~s/,*$//g;
    
    $genome->{chromosomes} = $chromosomes if $chromosomes;
    $genome->{plasmids} = $plasmids if $plasmids;
    $genome->{contigs} = $contigs if $contigs ;
    $genome->{sequences} = $sequences;
    $genome->{genome_length} = $genome_length;
    $genome->{gc_content} = sprintf("%.2f", ($gc_count*100/$genome_length));
    $genome->{genome_status} = ($contigs > 1)? "WGS": "Complete";

    return $genome;
}

#
# Use the --emit-contig-ids parameter to determine if we
# use the sequence_id (aka id field) or accession (aka original_id)
# field (if it exists) when emitting contig identifiers.
#
sub id_for_contig
{
    my($ctg) = @_;
    if ($opt->emit_contig_ids)
    {
	return $ctg->{id};
    }
    else
    {
	return $ctg->{original_id} // $ctg->{id};
    }
}

#
# In the case where we are emitting accessions, need
# to map the given sequence ID to an accession.
#
sub maybe_accession_of
{
    my($id, $gto) = @_;
    return $seq_to_accession{$id};
}
