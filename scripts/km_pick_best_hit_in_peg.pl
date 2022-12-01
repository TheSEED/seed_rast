use strict;
use Data::Dumper;

use Getopt::Long::Descriptive;
my $no_scores = 0;

my($opt, $usage) = describe_options("%c %o",
				    ["no-scores|n", "Do not write scores"],
				    ["help|h "=> "Show this help message"]);
print($usage->text), exit if $opt->help;
print($usage->text), exit 1 if (@ARGV != 0);

my $usage = "usage: km_pick_best_hit_in_peg [-n]\n";

use SeedUtils;
while (defined($_ = <STDIN>) && ($_ ne "---------------------\n")) { }
my $last = <STDIN>;
my %seen;
while ($last && ($last =~ /^(\S+)/))
{
    my $peg = $1;
    my %funcsW;
    my %funcs;
    my $z_sc = -100;
    #
    # Row format: 0 contig, 1 beg, 2 end, 3 strand, 4 frame, 5 hits, 6 func, 7 weighted-score, 8 funcIdx
    #
    #                    0     1    2    3    4    5      6            7       8     9
    #                  $1                         $2     $3            $4        $5 $6
    while (($last =~ /^(\S+)\t\S+\t\S+\t\S+\t\S+\t(\d+)\t(\S[^\t]*\S)\t(\S+)\t\S+(\t(\S+))?/) && ($peg eq $1))
    {
	my $hits = $2;
	my $func = $3;
	my $weighted_score = $4;
	my $hit_zsc = $6;
	
	$funcsW{$func} += $weighted_score;
	$funcs{$func} += $hits;
	# print Dumper($func, \%funcsW, \%funcs);
	$z_sc       = ($hit_zsc > $z_sc) ? $hit_zsc : $z_sc;
	$last = <STDIN>;
    }
    my @funcs = sort { $funcsW{$b} <=> $funcsW{$a} } keys(%funcs);
    my $best_func = $funcs[0];

    if ($opt->no_scores)
    {
	print "$peg\t$best_func\n";
    }
    else
    {
	print "$peg\t$best_func\t$funcs{$best_func}\t$funcsW{$best_func}\t$z_sc\n";
    }
}
