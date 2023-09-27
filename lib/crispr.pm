package crispr;

# This is a SAS component.

=head1 Find CRISPR repeat arrays in genomic DNA

All positions in this code are 0-based numbering from the start of the
sequence.  Adjustment is made when the output report is generated.

=head2 Revision notes

=head3 2023-09-08 - 13

Nearly complete rewrite of all key steps.  There are many tests of
intermediate results that screen out repeated sequences in DNA that
are not CRISPR arrays.

=cut

use strict;
use gjoseqlib;
use gjostat;
use Data::Dumper;
# use Time::HiRes qw(time);

require Exporter;
our @ISA    = qw( Exporter );
our @EXPORT = qw( find_crisprs );

=head2 find_crisprs()

Find CRISPR repeat regions in DNA sequences.

=head3 Usage

     @crisprs = find_crisprs( \@seq_entries, \%opts )
     @crisprs = find_crisprs(  $seq_entry,   \%opts )
    \@crisprs = find_crisprs( \@seq_entries, \%opts )
    \@crisprs = find_crisprs(  $seq_entry,   \%opts )

=head3 Output

     @crisprs = ( [ $array_loc, $repeat_consen, \@repeats, \@spacers ], ... )

     $array_loc = [ [ $contig, $beg, $dir, $len ] ]
     @repeats   = ( [ $rep_loc, $rep_seq ], ... )
     $rep_loc   = [ [ $contig, $beg, $dir, $len ] ]
     @spacers   = ( [ $spc_loc, $spc_seq ], ... )
     $spc_loc   = [ [ $contig, $beg, $dir, $len ] ]

=over 4

All locations are SEEDtk style feature locations.  All locations produced
are one continuous stretch of the contig, they are a list of one
component.  The string form of $loc would be:

=back

        join( '', $loc->[0]->[0], '_', @{$loc->[0]}[1,2,3] )

=head3 Options

    debug      => bool      #  Debugging information about candidate array
                            #      rejection tests.  Nearly all of these are
                            #      commented out due to the shear volume.
                            #      Search code for $debug.
    maxreplen  => int       #  Maximum repeat length (D = 55)
    maxsplen   => int       #  Maximum spacer length (D = 50)
    minmatchid => fraction  #  Minimum conservation to match repeat; default
                            #     is defined by pmatch
    minpartial => int       #  Array may start/end with a partial repeat this
                            #     short (D = 8)
    minrepeat  => int       #  Minimum nuber of repeats in an array (D = 3)
    minreplen  => int       #  Minimum repeat length (D = 24)
    minsplen   => int       #  Minimum spacer length (D = 30)
    mintimes   => int       #  Minimum nuber of repeats in an array (D = 4)
    pmatch     => fraction  #  P-value threshold for adding another repeat
                            #     (D = 0.001); this value alters minmatchid
    spacerdev  => int       #  Acceptable deviation of spacer length from median (D = 5)
    spclendev  => int       #  Acceptable deviation of spacer length from median (D = 5)
    start_cid  => id        #  Ignore all contigs before this contig id;
                            #     allows continuing an interrupted analysis.

=head3 Analysis Steps

The subroutine finds CRISPR arrays through the following series of steps and tests:

=over 4

=item 1

Find a minimal candidate array with minrepeat precisely matching
repeats of exactly minreplen nucleotides.

=item 2

Check that the candidate repeat sequence is has compositional
complexity.  If not, continue with step 1 from the end of
this candidate array.

=item 3

Try to get a longer array by progressively shifting the start
site 1 nt to the left.  If the array gets longer, record the
start location and new length.  Repeat this step until the
resulting candidate array gets shorter.

=item 4

Try extending the array to the left and the right by allowing
imperfect matching to the repeat.  In each step (each new
repeat), the trial repeat sequence is updated to the previous
matching sequence, allowing the repeat to walk.  These are
still minimum length repeats.

=item 5

Verify the the number of repeats is at least minrepeat (why this
should not always be so is not clear, but test does sometimes
fail.  If it fails, restart step 1 from the end of this
candidate array.

=item 6

This is the final set of full-length repeat regions.

=item 7

Calculate the median repeat period of the array, and that it falls
within the acceptable range defined by min and max repeat and
spacer lengths.  If not, restart step 1 from the end of this
candidate array.

=item 8

Count the number of repeat periods that deviate from the median
length by more than spclendev (spacer length deviation).  If
more that 1/4 of spacers exceed deviation, restart step 1 from
the end of this candidate array.

=item 9

To this point, all analyses have been based on minimum length
repeats.  Conceptually, these define an alignment, which might
be extended to the left and/or the right.  A candidate new column
is added if its nucleotide would change identity less once per
8 repeat regions.  This expands the repeat region to its final
length.

=item 10

If the resulting repeat length is greater than maxreplen, restart
step 1 from the end of this candidate array.

=item 11

If the resulting median spacer length is less than minsp, restart
step 1 from the end of this candidate array.

=item 12

Extract the repeat sequences.

=item 13

Calculate the consensus nucleotide at each position.

=item 14

Look for a partial repeat to the left of the first repeat, and to
the right of the last repeat.

=item 15

Sanity check the final set of repeat coordinates.  All failures at
this step are again, repeated sequences in the genome, not bad
analysis of a CRIPR array.  On failed check, restart step 1 from
the end of this candidate array.

=item 16

Construct the description of this CRISPR array.  This is the first
time that the spacer sequences are extracted.

=item 17

Examine the spacer sequences for repeated sequence content.  Compile
the length 6 kmers in the spacers.  If more than 5 of the
kmers occur more than 0.6 times the number of spacers, reject
the array as being repeated sequences.  Restart step 1 from
the end of this candidate array.

=item 18

Push the CRISPR array data on the output list of arrays found.

=back

=cut

sub find_crisprs
{
    my ( $seqs, $opts ) = @_;
    $opts ||= {};

    my $debug      = defined $opts->{ debug } ? $opts->{ debug } : 0;

    my $minrepeat  = $opts->{ minrepeat  } ||
                     $opts->{ mintimes   } ||  3;      #  Min number of repeats

    my $minreplen  = $opts->{ minreplen  } ||
                     $opts->{ repeatlen  } || 24;      #  Min repeat length
    my $maxreplen  = $opts->{ maxreplen  } || 55;      #  Max repeat length

    my $minsp      = $opts->{ minsplen   } || 25;      #  Min spacer length
    my $maxsp      = $opts->{ maxsplen   } || 54;      #  Max spacer length

    my $minpartial = $opts->{ minpartial } ||  8;      #  An array may start/end with a repeat this long

    my $min_mat_id = $opts->{ minmatchid } ||  0;      #  Minimum match identity (D = P-val based)
    my $p_match    = $opts->{ pmatch     } ||  0.0001; #  Maximum match P-value

    my $spclendev  = $opts->{ spacerdev  } ||
                     $opts->{ spclendev  } ||  5;      #  Acceptable deviation of space length from median

    my $start_cid  = $opts->{ start_cid  } || '';      #  Ignore all contig ids before this;
                                                       #      allows continuing an analysis

    #  Failure to make these ints can create position indexing issues.

    $minrepeat  = int( $minrepeat );
    $minreplen  = int( $minreplen );
    $maxreplen  = int( $maxreplen );
    $minsp      = int( $minsp );
    $maxsp      = int( $maxsp );
    $minpartial = int( $minpartial );
    $spclendev  = int( $spclendev );

    # $start_cid = '326424.13:NC_008278';  # $debug

    #  This parameter is the maximum space between the end of a candidate
    #  repeat sequence and the start of the next.  It is used in the regular
    #  expression that identifies candidate repeat regions.

    my $maxgap = $maxsp + ( $maxreplen - $minreplen );

    #  Before the precise boundaries between repeats and spacers is defined,
    #  we will test the periodicity of the array against these:

    my $minperiod = $minsp + $minreplen;
    my $maxperiod = $maxsp + $maxreplen;

    #  If a contig is not at least this long, it cannot have an array.

    my $min_alen  = $minrepeat * $minperiod - $minsp;

    my $min_nid = $min_mat_id ? int( $min_mat_id * $minreplen )
                              : gjostat::binomial_critical_value_m_ge( $minreplen, 0.3, $p_match );

    my @crisprs;
    foreach my $entry ( ref $seqs->[0] eq 'ARRAY' ? @$seqs : ( $seqs ) )
    {
        my $cid = $entry->[0];
        if ( $start_cid )
        {
            next  if $cid ne $start_cid;
            $start_cid = '';
        }
        my $seq    = lc $entry->[2];
        my $seqlen = length $seq;
        if ( $seqlen < $min_alen )
        {
            # print STDERR "Contig $cid: too short ($seqlen nt) to contain a valid array.\n"  if $debug;
            next;
        }

        #  Find each repeat array:

        #  The pattern minimizes spacer length within the defined range.  This
        #  will avoid missing repeats when the repeat length and spacer lengths
        #  are both short.

        my $minrep_m1 = $minrepeat - 1;  #  Repeat count in the regex must be a value, not an expression
        while ( $seq =~ m/([acgt]{$minreplen})(?:[acgt]{$minsp,$maxgap}?\1){$minrep_m1,}/g )
        {
            my $n1     = $-[0];       #  First pos of match
            my $n2     = $+[0];       #  First pos beyond match (the whole array)
            my $repseq = $1;          #  Candidate CRISPR repeat sequence
            my $alen   = $n2 - $n1;   #  Array length

            #  Low complexity test on repeat sequence.  This is very loose.

            my $na = $repseq =~ tr/a//;
            my $nc = $repseq =~ tr/c//;
            my $ng = $repseq =~ tr/g//;
            my $nt = $repseq =~ tr/t//;
            my ( $maxcnt ) = sort { $b <=> $a } ( $na, $nc, $ng, $nt );
            if ( $maxcnt > 0.6 * length($repseq) )
            {
                # print STDERR "Contig $cid, at $n1: low complexity repeat skipped: $repseq\n"  if $debug;
                pos($seq) = $n2 + 1;  # restart following failed repeat array
                next;
            }

            #  If we shift the repeat start position to the right, can we get
            #  a longer array of perfect matches?

            my $max_alen = $alen;
            my $best_n1  = $n1;
            my $trial_n1 = $n1 + 1;

            #  Stopping conditions are failure of regex match, or shorter
            #  array length.

            while ( 1 )
            {
                my $trial_rep = substr($seq,$trial_n1,$minreplen);
                #  Anchor search at  trial repeat position
                pos($seq) = $trial_n1;
                last  unless $seq =~ m/\G$trial_rep(?:[acgt]{$minsp,$maxgap}?$trial_rep){$minrep_m1,}/;

                my $trial_alen = $+[0] - $-[0];     #  Length of this array
                last  if $trial_alen < $max_alen;   #  Array got shorter, we're done

                if ( $trial_alen > $max_alen )
                {
                    $best_n1  = $trial_n1;
                    $max_alen = $trial_alen;
                }
                $trial_n1++;
            }

            #  Was a better framing of the minimum length repeat was found?

            if ( $max_alen > $alen )
            {
                # print STDERR "Contig $cid: start moved from $n1 to $best_n1 lengthening array from $alen to $max_alen nt.\n"  if $debug;
                $n1     = $best_n1;
                $n2     = $best_n1 + $max_alen;
                $repseq = substr($seq,$n1,$minreplen);
                $alen   = $max_alen;
            }

            #  We have at least $minrepeat perfect repeats of length $minreplen.
            #  Try to lengthen the array with imperfect matching.

            my $maxdiff  = 3;        #  max differences between consecutive repeats
            my $maxreach = 2 * ( $maxsp + $maxreplen );  #  Allow one missing repeat
            my @locs = ( extend_array_backward( \$seq, $repseq, $n1, $maxdiff, $maxreach ),
                         $n1,
                         extend_array_forward( \$seq, $repseq, $n1, $maxdiff, $maxreach )
                       );

            #  It is not clear how there could be too few repeats, but check.

            if ( @locs < $minrepeat )
            {
                # print STDERR "Contig $cid, at $n1: too few repeast in array\n"  if $debug;
                pos($seq) = $n2 + 1;  # restart following failed repeat array
                next;
            }

            #  We need to screen candidate repeat locations for spacing issues.
            #  This seems to be mostly due to repeated sequences in the DNA.

            my $prev_loc = 0;
            my ( undef, @periods ) = map { my $dist = $_ - $prev_loc; $prev_loc = $_; $dist }
                                     @locs;

            #  The median period of repeats.  Note that this can be fractional,
            #  so converting to an int is essential.
            my $medperiod = int( gjostat::median( @periods ) );
            my $n_bad = grep { abs($_ - $medperiod) > $spclendev } @periods;

            if ( $medperiod < $minperiod )
            {
                # print STDERR "Contig $cid, at $n1: median repeat period $medperiod is smaller than the minimum $minperiod.\n"  if $debug;
                pos($seq) = $n2 + 1;  # restart following failed repeat array
                next;
            }

            if ( $medperiod > $maxperiod )
            {
                # print STDERR "Contig $cid, at $n1: median repeat period $medperiod is greater than the maximum $maxperiod.\n"  if $debug;
                pos($seq) = $n2 + 1;  # restart following failed repeat array
                next;
            }

            if ( $n_bad > max(0.25*@periods, 2) )  #  The larger of 1/4 of spacers, or 2
            {
                #  Failure of this test is almost always due to repeated sequences in spacers
                if ( 0 && $debug )
                {
                    print STDERR "Contig $cid, at $n1: too many bad repeat spacings.\n";
                    for ( my $i = 0; $i < @periods; $i++ )
                    {
                        print STDERR substr($seq,$locs[$i],$periods[$i]), "\n";
                    }
                }
                pos($seq) = $n2 + 1;  # restart following failed repeat array
                next;
            }

            #  Extend repeat to the left and right

            #  This might be too stringent, one change per 8 sequences
            my $max_flip = int( @locs / 8 );

            my $beg_adj = 0;
            my $max_i = $maxreplen - $minreplen;
            for ( my $i = 1; $i <= $max_i; $i++ )
            {
                my $n_flip = 0;
                my $nt0;
                foreach ( @locs )
                {
                    my $pos = $_ - $i;
                    last  if $pos < 0;
                    my $nt = substr($seq,$pos,1);
                    if ( $nt0 && $nt ne $nt0 )
                    {
                        $n_flip++;
                        last  if $n_flip > $max_flip;
                    }
                    $nt0 = $nt;
                }
                last  if $n_flip > $max_flip;

                $beg_adj = $i;
            }

            my $replen = $minreplen + $beg_adj;
            for ( my $i = $minreplen; $i <= $maxreplen; $i++ )
            {
                my $n_flip = 0;
                my $nt0;
                foreach ( @locs )
                {
                    my $pos = $_ + $i;
                    last  if $pos >= $seqlen;
                    my $nt = substr($seq,$pos,1);
                    if ( $nt0 && $nt ne $nt0 )
                    {
                        $n_flip++;
                        last  if $n_flip > $max_flip;
                    }
                    $nt0 = $nt;
                }
                last  if $n_flip > $max_flip;

                $replen++;
            }

            if ( $replen > $maxreplen )
            {
                # print STDERR "Contig $cid, at $n1: repeat $repseq extended to length $replen, exceeding maximum of $maxreplen.\n"  if $debug;
                pos($seq) = $n2 + 1;  # restart following the failed array
                next;
            }

            if ( $medperiod < $replen + $minsp )
            {
                my $medsp = $medperiod - $replen;
                # print STDERR "Contig $cid, at $n1: median spacer length of $medsp is shorter than the minimum of $minsp.\n"  if $debug;
                pos($seq) = $n2 + 1;  # restart following the failed array
                next;
            }

            #  Process the locations, and count base frequencies:

            my @rep_data;  # [ $beg, $end, $seq ]; still 0-based numbering
            my @counts;
            foreach ( @locs )
            {
                my $beg = $_ - $beg_adj;
                my $end = $beg + $replen - 1;
                if ( $beg < 0 )
                {
                    $beg = 0;
                    my $repseq = substr($seq,$beg,$end-$beg+1);
                    push @rep_data, [ $beg, $end, $repseq ];
                }
                elsif ( $end >= $seqlen )
                {
                    $end = $seqlen - 1;
                    my $repseq = substr($seq,$beg,$end-$beg+1);
                    push @rep_data, [ $beg, $end, $repseq ];
                }
                else    #  full length repeat
                {
                    my $rep = substr($seq,$beg,$replen);
                    push @rep_data, [ $beg, $end, $rep ];
                    $rep =~ tr/acgt/\000/c;
                    $rep =~ tr/acgt/\001\002\003\004/;
                    for ( my $i = 0; $i < $replen; $i++ )
                    {
                        $counts[$i]->[ord substr($rep,$i,1)]++;
                    }
                }
            }

            #  Repeat consensus sequence

            my $consen = '';
            foreach my $pos_cnt ( @counts )
            {
                 my ( $nt ) = map  { $_->[0] }
                              sort { $b->[1] <=> $a->[1] }
                              ( [ a => $pos_cnt->[1] ],
                                [ c => $pos_cnt->[2] ],
                                [ g => $pos_cnt->[3] ],
                                [ t => $pos_cnt->[4] ]
                              );
                $consen .= $nt;
            }

            #  Are there partial repeats before or after these?  Demand at least
            #  $minpartial perfectly matching nt, at an appropriate spacing.  With
            #  $minpartial = 8, there is about a 1/5000 false positive rate.

            my $medspacer = $medperiod - $replen;      #  median spacer length
            my $minspacer = $medspacer - $spclendev;
            my $maxspacer = $medspacer + $spclendev;

            my ( $r_beg, $r_end, $r_seq ) = @{$rep_data[0]};
            if ( length($r_seq) == $replen && $r_beg > $minspacer + $minpartial )
            {
                # r_beg-(maxspacer+minpartial)
                #      |                                   r_beg
                #      |     |<--------- maxspacer --------->|
                #      |     |          |<--- minspacer ---->|
                #   ----------------------------------------------------------
                #   ...======                                =============...
                #              ...======                     first full repeat
                #      partial repeat
                #        end range
                #
                my $pos1 = $r_beg - ( $maxspacer + $minpartial );
                my $len  = $minpartial + 2 * $spclendev;
                if ( $pos1 < 0 )
                {
                    $len += $pos1;
                    $pos1 = 0;
                }
                my $probe = substr($r_seq,-$minpartial);   # last $minpartial nt
                if ( substr($seq,$pos1,$len) =~ /$probe/ )
                {
                    my $b = $-[0] + $pos1;       # index into $seq
                    my $e = $+[0] + $pos1;       # index into $seq
                    my $i = $replen - $minpartial;   # index into $r_seq
                    while ( $b > 0 && $i > 0 && substr($seq,$b-1,1) eq substr($r_seq,$i-1,1) )
                    {
                        $b--;
                        $i--;
                    }
                    unshift @rep_data, [ $b, $e-1, substr($r_seq,$i) ];
                }
            }

            ( $r_beg, $r_end, $r_seq ) = @{$rep_data[-1]};
            if ( length($r_seq) == $replen && $r_end + $minspacer + $minpartial < $seqlen )
            {
                #                                r_end+minspacer   r_end+maxspacer+minpartial
                #                r_end                  |                |
                #                  |<--------- maxspacer --------->|     |
                #                  |<--- minspacer ---->|          |     |
                #   ----------------------------------------------------------
                #   ...=============                     ======....
                #   last full repeat                                ======...
                #                                        partial repeat
                #                                        start range
                #
                my $pos1 = $r_end + $minspacer + 1;
                my $len  = $minpartial + 2 * $spclendev;
                if ( $pos1 + $len > $seqlen )
                {
                    $len  = $seqlen - $pos1;
                }
                my $probe = substr($r_seq,0,$minpartial);   # firat $minpartial nt
                if ( substr($seq,$pos1,$len) =~ /$probe/ )
                {
                    my $b = $-[0] + $pos1;       # index into $seq
                    my $e = $+[0] + $pos1 - 1;   # index into $seq
                    my $i = $minpartial - 1;         # index into $r_seq
                    while ( $e+1 < $seqlen && $i+1 < $replen && substr($seq,$e+1,1) eq substr($r_seq,$i+1,1) )
                    {
                        $e++;
                        $i++;
                    }
                    push @rep_data, [ $b, $e, substr($r_seq,0,$i+1) ];
                }
            }

            #  Sanity check on repeat data (this could be much more stringent)
            #  All identified cases of this happening are repeated sequences
            #  in DNA, not bad CRISPR elements.  They would have failed the
            #  spacer kmer frequency test, but they have invalid coordinates
            #  that mess up the extraction of spacer sequences.

            my $prev_end = -1;
            my $bad_locs;
            foreach ( @rep_data )
            {
                my ( $r_beg, $r_end, $r_seq ) = @$_;
                if ( $r_beg <= $prev_end
                  || $r_end <= $r_beg
                  || $r_end >= $seqlen
                  || $r_end - $r_beg + 1 != length($r_seq)
                   )
                {
                    $bad_locs = 1;
                    if ( 0 && $debug )
                    {
                        print STDERR "Contig $cid, at $n1: bad repeat data: previous end=$prev_end, beg=$r_beg, end=$r_end, seqlen=$seqlen, repeat=$r_seq\n";
                        for ( my $i = 0; $i < @rep_data; $i++ )
                        {
                            my ( $b, $e, $rep ) = @{$rep_data[$i]};
                            my $next_b = $i+1 < @rep_data ? $rep_data[$i+1]->[0] : -1;
                            my $sp = $next_b > $e ? substr($seq,$e+1,$next_b-$e) : '';
                            print STDERR join( "\t", $b, $e, $rep, $sp ? $sp : () ), "\n";
                        }
                    }
                    last;
                }
                $prev_end = $r_end;
            }
            if ( $bad_locs )
            {
                #  restart following failed array
                pos($seq) = $n2 + 1;
                next;
            }

            #  Describe the array:

            # @rep_data = ( [ $beg, $end, $seq ], ...); still 0-based numbering
            my $a_beg = $rep_data[ 0]->[0];
            my $a_end = $rep_data[-1]->[1];
            my $a_len = $a_end - $a_beg + 1;
            my $a_loc = [ [ $cid, $a_beg+1, '+', $a_len ] ];  #  Array location

            #  Describe repeats and spacers:

            my @reps;
            my @spcs;
            $prev_end = 0;
            foreach ( @rep_data )
            {
                # Repeat element:

                my ( $b, $e, $s ) = @$_;
                my $len = $e - $b + 1;
                #  [ $rep_location, $rep_sequence ]
                push @reps, [ [ [ $cid, $b+1, '+', $len ] ], $s ];

                # Spacer element before this repeat?

                if ( $prev_end )
                {
                    my $s_beg = $prev_end + 1;
                    my $s_len = $b - $s_beg;
                    if ( $s_len > 0 )
                    {
                        #  [ $sp_location, $sp_sequence ]
                        push @spcs, [ [ [ $cid, $s_beg+1, '+', $s_len ] ], substr($seq,$s_beg,$s_len) ];
                    }
                }

                $prev_end = $e;
            }

            #  The spacers are still coming through with repeated sequences.
            #  For the time being, let's do a clean up based on kmer counts.

            my %kmers;
            my $klen = 6;
            foreach my $sp ( map { $_->[1] } @spcs )
            {
                my %k = map { $_ => 1 }                           # hash
                        map { m/(.{$klen})/g }                    # kmers
                        map { substr($sp,$_) } ( 0 .. $klen-1 );  # frames
                foreach ( keys %k ) { $kmers{$_}++ }
            }

            #  Number of kmers in 60% or more of the spacers.  Threshold set
            #  here to catch presence in 2 of 3 spacer sequences.

            my $nmax   = 0.6 * @spcs;
            my $nrecur = grep { $_ >= $nmax } values %kmers;

            if ( $nrecur >= 5 )
            {
                if ( 0 && $debug )
                {
                    print STDERR "Contig $cid, at $n1: repeated sequences in the spacers.\n";
                    print STDERR map { "    $_->[1]\n" } @spcs;
                }
                #  restart following failed array
                pos($seq) = $n2 + 1;
                next;
            }

            #  Save location, repeat consensus, repeats and spacers:
            #  $crispr = [ $loc, $repseq, \@repeats, \@spacers ];

            push @crisprs, [ $a_loc, $consen, \@reps, \@spcs ];

            #  Move past this repeat array:

            pos($seq) = $a_end + 1;

        }
    }

    wantarray ? @crisprs : \@crisprs;
}


=head2 Extend the repeat array

Find locations that match a given repeat sequence, working from the position of
that repeat.  It updates the repeat sequence with each new match, allowing
the sequence to creep.  It does not do a partial match at the array ends (that
is done later in the analysis).

    @locs = extend_array_backward( \$seq, $rept, $pos, $maxdiff, $maxreach )
    @locs = extend_array_forward(  \$seq, $rept, $pos, $maxdiff, $maxreach )

        \$seq       reference to DNA sequence (aka, contig sequence)
         $rept      initial repeat sequence
         $pos       starting point in $seq of $rept
         $maxdiff   maximum differences between previous repeat and new repeat
         $maxreach  maximum distance between match sites

         @locs      offsets into $seq for newly identified repeats

=cut

sub extend_array_backward
{
    my ( $seqR, $rept, $pos, $maxdiff, $maxreach ) = @_;
    my $replen = length( $rept );
    my @locs = ();
    while ( 1 )
    {
        my $match = 0;
        for ( my $i = $replen; $i <= $maxreach; $i++ )
        {
            my $p = $pos - $i;  #  trial position in $seq
            last if ( $p < 0 );

            my $trial_seq = substr($$seqR,$p,$replen);
            if ( n_diff( $rept, $trial_seq ) <= $maxdiff )
            {
                unshift @locs, $p;
                $rept = $trial_seq;
                $pos = $p;
                $match = 1;
                last;
            }
        }
        last unless $match;
    }

    @locs;
}


sub extend_array_forward
{
    my ( $seqR, $rept, $pos, $maxdiff, $maxreach ) = @_;
    my $replen = length( $rept );
    my $maxp = length( $$seqR ) - $replen;
    my @locs = ();
    while ( 1 )
    {
        my $match = 0;
        for ( my $i = $replen; $i <= $maxreach; $i++ )
        {
            my $p = $pos + $i;    #  trial position in $seq
            last if ( $p > $maxp );

            my $trial_seq = substr($$seqR,$p,$replen);
            if ( n_diff( $rept, $trial_seq ) <= $maxdiff )
            {
                push @locs, $p;
                $rept = $trial_seq;
                $pos = $p;
                $match = 1;
                last;
            }
        }
        last unless $match;
    }

    @locs;
}


=head2 Find the consensus residue, and number of state changes at a given
offset from the start of the candidate repeat sequence

    ( $residue, $n_change, $ttl ) = consensus_at_offset( \$seq, \@locs, $offset )

This is not used in the current version of the find_crisprs().

=cut

sub consensus_at_offset
{
    my ( $seqR, $locL, $offset ) = @_;

    my $len  = length( $$seqR );
    my $prev = '';
    my $nchg = 0;
    my $ttl  = 0;
    my %cnt;
    foreach ( @$locL )
    {
        my $loc = $_ + $offset;
        next if ( $loc < 0 ) || ( $loc >= $len );
        my $nt = substr($$seqR,$loc,1);
        next if $nt !~ /[acgt]/;
        $cnt{ $nt }++;
        $ttl++;
        $nchg++ if $prev && $prev ne $nt;
        $prev = $nt;
    }

    my ( $nt ) = sort { $cnt{$b} <=> $cnt{$a} } keys %cnt;

    ( $nt, $nchg, $ttl );
}


=head2 Propose a pairing pattern within the repeat sequence

Produce a pairing proposal for CRISPR repeat (maximizes number of canonical
pairs with no bulges or asymmetrical loops)

     $pairing = crispr_palindrome( $seq )

Example:

     seq      'gtttcaatccctaatagggatttaagttaattgcaac'
     pairing  '  << <<<<<<<   >>>>>>> >>            '

This function is not used in the current version of the find_crisprs().

=cut

sub crispr_palindrome
{
    my $seq = shift;    my $cmp = gjoseqlib::complement_DNA_seq( $seq );
    my $len = length( $seq );

    my $max_id = 0;
    my $i_max;

    #  Explore alignments with sequence shifted to right ($i_max < 0)
    #
    #  seq     ------------------
    #  cmp  ------------------

    for ( my $i = 1; $i <= $len - 10; $i++ )
    {
        my $ident = ( gjoseqlib::interpret_nt_align( substr($seq,0,$len-$i),
                                                     substr($cmp,$i,$len-$i) ) )[1];
        if ( $ident > $max_id ) { ( $max_id, $i_max ) = ( $ident, -$i ) }
    }

    #  Explore alignments with complement shifted to right ($i_max >= 0)
    #
    #  seq  ------------------
    #  cmp     ------------------

    for ( my $i = 0; $i <= $len - 20; $i++ )
    {
        my $ident = ( gjoseqlib::interpret_nt_align( substr($seq,$i,$len-$i),
                                                     substr($cmp,0,$len-$i) ) )[1];
        if ( $ident > $max_id ) { ( $max_id, $i_max ) = ( $ident, $i ) }
    }

    my $pairs = ' ' x $len;
    if ( $i_max < 0 )
    {
        #  seq     ------------------
        #  cmp  ------------------
        $cmp = substr($cmp,-$i_max);        # Shift the complement of simplify comparison
        my $over_end  = $len + $i_max - 1;  # End coordinate (0-based) of alignment
        my $j_mid     = 0.5 * $over_end;    # Midpoint of alignment
        for ( my $j = 0; $j < $j_mid - 1; $j++ )
        {
            if ( substr($seq,$j,1) eq substr($cmp,$j,1) )
            {
                substr($pairs,$j,1)           = '<';
                substr($pairs,$over_end-$j,1) = '>';
            }
        }
    }
    else
    {
        #  seq  ------------------
        #  cmp     ------------------
        $cmp = scalar( ' ' x $i_max ) . $cmp;         # Shift the complement of simplify comparison
        my $over_end = $len - 1;                      # End coordinate (0-based) of alignment
        my $j_mid    = 0.5 * ( $over_end + $i_max );  # Midpoint of alignment
        for ( my $j = $i_max; $j < $j_mid - 1; $j++ )
        {
            if ( substr($seq,$j,1) eq substr($cmp,$j,1) )
            {
                substr($pairs,$j,1)                  = '<';
                substr($pairs,$over_end-$j+$i_max,1) = '>';
            }
        }
    }

    $pairs;
}


=head2 Utility functions

=head3 Number of identical nontrivial nucleotides in align

     $nid = identical_nt( $seq1, $seq2 )

=head3 Number of differing residues in two sequences

     $n_diff = n_diff( $seq1, $seq2 )

=head3 Maximum of 2 numbers

     $max = max( $n1, $n2 )

=cut

sub identical_nt
{
    my ( $s1, $s2 ) = @_;
    $s1 =~ tr/acgt/\377/c;                # Disallowed symbols to x'FF' byte
    $s2 =~ tr/acgt/\376/c;                # Disallowed symbols to x'FE' byte
    scalar ( ( $s1^$s2 ) =~ tr/\000// );  # Count the nulls (identical residues)
}


sub n_diff
{
    my ( $s1, $s2 ) = @_;
    $s1 =~ tr/acgt/\377/c;                # Disallowed symbols to x'FF' byte
    $s2 =~ tr/acgt/\376/c;                # Disallowed symbols to x'FE' byte
    length($s1) - ( ( $s1^$s2 ) =~ tr/\000// );  # Count the nulls (identical residues)
}



sub max { $_[0] > $_[1] ? $_[0] : $_[1] }


1;
