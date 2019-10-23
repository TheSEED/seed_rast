#
# Copyright (c) 2003-2015 University of Chicago and Fellowship
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
#

#
# This is a SAS component.
#


package SubsystemProjector;

    use strict;
    use warnings;
    use Stats;
    use ServicesUtils;
    use RoleParse;
use Data::Dumper;

=head1 Find Subsystems for a Genome

This object is used to perform a basic role-based projection of subsystems onto genomes. The
core method takes as input a pair of hashes and a table of feature IDs to function IDs. The
basic algorithm is then as follows:

=over 4

=item 1

Use the functions to compute the subsystems for each feature.

=item 2

For each subsystem found, match the set of roles therein to known variants. If all of the roles in a
variant are present, then it is a I<candidate variant> for the subsystem.

=item 3

Choose the candidate variant with the most roles and output the features in the subsystem based on
that variant.

=back

The object contains the following fields.

=over 4

=item roleMap

Reference to a hash mapping each role checksum to a list of subsystem names.

=item variantMap

Reference to a hash mapping each subsystem name to a list of variant tuples. Each variant tuple consists of
a variant code followed by a list of role checksums.

=item stats

A statistics object containing statistics for this process.

=back

=head2 Special Methods

    my $projector = SubsystemProjector->new($roleFile, $variantFile);

Create a projector object. The projector object contains a map of roles to subsystems and a map of
subsystems to variants. Each variant is a list of roles. Each role is stored in the form of the checksum
computed by L<RoleParse>.

=over 4

=item roleFile

Name of a tab-delimited file containing [role checksum, subsystem name] pairs.

=item variantFile

Name of a tab-delimited file containing in each record (0) a subsystem name, (1) a variant code, and
(2) a space-delimited list of role checksums.

=back

=cut

sub new {
    my ($class, $roleFile, $variantFile) = @_;
    # Get the statistics object.
    my $stats = Stats->new();
    # Start by processing the role file.
    my %roleHash;
    open(my $ih, '<', $roleFile) || die "Could not open role file $roleFile: $!";
    while (! eof $ih) {
        my $line = <$ih>;
        chomp $line;
        $stats->Add(roleFileIn => 1);
        my ($roleID, $ssName) = split /\t/, $line;
        push @{$roleHash{$roleID}}, $ssName;
    }
    close $ih; undef $ih;
    # Now process the variant map.
    my %variantMap;
    open($ih, '<', $variantFile) || die "Could not open variant file $variantFile: $!";
    while (! eof $ih) {
        my $line = <$ih>;
        chomp $line;
        $stats->Add(variantFileIn => 1);
        my ($ssName, $vCode, $roles) = split /\t/, $line;
        my @roleList = split /\s+/, $roles;
        push @{$variantMap{$ssName}}, [$vCode, @roleList];
    }
    # Create the projector object.
    my $retVal = {
        variantMap => \%variantMap,
        roleMap => \%roleHash,
        stats => $stats
    };
    # Bless and return it.
    bless $retVal, $class;
    return $retVal;
}

=head2 Query Methods

=head3 stats

    my $stats = $projector->stats;

Get the statistics object.

=cut

sub stats {
    return $_[0]->{stats};
}

=head2 Public Manipulation Methods

=head3 Project

    my $subsystemHash = $projector->Project(\%featureAssignments);

Compute the subsystems that occur for a specified set of annotated
features.

=over 4

=item featureAssignments

Reference to a hash mapping each feature ID to its functional assignment.

=item RETURN

Returns a reference to a hash mapping each subsystem name to a 2-tuple containing (0) a variant code, and (1) a
reference to a list of [role, fid] tuples.

=back

=cut

sub Project {
    # Get the parameters.
    my ($self, $featureAssignments) = @_;
    # Get the statistics object.
    my $stats = $self->{stats};
    # Declare the return variable.
    my %retVal;
    # Create a map of role checksums to role names.
    my %roleNames;
    # Create a hash of all the roles in the genome. Each role checksum will be mapped to a list of feature IDs.
    my %roleFids;
    for my $fid (keys %$featureAssignments) {
        $stats->Add(featuresIn => 1);
        my $function = $featureAssignments->{$fid};
        my @roles = SeedUtils::roles_of_function($function);
        for my $role (@roles) {
            $stats->Add(rolesIn => 1);
            my $roleID = RoleParse::Checksum($role);
            $roleNames{$roleID} = $role;
            push @{$roleFids{$roleID}}, $fid;
        }
    }

    # We need all the subsystems containing these roles. We map each subsystem to a hash of the
    # roles from that subsystem found in this genome.
    my $roleMap = $self->{roleMap};
    my %subs;
    for my $roleID (keys %roleNames) {
        my $subList = $roleMap->{$roleID};
        for my $sub (@$subList) {
            $subs{$sub}{$roleID} = 1;
        }
    }

    # Now process each of the subsystems. For each subsystem we have a sub-hash of all its roles
    # currently in the genome. We want the best match, that is, the variant whose roles are fully
    # represented and has the most roles in it.
    my $variantMap = $self->{variantMap};
    for my $sub (sort keys %subs) {
        # These variables will contain the best match so far.
        my ($bestVariant, $bestRoles);
        # This is the role count for the best match.
        my $bestCount = 0;
        # Get the hash of represented roles in this subsystem.
        my $subRolesH = $subs{$sub};
        my $represented = scalar keys %$subRolesH;
        # Get all the maps for this subsystem.
        my $maps = $variantMap->{$sub} // [];

        my @miss_info;

        for my $map (@$maps) {
            my ($variant, @roleIDs) = @$map;
            my $count = scalar @roleIDs;
            # Do we have enough represented roles to fill this variant?
            if ($count <= $represented) {
                # Yes. Count the roles found.
                my @found_ids = grep { $subRolesH->{$_} } @roleIDs;
                my $found = @found_ids;
                # print STDERR "Found=$found count=$count $variant @roleIDs\n";
                if ($found == $count) {
                    # Here all the roles in the map are represented in the genome.
                    if ($count > $bestCount) {
                        # Here this match is the best one found so far.
                        ($bestVariant, $bestRoles, $bestCount) = ($variant, \@roleIDs, $count);
                    }
                } else {
                    push(@miss_info, [$found, $count, [@roleIDs], [@found_ids]]);
                }
            } elsif (0) {
                if ($variant eq '1.2022')
                {
                    print "Miss for $variant with count=$count rep=$represented\n";
                    for my $rid (@roleIDs)
                    {
                        my $r = $roleNames{$rid};
                        print join("\t", "@{$roleFids{$rid}}", $rid, $r), "\n";
                    }
                }
            }
        }
        # Did we find a match?
        if ($bestCount) {
            # Yes. Create the variant description.
            my @variantRoles;
            # print STDERR Dumper(BEST => $bestRoles);
            for my $roleID (@$bestRoles) {
                my $rolePegs = $roleFids{$roleID};
                my $role = $roleNames{$roleID};
                for my $peg (@$rolePegs) {
                    push @variantRoles, [$role, $peg];
                }
            }
            $retVal{$sub} = [$bestVariant, \@variantRoles];
        }

        if (0) {
            no warnings;
            print STDERR "$sub mismatch: (" . scalar @miss_info . ") entries:\n";
            for my $miss (@miss_info)
            {
                my($found, $count, $role_ids, $found_ids) = @$miss;
                print STDERR "\tfound=$found count=$count\n";
                print STDERR "\troles=" . join(" ", map { "$_:'" . $roleNames{$_} . "'" } @$role_ids) . "\n";
                print STDERR "\tfound=" . join(" ", map { "$_:'" .  $roleNames{$_} . "'" } @$found_ids) . "\n";
            }
        }
    }
    # Return the result.
    return \%retVal;
}


=head3 ProjectForGto

    my $subsystemHash = $projector->ProjectForGto($gto, %options);

Compute the subsystems that occur in a genome defined by a L<GenomeTypeObject>. This method essentially
computes the feature assignment hash and then calls L</Project>.

=over 4

=item gto

L<GenomeTypeObject> for the genome on which the subsystems should be projected.

=item options

A hash of options, including zero or more of the following.

=over 8

=item store

If TRUE, the subsystems will be stored directory into the GenomeTypeObject. The default is FALSE.

=back

=item RETURN

Returns a reference to a hash mapping each subsystem name to a 2-tuple containing (0) a variant code, and (1) a
reference to a list of [role, fid] tuples.

=back

=cut

sub ProjectForGto {
    my ($self, $gto, %options) = @_;
    # Get the feature list.
    my $featureList = ServicesUtils::json_field($gto, 'features');
    # Loop through the features, creating the assignment hash.
    my %assigns;
    for my $featureData (@$featureList) {
        my $fid = $featureData->{id};
        $assigns{$fid} = $featureData->{function};
    }
    # Project the subsystems.
    my $retVal = $self->Project(\%assigns);
    # Store the results if needed.
    if ($options{store}) {
        my %subs;
        for my $sub (keys %$retVal) {
            my $projectionData = $retVal->{$sub};
            my ($variant, $subRow) = @$projectionData;
            my %cells;
            for my $subCell (@$subRow) {
                my ($role, $fid) = @$subCell;
                push @{$cells{$role}}, $fid;
            }
            $subs{$sub} = [$variant, \%cells];
        }
        $gto->{subsystems} = \%subs;
    }
    return $retVal;
}


1;
