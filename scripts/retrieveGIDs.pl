#!/usr/bin/perl -w

## read a list of sequence GI and retrieve the taxon ID for this sequence
## from the huge gi_taxid mapping 


use strict;
use Getopt::Long;
use File::Spec qw(catpath);
use File::Basename;
#use IO::Zlib; # TODO: read gzip files to avoid unpacking


my $hugefile="taxonomy/gi_taxid_nucl.dmp";
my $sep=",";
my $col=0;

usage() unless GetOptions("f=s"=> \$hugefile,
			  "d=s"=> \$sep,
			  "c=i"=> \$col,
			  "h"   => \&usage);

sub usage
{    
    printf STDERR "\n";
    printf STDERR "Reads  a list of sequence IDs, scans a large table file \n";
    printf STDERR "and prints lines where a specified field equals the ID\n\n";

    printf STDERR "\tUsage: @{[basename($0)]} [options] <ID file>\n";
    printf STDERR "\n";
    
    printf STDERR "\t-f=s\tfile to scan, current: $hugefile\n";
    printf STDERR "\t-d=s\tcolumn separator, current: $sep\n";
    printf STDERR "\t-c=s\tcolumn number where ID is searched, current: $col\n";
    printf STDERR "\n";
    exit;
}

my %IDs;
my $cnt=0;
while(<>) {
    chomp;
    print STDERR "$_ already registered\n" if ( defined($IDs{$_}) );
    $IDs{$_} = "find";
    #print STDERR $_." registered\n";
    $cnt++;
}
print STDERR (keys %IDs)." IDs registered ($cnt).\n";

open(HUGE, "<".$hugefile) or die "can't find $hugefile\n";

$cnt=0;
while(<HUGE>) {
    my @vals = split $sep;
    if ( defined($IDs{$vals[0]}) ) {
	print $_ ;
	$cnt++;
	$IDs{$vals[0]} = "found";
    }    
    last if $cnt == scalar keys %IDs;
}
print STDERR "done\n";
close(HUGE);

foreach my $id ( keys %IDs ) {
    print STDERR $id." not found!\n" if $IDs{$id} eq "find";
}
