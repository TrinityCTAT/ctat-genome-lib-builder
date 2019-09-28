#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use Data::Dumper;
use JSON::XS;
use FindBin;
use lib ("$FindBin::Bin/../lib");
use Fasta_reader;
use TiedHash;
use Process_cmd;
use JSON::XS;

my $usage = <<__EOUSAGE__;

##########################################################################
#
#  --pfam_domains <string>    : pfam output  (ok if .gz)
#
#  --genome_lib_dir <string>  : CTAT genome lib directory
#
##########################################################################


 
__EOUSAGE__

    ;


my $help_flag;
my $pfam_domains_file;
my $genome_lib_dir;

&GetOptions ( 'h' => \$help_flag,
              'pfam_domains=s' => \$pfam_domains_file,
              'genome_lib_dir=s' => \$genome_lib_dir,
    );

unless ($pfam_domains_file && $genome_lib_dir) {
    die $usage;
}



main: {

    my %pfam_hits = &parse_pfam($pfam_domains_file);



    print STDERR "-building pfam dbm\n";

    my $tied_hash = new TiedHash( { create => "$genome_lib_dir/pfam_domains.dbm" } );

    
    foreach my $isoform (keys %pfam_hits) {
        my @pfam_domains = @{$pfam_hits{$isoform}};

        my $json = encode_json(\@pfam_domains);

        #print "$isoform\t$json\n";
    
        $tied_hash->store_key_value($isoform, $json);

    }

    $tied_hash->store_key_value("geneABC--geneXYZ", "__placeholder_testval__");
    
    print STDERR "done.\n";
    
    exit(0);
}


=pfam_entry

0       BTG
1       PF07742.9
2       116
3       BTG2|ENST00000290551.4
4       -
5       158
6       3.8e-42
7       142.7
8       0.0
9       1
10      1
11      2.8e-46
12      4.5e-42
13      142.4
14      0.0
15      1
16      115
17      9
18      121
19      9
20      122
21      0.99
22      BTG
23      family

=cut


####
sub parse_pfam {
    my ($pfam_file) = @_;

    print STDERR "-parsing $pfam_file\n";
    
    my %model_to_domains;

    my $fh;
    if ($pfam_file =~ /\.gz$/) {
        open($fh, "gunzip -c $pfam_file | ") or die "Error, cannot open $pfam_file via gunzip -c";
    }
    else {
        open ($fh, $pfam_file) or die "Error, cannot open file $pfam_file";
    }

    while (<$fh>) {
        chomp;
        unless (/\w/) { next; }
        if (/^\#/) { next; }
        my @x = split(/\s+/);
        
        if (scalar @x < 22) {
            print STDERR "WARNING: Skipping line: $_ as likely corrupt.\n";
            next;
        }
        
        my $QueryProtID = $x[3];
        my $pfam_id = $x[1];
        my $HMMERDomain = $x[0];
        my $HMMERTDomainDescription = join(" ", @x[22..$#x]);
        my $QueryStartAlign = $x[17];
        my $QueryEndAlign = $x[18];
        my $PFAMStartAlign = $x[15];
        my $PFAMEndAlign = $x[16];
        my $FullSeqEvalue = $x[6];
        my $ThisDomainEvalue = $x[11];
        my $FullSeqScore = $x[7];
        my $FullDomainScore = $x[13];
        
        my $pfam_domain_struct = { cds_id => $QueryProtID,
                                   hmmer_domain => $HMMERDomain,
                                   query_start => $QueryStartAlign,
                                   query_end => $QueryEndAlign,
                                   domain_evalue => $ThisDomainEvalue,
        };
        
        push (@{$model_to_domains{$QueryProtID}}, $pfam_domain_struct);
    }
    close $fh;

    return(%model_to_domains);
}
        
        
