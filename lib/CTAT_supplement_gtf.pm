package CTAT_supplement_gtf;

## Parsing, gene assignment and normalization of supplemental transcript models (GTF), shared by
## util/supplement_ctat_genome_lib.pl (supplementing a built CTAT genome lib) and
## util/merge_supplemental_gtfs.pl (adding them to the annotation of a lib being built).

use strict;
use warnings;
use Carp;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(parse_attributes load_gene_spans add_prior_supplement_genes build_gene_lookups load_transcript_ids load_ref_exons overlapping_gene_keys parse_supplemental_gtf assign_and_normalize assign_gene gene_label normalized_attributes print_report);


####
sub parse_attributes {
    my ($info) = @_;

    my @atts;
    ## the final attribute's ';' is optional
    while ($info =~ /\s*([^\s;]+)\s+(\"[^\"]*\"|[^;]+?)\s*(?:;|$)/g) {
        my ($key, $val) = ($1, $2);
        $val =~ s/^\"|\"$//g;
        push (@atts, [$key, $val]);
    }

    return(@atts);
}


####
sub load_gene_spans {
    my ($gene_spans_file) = @_;

    my %genes;
    open(my $fh, $gene_spans_file) or die "Error, cannot open $gene_spans_file";
    while (<$fh>) {
        chomp;
        unless (/\w/) { next; }
        my ($gene_id, $chr, $lend, $rend, $orient, $gene_name, $gene_type) = split(/\t/);
        ## gtf_to_gene_spans.pl lists a gene lacking a gene_name under its gene_id: it has no name
        if (! defined $gene_name || $gene_name eq "" || $gene_name eq $gene_id) {
            $gene_name = ".";
        }
        my $key = "$gene_id\t$chr";
        $genes{$key} = { gene_id => $gene_id,
                         chr => $chr,
                         lend => $lend,
                         rend => $rend,
                         orient => $orient,
                         gene_name => $gene_name,
                         gene_type => $gene_type,
                         novel => 0,
        };
    }
    close $fh;

    return(\%genes);
}


####
sub add_prior_supplement_genes {
    my ($genes_href, $prior_gtfs_aref) = @_;

    foreach my $gtf (@$prior_gtfs_aref) {
        open(my $fh, $gtf) or die "Error, cannot open $gtf";
        while (<$fh>) {
            chomp;
            my @x = split(/\t/);
            unless ($x[2] eq "exon") { next; }
            my %att = map { $_->[0] => $_->[1] } &parse_attributes($x[8]);
            my $key = "$att{gene_id}\t$x[0]";
            if (my $gene = $genes_href->{$key}) {
                if ($gene->{novel}) {
                    $gene->{lend} = $x[3] if $x[3] < $gene->{lend};
                    $gene->{rend} = $x[4] if $x[4] > $gene->{rend};
                }
                next;
            }
            $genes_href->{$key} = { gene_id => $att{gene_id},
                                    chr => $x[0],
                                    lend => $x[3],
                                    rend => $x[4],
                                    orient => $x[6],
                                    gene_name => ($att{gene_name} && $att{gene_name} ne $att{gene_id}) ? $att{gene_name} : ".",
                                    gene_type => $att{gene_type} || ".",
                                    novel => 1,
            };
        }
        close $fh;
    }

    return;
}


####
sub build_gene_lookups {
    my ($genes_href) = @_;

    my %by_id;
    my %by_id_nover;
    my %by_name;

    foreach my $key (keys %$genes_href) {
        my $gene = $genes_href->{$key};
        push (@{$by_id{$gene->{gene_id}}}, $key);

        ## version-stripped ids only for reference genes: ids of genes added by supplements (eg.
        ## StringTie's STRG.100) are not versioned, and stripping them would conflate unrelated genes
        unless ($gene->{novel}) {
            my $id_nover = $gene->{gene_id};
            if ($id_nover =~ s/\.\d+$//) {
                push (@{$by_id_nover{$id_nover}}, $key);
            }
        }

        if ($gene->{gene_name} ne ".") {
            push (@{$by_name{$gene->{gene_name}}}, $key);
        }
    }

    return( { by_id => \%by_id,
              by_id_nover => \%by_id_nover,
              by_name => \%by_name } );
}


####
sub load_transcript_ids {
    my ($gtf, $tids_href, $source) = @_;

    open(my $fh, $gtf) or die "Error, cannot open $gtf";
    while (<$fh>) {
        if (/transcript_id \"([^\"]+)\"/) {
            $tids_href->{$1} = $source;
        }
    }
    close $fh;

    return;
}


####
sub load_ref_exons {
    my ($exon_record_files_aref, $genes_href) = @_;

    ## map NAME^ID tokens back to gene keys
    my %token_to_key;
    foreach my $key (keys %$genes_href) {
        my $gene = $genes_href->{$key};
        if ($gene->{gene_name} ne ".") {
            $token_to_key{"$gene->{gene_name}^$gene->{gene_id}\t$gene->{chr}"} = $key;
        }
        else {
            $token_to_key{"$gene->{gene_id}\t$gene->{chr}"} = $key;
            ## as written for novel genes by earlier versions of this utility
            $token_to_key{"$gene->{gene_id}^$gene->{gene_id}\t$gene->{chr}"} = $key;
        }
    }

    my %exons; # chr -> orient -> [ [lend, rend, gene_key], ... ] sorted by lend
    foreach my $exon_records_file (@$exon_record_files_aref) {
        open(my $fh, $exon_records_file) or die "Error, cannot open $exon_records_file";
        while (<$fh>) {
            chomp;
            my @x = split(/\t/);
            unless ($x[8] =~ /gene_id \"([^\"]+)\"/) { next; }
            my $gene_key = $token_to_key{"$1\t$x[0]"} or next;
            push (@{$exons{$x[0]}->{$x[6]}}, [$x[3], $x[4], $gene_key]);
        }
        close $fh;
    }

    my %sorted_exons; # chr -> orient -> { exons => [ sorted by lend ], max_len => longest exon }
    foreach my $chr (keys %exons) {
        foreach my $orient (keys %{$exons{$chr}}) {
            my @sorted = sort { $a->[0] <=> $b->[0] } @{$exons{$chr}->{$orient}};
            my $max_len = 0;
            foreach my $exon (@sorted) {
                my $len = $exon->[1] - $exon->[0] + 1;
                $max_len = $len if $len > $max_len;
            }
            $sorted_exons{$chr}->{$orient} = { exons => \@sorted, max_len => $max_len };
        }
    }

    return(\%sorted_exons);
}


####
sub overlapping_gene_keys {
    my ($exon_set, $lend, $rend) = @_;

    my $exons_aref = $exon_set->{exons};

    ## binary search for the first exon that could reach $lend
    my $min_start = $lend - $exon_set->{max_len};
    my ($lo, $hi) = (0, scalar(@$exons_aref));
    while ($lo < $hi) {
        my $mid = int(($lo + $hi) / 2);
        if ($exons_aref->[$mid]->[0] < $min_start) {
            $lo = $mid + 1;
        }
        else {
            $hi = $mid;
        }
    }

    my @gene_keys;
    for (my $i = $lo; $i <= $#$exons_aref; $i++) {
        my $ref_exon = $exons_aref->[$i];
        if ($ref_exon->[0] > $rend) { last; }
        if ($ref_exon->[1] >= $lend) {
            push (@gene_keys, $ref_exon->[2]);
        }
    }

    return(@gene_keys);
}


####
sub parse_supplemental_gtf {
    my ($gtf) = @_;

    my %transcripts;
    my %report = ( skipped_unstranded => [],
                   dropped_gene_records => 0,
                   errors => [] );
    my $errors = $report{errors};
    my %gene_id_conflict;

    open(my $fh, $gtf) or die "Error, cannot open $gtf";
    read($fh, my $magic, 2);
    if (defined $magic && $magic eq "\x1f\x8b") {
        die "Error, $gtf is gzip-compressed; please provide it uncompressed.\n";
    }
    seek($fh, 0, 0);
    while (<$fh>) {
        if (/^\#/) { next; }
        unless (/\w/) { next; }
        chomp;
        s/\r$//;
        my @x = split(/\t/);
        unless (scalar(@x) == 9) {
            push (@$errors, "line $. does not have 9 tab-delimited fields: $_");
            next;
        }
        if ($x[2] eq "gene") {
            ## gene records are taken from the reference (or derived for novel genes).
            $report{dropped_gene_records}++;
            next;
        }
        unless ($x[3] =~ /^\d+$/ && $x[4] =~ /^\d+$/ && $x[3] <= $x[4]) {
            push (@$errors, "line $. has invalid coordinates: $_");
            next;
        }
        my @atts = &parse_attributes($x[8]);
        my %att = map { $_->[0] => $_->[1] } @atts;
        my $tid = $att{transcript_id};
        unless (defined $tid && $tid ne "") {
            push (@$errors, "line $. lacks a transcript_id: $_");
            next;
        }
        my $transcript = $transcripts{$tid};
        unless ($transcript) {
            $transcript = $transcripts{$tid} = { transcript_id => $tid,
                                                 chr => $x[0],
                                                 orient => $x[6],
                                                 atts => \%att,
                                                 lines => [],
                                                 exons => [],
                                                 order => scalar(keys %transcripts),
            };
        }
        if ($x[0] ne $transcript->{chr} || $x[6] ne $transcript->{orient}) {
            push (@$errors, "transcript $tid has records on more than one chromosome/strand (line $.)");
            $transcript->{bad} = 1;
            next;
        }
        ## the same transcript_id under different genes (eg. merged GTFs) must not be combined
        my $gene_id = $att{gene_id} // "";
        my $first_gene_id = $transcript->{atts}->{gene_id} // "";
        if ($gene_id ne "" && $first_gene_id ne "" && $gene_id ne $first_gene_id && ! $gene_id_conflict{$tid}) {
            push (@$errors, "transcript $tid has records with different gene_ids ('$first_gene_id', '$gene_id'; line $.)");
            $gene_id_conflict{$tid} = 1;
            $transcript->{bad} = 1;
        }
        foreach my $key (keys %att) {
            unless (exists $transcript->{atts}->{$key}) {
                $transcript->{atts}->{$key} = $att{$key};
            }
        }
        push (@{$transcript->{lines}}, [ @x[0..7], \@atts ]);
        if ($x[2] eq "exon") {
            push (@{$transcript->{exons}}, [$x[3], $x[4]]);
        }
    }
    close $fh;

    foreach my $tid (sort { $transcripts{$a}->{order} <=> $transcripts{$b}->{order} } keys %transcripts) {
        my $transcript = $transcripts{$tid};
        if ($transcript->{bad}) {
            delete $transcripts{$tid};
            next;
        }
        unless (@{$transcript->{exons}}) {
            push (@$errors, "transcript $tid has no exon records");
            delete $transcripts{$tid};
            next;
        }
        if ($transcript->{orient} ne '+' && $transcript->{orient} ne '-') {
            push (@{$report{skipped_unstranded}}, $tid);
            delete $transcripts{$tid};
        }
    }

    return(\%transcripts, \%report);
}


####
sub assign_and_normalize {
    my ($transcripts_href, $genes_href, $lookup, $existing_tids_href, $ref_exons_href, $opts) = @_;

    ## $opts: { source_tag, allow_novel_genes, novel_gene_type, allow_nonoverlapping }
    foreach my $opt ("source_tag", "novel_gene_type") {
        unless (defined $opts->{$opt}) {
            confess "Error, assign_and_normalize requires option: $opt";
        }
    }

    my %report = ( errors => [],
                   warnings => [],
                   num_transcripts_kept => 0,
                   existing_genes => {},
                   earlier_supplement_genes => {},
                   novel_genes => {},
                   assigned_by => {},
        );

    my %novel_gene_info; # gene_id -> { chr, orient }
    my %gene_to_transcripts;

    foreach my $tid (sort { $transcripts_href->{$a}->{order} <=> $transcripts_href->{$b}->{order} } keys %$transcripts_href) {
        my $transcript = $transcripts_href->{$tid};
        my $chr = $transcript->{chr};
        my $orient = $transcript->{orient};
        my $atts = $transcript->{atts};

        if ($tid =~ /[\|\^\s]/) {
            push (@{$report{errors}}, "transcript_id '$tid' contains a disallowed character ('|', '^' or whitespace)");
            next;
        }
        if (my $source = $existing_tids_href->{$tid}) {
            push (@{$report{errors}}, "transcript_id $tid already exists in the $source annotation");
            next;
        }

        my ($gene_key, $method, $problem) = &assign_gene($transcript, $genes_href, $lookup, $ref_exons_href);

        if ($problem) {
            push (@{$report{errors}}, "transcript $tid: $problem");
            next;
        }

        if ($gene_key) {
            my $gene = $genes_href->{$gene_key};
            if ($gene->{chr} ne $chr || $gene->{orient} ne $orient) {
                push (@{$report{errors}}, "transcript $tid ($chr:$orient) assigned by $method to gene "
                      . &gene_label($gene) . " which is on $gene->{chr}:$gene->{orient}");
                next;
            }
            ## catch likely mismatches (eg. a different genome build, or a reused id or name)
            my @coords = sort { $a <=> $b } map { @$_ } @{$transcript->{exons}};
            if (defined $gene->{lend} && ($coords[-1] < $gene->{lend} || $coords[0] > $gene->{rend})) {
                my $msg = "transcript $tid ($chr:$coords[0]-$coords[-1]) does not overlap the gene it was assigned to by $method, "
                    . &gene_label($gene) . " ($gene->{chr}:$gene->{lend}-$gene->{rend})";
                if ($opts->{allow_nonoverlapping}) {
                    push (@{$report{warnings}}, "$msg; the gene span will be extended to include it");
                }
                else {
                    push (@{$report{errors}}, "$msg (see --allow_nonoverlapping)");
                    next;
                }
            }

            $report{assigned_by}->{$method}++;
            if ($gene->{novel}) {
                $report{earlier_supplement_genes}->{$gene_key} = 1;
            }
            else {
                $report{existing_genes}->{$gene_key} = 1;
            }
        }
        else {
            unless ($opts->{allow_novel_genes}) {
                push (@{$report{errors}}, "transcript $tid (gene_id: " . ($atts->{gene_id} || "none") . ") "
                      . "could not be assigned to a reference gene (see --assign_by_overlap, --allow_novel_genes)");
                next;
            }
            my $gene_id = $atts->{gene_id};
            unless ($gene_id) {
                push (@{$report{errors}}, "transcript $tid has no gene_id to use for a novel gene");
                next;
            }
            if ($gene_id =~ /[\|\^\s]/ || $gene_id =~ /--/) {
                push (@{$report{errors}}, "novel gene_id '$gene_id' (transcript $tid) contains a disallowed character ('|', '^', '--' or whitespace)");
                next;
            }
            if ($lookup->{by_name}->{$gene_id}) {
                push (@{$report{errors}}, "novel gene_id '$gene_id' (transcript $tid) collides with an existing gene_name");
                next;
            }
            if (my $info = $novel_gene_info{$gene_id}) {
                if ($info->{chr} ne $chr || $info->{orient} ne $orient) {
                    push (@{$report{errors}}, "novel gene $gene_id has transcripts on more than one chromosome/strand");
                    next;
                }
            }
            $novel_gene_info{$gene_id} = { chr => $chr, orient => $orient };
            $gene_key = "$gene_id\t$chr";
            $genes_href->{$gene_key} = { gene_id => $gene_id,
                                         chr => $chr,
                                         orient => $orient,
                                         gene_name => ".", # listed under its gene_id, as for any gene without a name
                                         gene_type => $opts->{novel_gene_type},
                                         novel => 1,
            };
            $report{novel_genes}->{$gene_id} = 1;
            $report{assigned_by}->{"novel gene"}++;
        }

        push (@{$gene_to_transcripts{$gene_key}}, $tid);
    }

    ## write normalized records, grouped by gene and transcript
    my @lines;
    foreach my $gene_key (sort keys %gene_to_transcripts) {
        my $gene = $genes_href->{$gene_key};
        foreach my $tid (@{$gene_to_transcripts{$gene_key}}) {
            my $transcript = $transcripts_href->{$tid};
            $report{num_transcripts_kept}++;

            my @recs = sort { $a->[3] <=> $b->[3]
                                  || ($a->[2] eq "transcript" ? 0 : 1) <=> ($b->[2] eq "transcript" ? 0 : 1)
                                  || $a->[2] cmp $b->[2] } @{$transcript->{lines}};

            foreach my $rec (@recs) {
                my @x = @$rec[0..7];
                $x[1] = $opts->{source_tag};
                $x[8] = &normalized_attributes($gene, $tid, $rec->[8]);
                push (@lines, join("\t", @x));
            }
        }
    }

    return(\@lines, \%report);
}


####
sub assign_gene {
    my ($transcript, $genes_href, $lookup, $ref_exons_href) = @_;

    my $chr = $transcript->{chr};
    my $orient = $transcript->{orient};
    my $atts = $transcript->{atts};

    my $pick = sub {
        my ($keys_aref) = @_;
        ## prefer candidates on the transcript's chromosome
        my @on_chr = grep { $genes_href->{$_}->{chr} eq $chr } @$keys_aref;
        if (scalar(@on_chr) == 1) {
            return($on_chr[0]);
        }
        elsif (scalar(@$keys_aref) == 1) {
            return($keys_aref->[0]);
        }
        return(undef);
    };

    ## by gene_id, then ref_gene_id
    foreach my $att ("gene_id", "ref_gene_id") {
        my $val = $atts->{$att} or next;
        if (my $keys_aref = $lookup->{by_id}->{$val}) {
            my $key = &$pick($keys_aref);
            return($key, $att) if $key;
            my @chrs = sort map { $genes_href->{$_}->{chr} } @$keys_aref;
            return(undef, undef, "$att '$val' matches genes on more than one other chromosome (@chrs), none on $chr");
        }
    }

    ## by version-stripped gene_id / ref_gene_id (reference genes only), requiring chr and strand agreement
    foreach my $att ("gene_id", "ref_gene_id") {
        my $val = $atts->{$att} or next;
        $val =~ s/\.\d+$//;
        if (my $keys_aref = $lookup->{by_id_nover}->{$val}) {
            my @candidates = grep { $genes_href->{$_}->{chr} eq $chr
                                        && $genes_href->{$_}->{orient} eq $orient } @$keys_aref;
            if (scalar(@candidates) == 1) {
                return($candidates[0], "$att (version-stripped)");
            }
            elsif (scalar(@candidates) > 1) {
                return(undef, undef, "$att '$atts->{$att}' (version-stripped) matches more than one reference gene on $chr:$orient");
            }
            ## otherwise fall through to gene_name matching
        }
    }

    ## by gene name, requiring chr and strand agreement
    foreach my $att ("gene_name", "ref_gene_name") {
        my $val = $atts->{$att} or next;
        if (my $keys_aref = $lookup->{by_name}->{$val}) {
            my @candidates = grep { $genes_href->{$_}->{chr} eq $chr
                                        && $genes_href->{$_}->{orient} eq $orient } @$keys_aref;
            if (scalar(@candidates) == 1) {
                return($candidates[0], $att);
            }
            elsif (scalar(@candidates) > 1) {
                return(undef, undef, "$att '$val' matches more than one reference gene on $chr:$orient");
            }
            else {
                return(undef, undef, "$att '$val' matches a reference gene, but not on $chr:$orient");
            }
        }
    }

    ## by exon overlap
    if ($ref_exons_href) {
        my %overlapping_genes;
        if (my $exon_set = $ref_exons_href->{$chr}->{$orient}) {
            foreach my $exon (@{$transcript->{exons}}) {
                foreach my $gene_key (&overlapping_gene_keys($exon_set, @$exon)) {
                    $overlapping_genes{$gene_key} = 1;
                }
            }
        }
        my @candidates = keys %overlapping_genes;
        if (scalar(@candidates) == 1) {
            return($candidates[0], "exon overlap");
        }
        elsif (scalar(@candidates) > 1) {
            my @names = sort map { &gene_label($genes_href->{$_}) } @candidates;
            return(undef, undef, "exons overlap more than one reference gene (@names)");
        }
    }

    return(undef);
}


####
sub gene_label {
    my ($gene) = @_;
    return( ($gene->{gene_name} ne ".") ? "$gene->{gene_name} ($gene->{gene_id})" : $gene->{gene_id} );
}


####
sub normalized_attributes {
    my ($gene, $tid, $atts_aref) = @_;

    my @parts = ("gene_id \"$gene->{gene_id}\"",
                 "transcript_id \"$tid\"");
    if ($gene->{gene_name} ne ".") {
        push (@parts, "gene_name \"$gene->{gene_name}\"");
    }
    if ($gene->{gene_type} ne ".") {
        push (@parts, "gene_type \"$gene->{gene_type}\"");
    }

    my %skip = map { $_ => 1 } qw(gene_id transcript_id gene_name gene_type gene_biotype
                                  ref_gene_id ref_gene_name);
    my $have_transcript_type = 0;
    foreach my $att (@$atts_aref) {
        my ($key, $val) = @$att;
        if ($skip{$key}) { next; }
        if ($key eq "transcript_type") { $have_transcript_type = 1; }
        push (@parts, "$key \"$val\"");
    }
    unless ($have_transcript_type) {
        push (@parts, "transcript_type \"supplemental\"");
    }

    return(join("; ", @parts) . ";");
}


####
sub print_report {
    my ($report) = @_;

    print STDERR "\n#### supplemental annotation summary ####\n";
    print STDERR "  transcripts kept: $report->{num_transcripts_kept}\n";
    print STDERR "  existing genes receiving isoforms: " . scalar(keys %{$report->{existing_genes}}) . "\n";
    if (my $num = scalar(keys %{$report->{earlier_supplement_genes}})) {
        print STDERR "  genes from earlier supplements receiving isoforms: $num\n";
    }
    print STDERR "  novel genes: " . scalar(keys %{$report->{novel_genes}}) . "\n";
    foreach my $method (sort keys %{$report->{assigned_by}}) {
        print STDERR "    assigned by $method: $report->{assigned_by}->{$method}\n";
    }
    if (my $num = scalar(@{$report->{skipped_unstranded}})) {
        print STDERR "  skipped (no strand): $num\n";
    }
    if ($report->{dropped_gene_records}) {
        print STDERR "  'gene' records ignored: $report->{dropped_gene_records}\n";
    }
    if (my @warnings = @{$report->{warnings}}) {
        print STDERR "  WARNINGS: " . scalar(@warnings) . "\n";
        foreach my $warning (@warnings[0 .. ($#warnings < 24 ? $#warnings : 24)]) {
            print STDERR "    $warning\n";
        }
        print STDERR "    ...\n" if (@warnings > 25);
    }
    print STDERR "  problems: " . scalar(@{$report->{errors}}) . "\n\n";

    return;
}

1;
