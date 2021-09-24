#!/usr/bin/perl -w

use strict;


my $data;
my $db;
my $aligner;
my %files;
my %functions;
my $headings;
foreach my $file (@ARGV) {
	open(IN, $file) || die "$! $file";
	my $query = <IN>;
	my $thisdb = <IN>;
	my $thisaligner = <IN>;
	if (!defined $db) {$db = $thisdb}
	if ($thisdb ne $db) {
		print STDERR "WARNING: One file used $db and another used $thisdb. Are they the same?\n";
	}
	if (!defined $aligner) {$aligner = $thisaligner}
	if ($aligner ne $thisaligner) {
		print STDERR "WARNING: One file used $aligner and another used $thisaligner. Are they the same?\n";
	}
	my $header = <IN>;
	$header = <IN>;
	chomp($header);
 	my @header = split /\t/, $header;
	my $perh = pop(@header);
	my $coh = pop(@header);
	my $filename = "$coh\t$perh";
	$headings = join("\t", @header);
	$files{$filename}=1;
	while (<IN>) {
		chomp;
		my @a=split /\t/;
		my $per = pop(@a);
		my $co = pop(@a);
		my $key = join("\t", @a);
		$data->{$filename}->{$key} = "$co\t$per";
		$functions{$key}=1;
	}
	close IN;
}

my @filenames = sort {$a cmp $b} keys %files;
print "$headings\t", join("\t", @filenames), "\n";
foreach my $k (sort {$a cmp $b} keys %functions) {
	print $k;
	foreach my $f (@filenames) {
		if ($data->{$f}->{$k}) {
			print "\t", $data->{$f}->{$k};
		} else {
			print "\t0";
		}
	}
	print "\n";
}

		
