#!/usr/bin/perl -w
use strict;
use DBI;
use lib "lib/";
use lib "lib/pangu";
use pangu;

## this script is to delete alignment files

my $commoncfg = readConfig("main.conf");

print "Please type a sequence ID:";
chop(my $sequenceId=<STDIN>);

my $queryDir;
for (my $position = 0; $position < length($sequenceId); $position += 2)
{
	$queryDir .= "/q". substr($sequenceId,$position,2);
}
if(-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir")
{
	my @dirs;
	push @dirs, "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir";

	do{
		my $dir = shift @dirs;
		opendir(DIR, $dir) or die "can't opendir $dir: $!";
		my @files = readdir(DIR);
		closedir DIR;
		foreach my $file (sort @files)
		{
			next if ($file =~ /^\./);
			if(-f "$dir/$file")
			{
				if ($file =~ /(\d+)-(\d+).tbl$/)
				{
					print "Deleting alignments in $dir/$file ... ";
					#unlink ("$file"); #delete old alignments
					print "done!\n";
					my $queryId = $1;
					my $subjectId = $2;
					my $queryDirSwitched;
					for (my $position = 0; $position < length($subjectId); $position += 2)
					{
						$queryDirSwitched .= "/q". substr($subjectId,$position,2);
					}

					my $subjectDirSwitched;
					for (my $position = 0; $position < length($queryId); $position += 2)
					{
						$subjectDirSwitched .= "/s". substr($queryId,$position,2);
					}

					if(-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$subjectId-$queryId.tbl")
					{
						print "Deleting alignments in $commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$subjectId-$queryId.tbl ... ";
						#unlink ("$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$subjectId-$queryId.tbl"); #delete old alignments
						print "done!\n";
					}
				}
			}
			elsif(-d "$dir/$file")
			{
				unshift @dirs, "$dir/$file";
			}
		}
	}while (@dirs);
}
