#!/usr/bin/perl -w
use strict;
use DBI;
use lib "lib/";
use lib "lib/pangu";
use pangu;

## this script is to export records in table alignment into files

my $commoncfg = readConfig("main.conf");
my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});

## sequences will be saved to $commoncfg->{DATADIR}/alignments
until (-e "$commoncfg->{DATADIR}/alignments")
{
	print "Creating directory $commoncfg->{DATADIR}/alignments ... ";
	mkdir "$commoncfg->{DATADIR}/alignments";
	print "done!\n";
}
until (-e "$commoncfg->{DATADIR}/alignments/seqToSeq")
{
	print "Creating directory $commoncfg->{DATADIR}/alignments/seqToSeq ... ";
	mkdir "$commoncfg->{DATADIR}/alignments/seqToSeq";
	print "done!\n";
}
until (-e "$commoncfg->{DATADIR}/alignments/seqToSet")
{
	print "Creating directory $commoncfg->{DATADIR}/alignments/seqToSet ... ";
	mkdir "$commoncfg->{DATADIR}/alignments/seqToSet";
	print "done!\n";
}
until (-e "$commoncfg->{DATADIR}/alignments/setToSet")
{
	print "Creating directory $commoncfg->{DATADIR}/alignments/setToSet ... ";
	mkdir "$commoncfg->{DATADIR}/alignments/setToSet";
	print "done!\n";
}

print "Please check $commoncfg->{DATADIR}/alignments/alignmentToFile.log for exporting details.\n";
open (LOG,">$commoncfg->{DATADIR}/alignments/alignmentToFile.log") or die "can't open file: $commoncfg->{DATADIR}/alignments/alignmentToFile.log";

my $setId;
my $getSequences = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence'");
$getSequences->execute();
while(my @getSequences = $getSequences->fetchrow_array())
{
	$setId->{$getSequences[0]} = $getSequences[4];
}

my $seqToSeqNumber = 0;
my $seqToSet;
my $setToSet;
my $alnNumber = 0;
my $getAlignments = $dbh->prepare("SELECT * FROM alignment");
$getAlignments->execute();
while(my @getAlignments = $getAlignments->fetchrow_array())
{
	$alnNumber++;
	my $queryDir;
	for (my $position = 0; $position < length($getAlignments[2]); $position += 2)
	{
		$queryDir .= "/q". substr($getAlignments[2],$position,2);
		until (-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir")
		{
			print LOG "Creating directory $commoncfg->{DATADIR}/alignments/seqToSeq$queryDir ... ";
			mkdir "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir";
			print LOG "done!\n";
		}
	}

	my $subjectDir;
	for (my $position = 0; $position < length($getAlignments[3]); $position += 2)
	{
		$subjectDir .= "/s". substr($getAlignments[3],$position,2);
		until (-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir")
		{
			print LOG "Creating directory $commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir ... ";
			mkdir "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir";
			print LOG "done!\n";
		}
	}

	until (-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$getAlignments[2]-$getAlignments[3].tbl")
	{
		print LOG "Creating alignment file $commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$getAlignments[2]-$getAlignments[3].tbl ... ";
		open (ALN,">$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$getAlignments[2]-$getAlignments[3].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$getAlignments[2]-$getAlignments[3].tbl";
		print ALN "#$getAlignments[1]\n";
		print ALN "#query\tsubject\tperc_indentity\talign_length\tmismatches\tgaps\tq_start\tq_end\ts_start\ts_end\te_val\tbit_score\thidden\n";
		close(ALN);
		print LOG "done!\n";
		$seqToSeqNumber++;
		if (exists $seqToSet->{$getAlignments[2]}->{$setId->{$getAlignments[3]}})
		{
			$seqToSet->{$getAlignments[2]}->{$setId->{$getAlignments[3]}} .= ",alignments/seqToSeq$queryDir$subjectDir/$getAlignments[2]-$getAlignments[3].tbl";
		}
		else
		{
			$seqToSet->{$getAlignments[2]}->{$setId->{$getAlignments[3]}} = "alignments/seqToSeq$queryDir$subjectDir/$getAlignments[2]-$getAlignments[3].tbl";
		}
		if (exists $setToSet->{$setId->{$getAlignments[2]}}->{$setId->{$getAlignments[3]}})
		{
			$setToSet->{$setId->{$getAlignments[2]}}->{$setId->{$getAlignments[3]}} .= ",alignments/seqToSeq$queryDir$subjectDir/$getAlignments[2]-$getAlignments[3].tbl";
		}
		else
		{
			$setToSet->{$setId->{$getAlignments[2]}}->{$setId->{$getAlignments[3]}} = "alignments/seqToSeq$queryDir$subjectDir/$getAlignments[2]-$getAlignments[3].tbl";
		}
	}

	$getAlignments[13] =~ s/\W//g;
	print LOG "Exporting No. $alnNumber alignment '$getAlignments[2]-$getAlignments[3]' ... ";
	open (ALN,">>$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$getAlignments[2]-$getAlignments[3].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$getAlignments[2]-$getAlignments[3].tbl";
	print ALN "$getAlignments[2]\t$getAlignments[3]\t$getAlignments[4]\t$getAlignments[5]\t$getAlignments[6]\t$getAlignments[7]\t$getAlignments[8]\t$getAlignments[9]\t$getAlignments[10]\t$getAlignments[11]\t$getAlignments[12]\t$getAlignments[13]\t$getAlignments[14]\n";
	close(ALN);
	print LOG "done!\n";
}

my $seqToSetNumber = 0;
foreach my $sequenceId (keys %$seqToSet)
{
	my $queryDir;
	for (my $position = 0; $position < length($sequenceId); $position += 2)
	{
		$queryDir .= "/q". substr($sequenceId,$position,2);
		until (-e "$commoncfg->{DATADIR}/alignments/seqToSet$queryDir")
		{
			print LOG "Creating directory $commoncfg->{DATADIR}/alignments/seqToSet$queryDir ... ";
			mkdir "$commoncfg->{DATADIR}/alignments/seqToSet$queryDir";
			print LOG "done!\n";
		}
	}
	foreach my $subjectSetId (keys %{$seqToSet->{$sequenceId}})
	{

		my $subjectDir;
		for (my $position = 0; $position < length($subjectSetId); $position += 2)
		{
			$subjectDir .= "/s". substr($subjectSetId,$position,2);
			until (-e "$commoncfg->{DATADIR}/alignments/seqToSet$queryDir$subjectDir")
			{
				print LOG "Creating directory $commoncfg->{DATADIR}/alignments/seqToSet$queryDir$subjectDir ... ";
				mkdir "$commoncfg->{DATADIR}/alignments/seqToSet$queryDir$subjectDir";
				print LOG "done!\n";
			}
		}
		
		until (-e "$commoncfg->{DATADIR}/alignments/seqToSet$queryDir$subjectDir/$sequenceId-$subjectSetId.list")
		{
			print LOG "Creating alignment seqToSet list $commoncfg->{DATADIR}/alignments/seqToSet$queryDir$subjectDir/$sequenceId-$subjectSetId.list ... ";
			open (LIST,">$commoncfg->{DATADIR}/alignments/seqToSet$queryDir$subjectDir/$sequenceId-$subjectSetId.list") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSet$queryDir$subjectDir/$sequenceId-$subjectSetId.list";
			foreach (split ",", $seqToSet->{$sequenceId}->{$subjectSetId})
			{
				print LIST "$_\n";
			}
			close(LIST);
			print LOG "done!\n";
			$seqToSetNumber++;
		}
	}
}

my $setToSetNumber = 0;
foreach my $querySetId (keys %$setToSet)
{
	my $queryDir;
	for (my $position = 0; $position < length($querySetId); $position += 2)
	{
		$queryDir .= "/q". substr($querySetId,$position,2);
		until (-e "$commoncfg->{DATADIR}/alignments/setToSet$queryDir")
		{
			print LOG "Creating directory $commoncfg->{DATADIR}/alignments/setToSet$queryDir ... ";
			mkdir "$commoncfg->{DATADIR}/alignments/setToSet$queryDir";
			print LOG "done!\n";
		}
	}
	foreach my $subjectSetId (keys %{$setToSet->{$querySetId}})
	{

		my $subjectDir;
		for (my $position = 0; $position < length($subjectSetId); $position += 2)
		{
			$subjectDir .= "/s". substr($subjectSetId,$position,2);
			until (-e "$commoncfg->{DATADIR}/alignments/setToSet$queryDir$subjectDir")
			{
				print LOG "Creating directory $commoncfg->{DATADIR}/alignments/setToSet$queryDir$subjectDir ... ";
				mkdir "$commoncfg->{DATADIR}/alignments/setToSet$queryDir$subjectDir";
				print LOG "done!\n";
			}
		}
		
		until (-e "$commoncfg->{DATADIR}/alignments/setToSet$queryDir$subjectDir/$querySetId-$subjectSetId.list")
		{
			print LOG "Creating alignment setToSet list $commoncfg->{DATADIR}/alignments/setToSet$queryDir$subjectDir/$querySetId-$subjectSetId.list ... ";
			open (LIST,">$commoncfg->{DATADIR}/alignments/setToSet$queryDir$subjectDir/$querySetId-$subjectSetId.list") or die "can't open file: $commoncfg->{DATADIR}/alignments/setToSet$queryDir$subjectDir/$querySetId-$subjectSetId.list";
			foreach (split ",", $setToSet->{$querySetId}->{$subjectSetId})
			{
				print LIST "$_\n";
			}
			close(LIST);
			print LOG "done!\n";
			$setToSetNumber++;
		}
	}
}

print LOG "$seqToSeqNumber seqToSeq alignments exported.\n";
print LOG "$seqToSetNumber seqToSet lists exported.\n";
print LOG "$setToSetNumber setToSet lists exported.\n";
close(LOG);
