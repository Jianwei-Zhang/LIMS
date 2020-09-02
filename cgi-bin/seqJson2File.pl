#!/usr/bin/perl -w
use strict;
use JSON::XS; #JSON::XS is recommended to be installed for handling JSON string of big size 
use DBI;
use lib "lib/";
use lib "lib/pangu";
use pangu;

## this script is to convert sequences in Matrix into files

my $json = JSON::XS->new->allow_nonref;
my $commoncfg = readConfig("main.conf");
my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});

## sequences will be saved to $commoncfg->{DATADIR}/sequences
until (-e "$commoncfg->{DATADIR}/sequences")
{
	print "Creating directory $commoncfg->{DATADIR}/sequences ... ";
	mkdir "$commoncfg->{DATADIR}/sequences";
	print "done!\n";
}

my $seqNumber = 0;
my $getSequences = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence'");
$getSequences->execute();
while(my @getSequences = $getSequences->fetchrow_array())
{
	$seqNumber++;
	my $sequenceDetails = decode_json $getSequences[8];
	$sequenceDetails->{'sequence'} = '' unless (exists $sequenceDetails->{'sequence'});
	my $seqDir;
	for (my $position = 0; $position < length($getSequences[4]); $position += 2)
	{
		$seqDir .= "/g". substr($getSequences[4],$position,2);
		until (-e "$commoncfg->{DATADIR}/sequences$seqDir")
		{
			print "Creating directory $commoncfg->{DATADIR}/sequences$seqDir ... ";
			mkdir "$commoncfg->{DATADIR}/sequences$seqDir";
			print "done!\n";
		}
	}

	my $seqFile;
	for (my $position = 0; $position < length($getSequences[0]); $position += 2)
	{
		$seqFile .= "/s". substr($getSequences[0],$position,2);
		if (length ($seqFile) % 4 == 0)
		{
			until (-e "$commoncfg->{DATADIR}/sequences$seqDir$seqFile")
			{
				print "Creating directory $commoncfg->{DATADIR}/sequences$seqDir$seqFile ... ";
				mkdir "$commoncfg->{DATADIR}/sequences$seqDir$seqFile";
				print "done!\n";
			}
		}
	}

	until (-e "$commoncfg->{DATADIR}/sequences$seqDir$seqFile.fa")
	{
		print "Exporting No. $seqNumber sequence '$getSequences[4]-$getSequences[0]' ... ";
		open (SEQ,">$commoncfg->{DATADIR}/sequences$seqDir$seqFile.fa") or die "can't open file: $commoncfg->{DATADIR}/sequences$seqDir$seqFile.fa";
		print SEQ ">$getSequences[0]\n";
		print SEQ &multiLineSeq($sequenceDetails->{'sequence'},80);
		close(SEQ);
		print "done!\n";
	}

	print "Updating JSON of sequence '$getSequences[4]-$getSequences[0]' ... ";	
	$sequenceDetails->{'sequence'} = "sequences$seqDir$seqFile.fa";
	my $seqDetailsEncoded = $json->encode($sequenceDetails);
	my $updateSequence = $dbh->prepare("UPDATE matrix SET note = ? WHERE id = ?");
	$updateSequence->execute($seqDetailsEncoded,$getSequences[0]);
	print "done!\n";
}
