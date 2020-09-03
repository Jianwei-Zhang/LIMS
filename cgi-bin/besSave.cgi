#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use JSON::XS;
use Bio::SeqIO;
use DBI;
use lib "lib/";
use lib "lib/pangu";
use pangu;
use user;
use userCookie;

my $userCookie = new userCookie;
my $userId = (cookie('cid')) ? $userCookie->checkCookie(cookie('cid')) : 0;
exit if (!$userId);

my $user = new user;
my $userDetail = $user->getAllFieldsWithUserId($userId);
my $userName = $userDetail->{"userName"};

my $commoncfg = readConfig("main.conf");
my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});
my $libraryId = param ('libraryId') || '';
my $besFile = upload ('besFile');
my $besFilePath = param ('besFilePath') || '';
my $replace = param ('replace') || '0';
my $infile = "$commoncfg->{TMPDIR}/$$.bes";
my $json = JSON::XS->new->allow_nonref;
my $seqDir;
$seqDir->{'NA'} = 0;
$seqDir->{'f'} = 1;
$seqDir->{'r'} = 2;

print header;

if($besFile || $besFilePath)
{
	my $pid = fork();
	if ($pid) {
		print <<END;
<script>
parent.closeDialog();
parent.refresh("menu");
</script>	
END
	}
	elsif($pid == 0){
		close (STDOUT);
		if($besFilePath)
		{
			$infile = $besFilePath;
		}
		else
		{
			open (FILE, ">$infile");
			while (read ($besFile, my $Buffer, 1024)) {
				print FILE $Buffer;
			}
			close FILE;
		}		
		my $library=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$library->execute($libraryId);
		my @library = $library->fetchrow_array();
		if($replace > 0)
		{
			my $deleteBes=$dbh->do("DELETE FROM matrix WHERE container LIKE 'sequence' AND o = 98 AND x = $libraryId");	
		}
		#loading sequence
		my @allSequenceId;
		my $in = Bio::SeqIO->new(-file => $infile,
								-format => 'Fasta');
		while ( my $seq = $in->next_seq() )
		{
			my $sequenceDetails;
			$sequenceDetails->{'id'} = $seq->id;
			$sequenceDetails->{'description'} = $seq->desc || '';
			$sequenceDetails->{'sequence'} = '';
			$sequenceDetails->{'gapList'} = '';
			my $sequence = $seq->seq;
			$sequence =~ tr/a-zA-Z/N/c; #replace nonword characters.
			my $cloneName = $seq->id;
			my $seqDirNumber = 0;
			if($cloneName =~ /(.+)\.(.*)$/)
			{
				$seqDirNumber = $seqDir->{$2};
				$cloneName = $1;
			}
			$cloneName =~ /(\d+)(\D+)(\d+)$/;
			my $plateName =  sprintf "%0*d", 4, $1;
			my $row =  uc $2;
			my $col =  sprintf "%0*d", 2, $3;
			$cloneName = "$library[2]$plateName$row$col";
			my $seqLength = $seq->length();
			my $sequenceDetailsEncoded = $json->encode($sequenceDetails);
			my $insertSequence=$dbh->prepare("INSERT INTO matrix VALUES ('', 'sequence', ?, 98, ?, ?, ?, 0, ?, ?, NOW())");
			$insertSequence->execute($cloneName,$libraryId,$seqLength,$seqDirNumber,$sequenceDetailsEncoded,$userName);
			my $sequenceId = $dbh->{mysql_insertid};
			my $seqDir;
			for (my $position = 0; $position < length($libraryId); $position += 2)
			{
				$seqDir .= "/l". substr($libraryId,$position,2);
				until (-e "$commoncfg->{DATADIR}/sequences$seqDir")
				{
					mkdir "$commoncfg->{DATADIR}/sequences$seqDir";
				}
			}

			my $seqFile;
			for (my $position = 0; $position < length($sequenceId); $position += 2)
			{
				$seqFile .= "/s". substr($sequenceId,$position,2);
				if (length ($seqFile) % 4 == 0)
				{
					until (-e "$commoncfg->{DATADIR}/sequences$seqDir$seqFile")
					{
						mkdir "$commoncfg->{DATADIR}/sequences$seqDir$seqFile";
					}
				}
			}

			until (-e "$commoncfg->{DATADIR}/sequences$seqDir$seqFile.fa")
			{
				open (SEQ,">$commoncfg->{DATADIR}/sequences$seqDir$seqFile.fa") or die "can't open file: $commoncfg->{DATADIR}/sequences$seqDir$seqFile.fa";
				print SEQ ">$sequenceId\n";
				print SEQ &multiLineSeq($sequence,80);
				close(SEQ);
			}

			$sequenceDetails->{'sequence'} = "sequences$seqDir$seqFile.fa";
			$sequenceDetailsEncoded = $json->encode($sequenceDetails);
			my $updateSequence = $dbh->prepare("UPDATE matrix SET note = ? WHERE id = ?");
			$updateSequence->execute($sequenceDetailsEncoded,$sequenceId);
		}
		unlink ($infile) if (!$besFilePath);
		exit 0;
	}
	else{
		die "couldn't fork: $!\n";
	} 
}
else
{
	print <<END;
<script>
	parent.errorPop("No BES file found!");
</script>
END
	exit;
}
