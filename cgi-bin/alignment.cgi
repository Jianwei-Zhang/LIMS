#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use JSON::XS; #JSON::XS is recommended to be installed for handling JSON string of big size 
use Bio::SeqIO;
use DBI;
use lib "lib/";
use lib "lib/pangu";
use pangu;
use user;
use config;
use userConfig;
use userCookie;

my $userCookie = new userCookie;
my $userId = (cookie('cid')) ? $userCookie->checkCookie(cookie('cid')) : 0;
exit if (!$userId);

my $commoncfg = readConfig("main.conf");
my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});

my $user = new user;
my $config = new config;
my $userConfig = new userConfig;
my $author = $config->getFieldValueWithFieldName('AUTHOR');
my $siteName = $config->getFieldValueWithFieldName("SITENAME");
my $userDetail = $user->getAllFieldsWithUserId($userId);
my $userName = $userDetail->{"userName"};
my $userEmail = $userConfig->getFieldValueWithUserIdAndFieldName($userId,"email") if ($userId);
my $userFullName = $userConfig->getFieldValueWithUserIdAndFieldName($userId,"firstName")." ".$userConfig->getFieldValueWithUserIdAndFieldName($userId,"lastName") if ($userId);

my $alignEngineList;
$alignEngineList->{'blastn'} = "blast+/bin/blastn";
$alignEngineList->{'BLAT'} = "blat";
my $windowmasker = 'blast+/bin/windowmasker';
my $makeblastdb = 'blast+/bin/makeblastdb';
my $numThreads = 16;

my $queryGenomeId = param ('queryGenomeId') || '';
my $subjectGenomeId = param ('subjectGenomeId') || '';
my $identityAlignment = param ('identityAlignment') || $userConfig->getFieldValueWithUserIdAndFieldName($userId,"SEQTOGNMIDENTITY");
my $minOverlapAlignment = param ('minOverlapAlignment') || $userConfig->getFieldValueWithUserIdAndFieldName($userId,"SEQTOGNMMINOVERLAP");
my $alignEngine = param ('alignEngine') || 'blastn';
my $speedyMode = param ('speedyMode') || '0';
my $checkGood = param ('checkGood') || '0';
my $task = param ('megablast') || 'blastn';
my $softMasking = param ('softMasking') || '0';
my $markRepeatRegion = param ('markRepeatRegion') || '0';
my $emailNotification = param ('emailNotification') || '0';

print header;

if($queryGenomeId && $subjectGenomeId)
{
	my $pid = fork();
	if ($pid) {
		print <<END;
<script>
	parent.closeDialog();
	parent.errorPop("It's running! This processing might take a while.");
</script>	
END
	}
	elsif($pid == 0){
		close (STDOUT);
		#connect to the mysql server
		my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});
		
		my $queryFile = "$commoncfg->{TMPDIR}/$queryGenomeId.$$.seq";
		my $sequenceLength;

		my $queryGenome=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$queryGenome->execute($queryGenomeId);
		my @queryGenome = $queryGenome->fetchrow_array();

		open (SEQALL,">$queryFile") or die "can't open file: $queryFile";
		if($queryGenome[1] eq 'library')
		{
			my $getClones = $dbh->prepare("SELECT * FROM clones WHERE sequenced > 0 AND libraryId = ?");
			$getClones->execute($queryGenomeId);
			while(my @getClones = $getClones->fetchrow_array())
			{
				my $getSequences = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o < 50 AND name LIKE ?");
				$getSequences->execute($getClones[1]);
				while(my @getSequences = $getSequences->fetchrow_array())
				{
					$sequenceLength->{$getSequences[0]} = $getSequences[5];
					my $sequenceDetails = decode_json $getSequences[8];
					$sequenceDetails->{'id'} = '' unless (exists $sequenceDetails->{'id'});
					$sequenceDetails->{'description'} = '' unless (exists $sequenceDetails->{'description'});
					$sequenceDetails->{'gapList'} = '' unless (exists $sequenceDetails->{'gapList'});
					my $sequence = 'ERROR: NO SEQUENCE FOUND! PLEASE CONTACT YOUR ADMINISTRATOR.';
					my $in = Bio::SeqIO->new(-file => "$commoncfg->{DATADIR}/$sequenceDetails->{'sequence'}",
											-format => 'Fasta');
					while ( my $seq = $in->next_seq() )
					{
						$sequence = $seq->seq;
					}

					print SEQALL ">$getSequences[0]\n$sequence\n";
				}
			}
		}
		if($queryGenome[1] eq 'genome')
		{
			my $getSequences = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 99 AND x = ?");
			$getSequences->execute($queryGenomeId);
			while(my @getSequences = $getSequences->fetchrow_array())
			{
				$sequenceLength->{$getSequences[0]} = $getSequences[5];
				my $sequenceDetails = decode_json $getSequences[8];
				$sequenceDetails->{'id'} = '' unless (exists $sequenceDetails->{'id'});
				$sequenceDetails->{'description'} = '' unless (exists $sequenceDetails->{'description'});
				$sequenceDetails->{'gapList'} = '' unless (exists $sequenceDetails->{'gapList'});
				my $sequence = 'ERROR: NO SEQUENCE FOUND! PLEASE CONTACT YOUR ADMINISTRATOR.';
				my $in = Bio::SeqIO->new(-file => "$commoncfg->{DATADIR}/$sequenceDetails->{'sequence'}",
										-format => 'Fasta');
				while ( my $seq = $in->next_seq() )
				{
					$sequence = $seq->seq;
				}
				print SEQALL ">$getSequences[0]\n$sequence\n";
			}
		}
		close(SEQALL);

		my $subjectGenome=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$subjectGenome->execute($subjectGenomeId);
		my @subjectGenome = $subjectGenome->fetchrow_array();
		my $subjectFile = "";
		if($queryGenomeId eq $subjectGenomeId)
		{
			$subjectFile = $queryFile;
		}
		else
		{
			$subjectFile = "$commoncfg->{TMPDIR}/$subjectGenomeId.$$.seq";

			open (SEQALL,">$subjectFile") or die "can't open file: $subjectFile";
			if($subjectGenome[1] eq 'library')
			{
				my $getClones = $dbh->prepare("SELECT * FROM clones WHERE sequenced > 0 AND libraryId = ?");
				$getClones->execute($subjectGenomeId);
				while(my @getClones = $getClones->fetchrow_array())
				{
					my $getSequences = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o < 50 AND name LIKE ?");
					$getSequences->execute($getClones[1]);
					while(my @getSequences = $getSequences->fetchrow_array())
					{
						$sequenceLength->{$getSequences[0]} = $getSequences[5];
						my $sequenceDetails = decode_json $getSequences[8];
						$sequenceDetails->{'id'} = '' unless (exists $sequenceDetails->{'id'});
						$sequenceDetails->{'description'} = '' unless (exists $sequenceDetails->{'description'});
						$sequenceDetails->{'gapList'} = '' unless (exists $sequenceDetails->{'gapList'});
						my $sequence = 'ERROR: NO SEQUENCE FOUND! PLEASE CONTACT YOUR ADMINISTRATOR.';
						my $in = Bio::SeqIO->new(-file => "$commoncfg->{DATADIR}/$sequenceDetails->{'sequence'}",
												-format => 'Fasta');
						while ( my $seq = $in->next_seq() )
						{
							$sequence = $seq->seq;
						}
						print SEQALL ">$getSequences[0]\n$sequence\n";
					}
				}
			}
			if($subjectGenome[1] eq 'genome')
			{
				my $getSequences = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 99 AND x = ?");
				$getSequences->execute($subjectGenomeId);
				while(my @getSequences = $getSequences->fetchrow_array())
				{
					$sequenceLength->{$getSequences[0]} = $getSequences[5];
					my $sequenceDetails = decode_json $getSequences[8];
					$sequenceDetails->{'id'} = '' unless (exists $sequenceDetails->{'id'});
					$sequenceDetails->{'description'} = '' unless (exists $sequenceDetails->{'description'});
					$sequenceDetails->{'gapList'} = '' unless (exists $sequenceDetails->{'gapList'});
					my $sequence = 'ERROR: NO SEQUENCE FOUND! PLEASE CONTACT YOUR ADMINISTRATOR.';
					my $in = Bio::SeqIO->new(-file => "$commoncfg->{DATADIR}/$sequenceDetails->{'sequence'}",
											-format => 'Fasta');
					while ( my $seq = $in->next_seq() )
					{
						$sequence = $seq->seq;
					}
					print SEQALL ">$getSequences[0]\n$sequence\n";
				}
			}
			close(SEQALL);
		}

		if($alignEngine eq 'blastn')
		{
			if($softMasking)
			{
				system( "$windowmasker -in $subjectFile -infmt fasta -mk_counts -parse_seqids -out $subjectFile.mask.counts" );
				system( "$windowmasker -in $subjectFile -infmt fasta -ustat $subjectFile.mask.counts -outfmt maskinfo_asn1_bin -parse_seqids -out $subjectFile.mask.asnb" );
				system( "$makeblastdb -in $subjectFile -inputtype fasta -dbtype nucl -parse_seqids -mask_data $subjectFile.mask.asnb" );
			}
			else
			{
				system( "$makeblastdb -in $subjectFile -dbtype nucl" );
			}
		}

		my $goodSequenceId;
		if($alignEngine eq 'blastn')
		{
			if($softMasking)
			{
				open (CMD,"$alignEngineList->{$alignEngine} -query $queryFile -task $task -db $subjectFile -db_soft_mask 30 -evalue 1e-200 -perc_identity $identityAlignment -num_threads $numThreads -outfmt 6 |") or die "can't open CMD: $!";
			}
			else
			{
				open (CMD,"$alignEngineList->{$alignEngine} -query $queryFile -task $task -db $subjectFile -evalue 1e-200 -perc_identity $identityAlignment -num_threads $numThreads -outfmt 6 |") or die "can't open CMD: $!";
			}
		}
		else
		{
			open (CMD,"$alignEngineList->{$alignEngine} $subjectFile $queryFile -out=blast8 -minIdentity=$identityAlignment |") or die "can't open CMD: $!";
		}
		while(<CMD>)
		{
			/^#/ and next;
			my @hit = split("\t",$_);
			next if($hit[0] eq $hit[1]);
			next if($hit[3] < $minOverlapAlignment);
			if($speedyMode)
			{

				my $deleteAlignmentFlag = 0;
				if($hit[0] < $hit[1])
				{
					unless(exists $goodSequenceId->{$hit[0]}->{$hit[1]})
					{
						$goodSequenceId->{$hit[0]}->{$hit[1]} = 1;
						$deleteAlignmentFlag = 1;
					}
				}
				else
				{
					unless(exists $goodSequenceId->{$hit[1]}->{$hit[0]})
					{
						$goodSequenceId->{$hit[1]}->{$hit[0]} = 1;
						$deleteAlignmentFlag = 1;
					}
				}
				if($deleteAlignmentFlag)
				{
					my $deleteAlignmentA = $dbh->do("DELETE FROM alignment WHERE query = $hit[0] AND subject = $hit[1]");
					my $deleteAlignmentB = $dbh->do("DELETE FROM alignment WHERE query = $hit[1] AND subject = $hit[0]");
				}

				#write to alignment
				my $insertAlignmentA=$dbh->prepare("INSERT INTO alignment VALUES ('', '$alignEngine\_1e-200\_$identityAlignment\_$minOverlapAlignment', ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 0)");
				$insertAlignmentA->execute(@hit);

				#switch query and subject
				if($hit[8] < $hit[9])
				{
					my $exchange = $hit[8];
					$hit[8] = $hit[6];
					$hit[6] = $exchange;
					$exchange = $hit[9];
					$hit[9] = $hit[7];
					$hit[7] = $exchange;
					$exchange = $hit[1];
					$hit[1] = $hit[0];
					$hit[0] = $exchange;
				}
				else
				{
					my $exchange = $hit[8];
					$hit[8] = $hit[7];
					$hit[7] = $exchange;
					$exchange = $hit[9];
					$hit[9] = $hit[6];
					$hit[6] = $exchange;
					$exchange = $hit[1];
					$hit[1] = $hit[0];
					$hit[0] = $exchange;
				}

				my $insertAlignmentB=$dbh->prepare("INSERT INTO alignment VALUES ('', '$alignEngine\_1e-200\_$identityAlignment\_$minOverlapAlignment', ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 0)");
				$insertAlignmentB->execute(@hit);
			}
			else
			{
				my $rerunBlastTwo = 0;
				if($hit[0] < $hit[1])
				{
					unless(exists $goodSequenceId->{$hit[0]}->{$hit[1]})
					{
						$goodSequenceId->{$hit[0]}->{$hit[1]} = 1;
						$rerunBlastTwo = 1;
					}
				}
				else
				{
					unless(exists $goodSequenceId->{$hit[1]}->{$hit[0]})
					{
						$goodSequenceId->{$hit[1]}->{$hit[0]} = 1;
						$rerunBlastTwo = 1;
					}
				}
				if($rerunBlastTwo)
				{
					my $deleteAlignmentA = $dbh->do("DELETE FROM alignment WHERE query = $hit[0] AND subject = $hit[1]");
					my $deleteAlignmentB = $dbh->do("DELETE FROM alignment WHERE query = $hit[1] AND subject = $hit[0]");

					unless(-e "$commoncfg->{TMPDIR}/$hit[0].$$.seq")
					{
						my $getSequenceA = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
						$getSequenceA->execute($hit[0]);
						my @getSequenceA =  $getSequenceA->fetchrow_array();
						open (SEQA,">$commoncfg->{TMPDIR}/$hit[0].$$.seq") or die "can't open file: $commoncfg->{TMPDIR}/$hit[0].$$.seq";
						my $sequenceDetailsA = decode_json $getSequenceA[8];
						$sequenceDetailsA->{'id'} = '' unless (exists $sequenceDetailsA->{'id'});
						$sequenceDetailsA->{'description'} = '' unless (exists $sequenceDetailsA->{'description'});
						$sequenceDetailsA->{'gapList'} = '' unless (exists $sequenceDetailsA->{'gapList'});
						my $sequenceA = 'ERROR: NO SEQUENCE FOUND! PLEASE CONTACT YOUR ADMINISTRATOR.';
						my $in = Bio::SeqIO->new(-file => "$commoncfg->{DATADIR}/$sequenceDetailsA->{'sequence'}",
												-format => 'Fasta');
						while ( my $seq = $in->next_seq() )
						{
							$sequenceA = $seq->seq;
						}
						print SEQA ">$getSequenceA[0]\n$sequenceA\n";
						close(SEQA);
					}
					unless(-e "$commoncfg->{TMPDIR}/$hit[1].$$.seq")
					{
						my $getSequenceB = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
						$getSequenceB->execute($hit[1]);
						my @getSequenceB =  $getSequenceB->fetchrow_array();
						open (SEQB,">$commoncfg->{TMPDIR}/$hit[1].$$.seq") or die "can't open file: $commoncfg->{TMPDIR}/$hit[1].$$.seq";
						my $sequenceDetailsB = decode_json $getSequenceB[8];
						$sequenceDetailsB->{'id'} = '' unless (exists $sequenceDetailsB->{'id'});
						$sequenceDetailsB->{'description'} = '' unless (exists $sequenceDetailsB->{'description'});
						$sequenceDetailsB->{'gapList'} = '' unless (exists $sequenceDetailsB->{'gapList'});
						my $sequenceB = 'ERROR: NO SEQUENCE FOUND! PLEASE CONTACT YOUR ADMINISTRATOR.';
						my $in = Bio::SeqIO->new(-file => "$commoncfg->{DATADIR}/$sequenceDetailsB->{'sequence'}",
												-format => 'Fasta');
						while ( my $seq = $in->next_seq() )
						{
							$sequenceB = $seq->seq;
						}
						print SEQB ">$getSequenceB[0]\n$sequenceB\n";
						close(SEQB);
					}
					my @alignments;
					my $goodOverlap = ($checkGood) ? 0 : 1;
					open (CMDA,"$alignEngineList->{'blastn'} -query $commoncfg->{TMPDIR}/$hit[0].$$.seq -subject $commoncfg->{TMPDIR}/$hit[1].$$.seq -dust no -evalue 1e-200 -perc_identity $identityAlignment -outfmt 6 |") or die "can't open CMD: $!";
					while(<CMDA>)
					{
						/^#/ and next;
						my @detailedHit = split("\t",$_);
						if($detailedHit[3] >= $minOverlapAlignment)
						{
							push @alignments, $_;
							if($detailedHit[6] == 1 || $detailedHit[7] == $sequenceLength->{$detailedHit[0]})
							{
								$goodOverlap = 1;
							}
							#switch query and subject
							if($detailedHit[8] < $detailedHit[9])
							{
								my $exchange = $detailedHit[8];
								$detailedHit[8] = $detailedHit[6];
								$detailedHit[6] = $exchange;
								$exchange = $detailedHit[9];
								$detailedHit[9] = $detailedHit[7];
								$detailedHit[7] = $exchange;
								$exchange = $detailedHit[1];
								$detailedHit[1] = $detailedHit[0];
								$detailedHit[0] = $exchange;
							}
							else
							{
								my $exchange = $detailedHit[8];
								$detailedHit[8] = $detailedHit[7];
								$detailedHit[7] = $exchange;
								$exchange = $detailedHit[9];
								$detailedHit[9] = $detailedHit[6];
								$detailedHit[6] = $exchange;
								$exchange = $detailedHit[1];
								$detailedHit[1] = $detailedHit[0];
								$detailedHit[0] = $exchange;
							}

							if($detailedHit[6] == 1 || $detailedHit[7] == $sequenceLength->{$detailedHit[0]})
							{
								$goodOverlap = 1;
							}
							my $reverseBlast = join "\t",@detailedHit;
							push @alignments, $reverseBlast;							
						}									
					}
					close(CMDA);
					if($goodOverlap)
					{
						foreach (@alignments)
						{
							my @detailedHit = split("\t",$_);
							#write to alignment
							my $insertAlignment=$dbh->prepare("INSERT INTO alignment VALUES ('', '$alignEngine\_1e-200\_$identityAlignment\_$minOverlapAlignment', ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 0)");
							$insertAlignment->execute(@detailedHit);
						}
					}
				}
			}
		}
		close(CMD);
		unlink("$queryFile");
		unlink("$subjectFile");
		unlink("$subjectFile.nhr");
		unlink("$subjectFile.nin");
		unlink("$subjectFile.nsq");
		`rm $commoncfg->{TMPDIR}/*.aln.html`; #delete cached files

		foreach my $queryId (keys %$goodSequenceId)
		{
			unlink("$commoncfg->{TMPDIR}/$queryId.$$.seq");
			foreach my $subjectId (keys %{$goodSequenceId->{$queryId}})
			{
				unlink("$commoncfg->{TMPDIR}/$subjectId.$$.seq");
			}
		}

		if ($markRepeatRegion)
		{
			my $converageCutoff = 0.95;
			my $todo = 1;
			do{
				$todo = 0;
				my $lastQuery = 0;
				my $lastQend = 0;
				my $lastAlignmentLength = 0;
				my $lastAlignmentId = 0;
				my $getAlignment = $dbh->prepare("SELECT * FROM alignment WHERE program = '$alignEngine\_1e-200\_$identityAlignment\_$minOverlapAlignment' AND align_length < 5000 AND hidden = 0 ORDER BY query,q_start,q_end");
				$getAlignment->execute();
				while (my @getAlignment = $getAlignment->fetchrow_array())
				{
					if($lastQuery == $getAlignment[2])
					{
						if($lastQend > $getAlignment[8])
						{
							my $converage = $lastQend - $getAlignment[8];
							if ( $converage/$lastAlignmentLength >= $converageCutoff)
							{
								my $hideLastAlignment=$dbh->do("UPDATE alignment SET hidden = 1 WHERE id = $lastAlignmentId");
								$todo = 1;
							}
							if ( $converage/$getAlignment[5] >= $converageCutoff)
							{
								my $hideAlignment=$dbh->do("UPDATE alignment SET hidden = 1 WHERE id = $getAlignment[0]");
								$todo = 1;
							}
						}
					}
					$lastQuery = $getAlignment[2];
					$lastQend = $getAlignment[9];
					$lastAlignmentLength = $getAlignment[5];
					$lastAlignmentId = $getAlignment[0];
				}
			} while($todo);
		}
		
		if($emailNotification)
		{
			#email to user after alignment finishes.
			open(MAIL,"|/usr/sbin/sendmail -t -oi");
			print MAIL "To: $userEmail\n";
			print MAIL "From: $author\n";
			print MAIL "Subject: Alignment Successfully Completed\n\n";
			print MAIL <<eof;
Dear $userFullName ($userName),

Your alignment job has successfully completed.

$queryGenome[2] vs $subjectGenome[2] (-minOverlap $minOverlapAlignment -perc_identity $identityAlignment)

Best regards,
Dev Team
$siteName
eof
		}
	}
	else{
		die "couldn't fork: $!\n";
	} 
}
else
{
	print <<END;
<script>
	parent.errorPop("Please provide required information!");
</script>	
END
}