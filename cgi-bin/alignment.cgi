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
## alignments will be saved to $commoncfg->{DATADIR}/alignments
unless (-e "$commoncfg->{DATADIR}/alignments")
{
	mkdir "$commoncfg->{DATADIR}/alignments";
}
unless (-e "$commoncfg->{DATADIR}/alignments/seqToSeq")
{
	mkdir "$commoncfg->{DATADIR}/alignments/seqToSeq";
}
unless (-e "$commoncfg->{DATADIR}/alignments/seqToSet")
{
	mkdir "$commoncfg->{DATADIR}/alignments/seqToSet";
}
unless (-e "$commoncfg->{DATADIR}/alignments/setToSet")
{
	mkdir "$commoncfg->{DATADIR}/alignments/setToSet";
}

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
my $redo = param ('redo') || '0';
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
	if ($pid)
	{
		print <<END;
<script>
	parent.closeDialog();
	parent.errorPop("It's running! This processing might take a while.");
</script>	
END
	}
	elsif($pid == 0)
	{
		close (STDOUT);
		#connect to the mysql server
		my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});
		my $setId;
		my @queryIdList;
		my @subjectIdList;
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
					$setId->{$getSequences[0]} = $getSequences[4];
					push @queryIdList,$getSequences[0];
					$sequenceLength->{$getSequences[0]} = $getSequences[5];
					my $sequenceDetails = decode_json $getSequences[8];
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
				$setId->{$getSequences[0]} = $getSequences[4];
				push @queryIdList,$getSequences[0];
				$sequenceLength->{$getSequences[0]} = $getSequences[5];
				my $sequenceDetails = decode_json $getSequences[8];
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
			@subjectIdList = @queryIdList;
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
						$setId->{$getSequences[0]} = $getSequences[4];
						push @subjectIdList,$getSequences[0];
						$sequenceLength->{$getSequences[0]} = $getSequences[5];
						my $sequenceDetails = decode_json $getSequences[8];
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
					$setId->{$getSequences[0]} = $getSequences[4];
					push @subjectIdList,$getSequences[0];
					$sequenceLength->{$getSequences[0]} = $getSequences[5];
					my $sequenceDetails = decode_json $getSequences[8];
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

		if ($redo)
		{
			foreach my $querySequenceId (@queryIdList)
			{
				my $queryAsQueryDir;
				for (my $position = 0; $position < length($querySequenceId); $position += 2)
				{
					$queryAsQueryDir .= "/q". substr($querySequenceId,$position,2);
				}	
				my $queryAsSubjectDir;
				for (my $position = 0; $position < length($querySequenceId); $position += 2)
				{
					$queryAsSubjectDir .= "/s". substr($querySequenceId,$position,2);
				}
				foreach my $subjectSequenceId (@subjectIdList)
				{
					my $subjectAsQueryDir;
					for (my $position = 0; $position < length($subjectSequenceId); $position += 2)
					{
						$subjectAsQueryDir .= "/q". substr($subjectSequenceId,$position,2);
					}	
					my $subjectAsSubjectDir;
					for (my $position = 0; $position < length($subjectSequenceId); $position += 2)
					{
						$subjectAsSubjectDir .= "/s". substr($subjectSequenceId,$position,2);
					}	
					unlink("$commoncfg->{DATADIR}/alignments/seqToSeq$queryAsQueryDir$subjectAsSubjectDir/$querySequenceId-$subjectSequenceId.tbl") if (-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryAsQueryDir$subjectAsSubjectDir/$querySequenceId-$subjectSequenceId.tbl");
					unlink("$commoncfg->{DATADIR}/alignments/seqToSeq$subjectAsQueryDir$queryAsSubjectDir/$subjectSequenceId-$querySequenceId.tbl") if (-e "$commoncfg->{DATADIR}/alignments/seqToSeq$subjectAsQueryDir$queryAsSubjectDir/$subjectSequenceId-$querySequenceId.tbl");
				}
			}
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

		my $seqToSeq;
		my $seqToSet;
		my $seqToSetSwitched;
		my $setToSet;
		my $setToSetSwitched;

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
			chop;
			/^#/ and next;
			my @hit = split("\t",$_);
			$hit[12] = 0; #add a hidden column
			next if($hit[0] eq $hit[1]);
			next if($hit[3] < $minOverlapAlignment);

			my $queryDir;
			my $queryDirSwitched;
			my $subjectDir;
			my $subjectDirSwitched;			
			if (exists $seqToSeq->{$hit[0]}->{$hit[1]})
			{
				for (my $position = 0; $position < length($hit[0]); $position += 2)
				{
					$queryDir .= "/q". substr($hit[0],$position,2);
				}
				for (my $position = 0; $position < length($hit[1]); $position += 2)
				{
					$subjectDir .= "/s". substr($hit[1],$position,2);
				}
				for (my $position = 0; $position < length($hit[1]); $position += 2)
				{
					$queryDirSwitched .= "/q". substr($hit[1],$position,2);
				}
				for (my $position = 0; $position < length($hit[0]); $position += 2)
				{
					$subjectDirSwitched .= "/s". substr($hit[0],$position,2);
				}
				$seqToSeq->{$hit[0]}->{$hit[1]}++;
			}
			else # clean old data first
			{
				for (my $position = 0; $position < length($hit[0]); $position += 2)
				{
					$queryDir .= "/q". substr($hit[0],$position,2);
					unless (-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir")
					{
						mkdir "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir";
					}
				}
				for (my $position = 0; $position < length($hit[1]); $position += 2)
				{
					$subjectDir .= "/s". substr($hit[1],$position,2);
					unless (-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir")
					{
						mkdir "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir";
					}
				}
				open (ALN,">$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$hit[0]-$hit[1].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$hit[0]-$hit[1].tbl";
				print ALN "#$alignEngine\_1e-200\_$identityAlignment\_$minOverlapAlignment\n";
				print ALN "#query\tsubject\tperc_indentity\talign_length\tmismatches\tgaps\tq_start\tq_end\ts_start\ts_end\te_val\tbit_score\thidden\n";
				close(ALN);

				for (my $position = 0; $position < length($hit[1]); $position += 2)
				{
					$queryDirSwitched .= "/q". substr($hit[1],$position,2);
					unless (-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched")
					{
						mkdir "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched";
					}
				}
				for (my $position = 0; $position < length($hit[0]); $position += 2)
				{
					$subjectDirSwitched .= "/s". substr($hit[0],$position,2);
					unless (-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched")
					{
						mkdir "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched";
					}
				}
				open (ALN,">$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$hit[1]-$hit[0].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$hit[1]-$hit[0].tbl";
				print ALN "#$alignEngine\_1e-200\_$identityAlignment\_$minOverlapAlignment\n";
				print ALN "#query\tsubject\tperc_indentity\talign_length\tmismatches\tgaps\tq_start\tq_end\ts_start\ts_end\te_val\tbit_score\thidden\n";
				close(ALN);
				$seqToSeq->{$hit[0]}->{$hit[1]} = 1;
			}

			if($speedyMode)
			{
				if ($seqToSeq->{$hit[0]}->{$hit[1]} == 1) #check if this is the first hit
				{
					if (exists $seqToSet->{$hit[0]}->{$setId->{$hit[1]}})
					{
						$seqToSet->{$hit[0]}->{$setId->{$hit[1]}} .= ",alignments/seqToSeq$queryDir$subjectDir/$hit[0]-$hit[1].tbl";
					}
					else
					{
						$seqToSet->{$hit[0]}->{$setId->{$hit[1]}} = "alignments/seqToSeq$queryDir$subjectDir/$hit[0]-$hit[1].tbl";
					}
					if (exists $setToSet->{$setId->{$hit[0]}}->{$setId->{$hit[1]}})
					{
						$setToSet->{$setId->{$hit[0]}}->{$setId->{$hit[1]}} .= ",alignments/seqToSeq$queryDir$subjectDir/$hit[0]-$hit[1].tbl";
					}
					else
					{
						$setToSet->{$setId->{$hit[0]}}->{$setId->{$hit[1]}} = "alignments/seqToSeq$queryDir$subjectDir/$hit[0]-$hit[1].tbl";
					}

					if (exists $seqToSetSwitched->{$hit[1]}->{$setId->{$hit[0]}})
					{
						$seqToSetSwitched->{$hit[1]}->{$setId->{$hit[0]}} .= ",alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$hit[1]-$hit[0].tbl";
					}
					else
					{
						$seqToSetSwitched->{$hit[1]}->{$setId->{$hit[0]}} = "alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$hit[1]-$hit[0].tbl";
					}
					if (exists $setToSetSwitched->{$setId->{$hit[1]}}->{$setId->{$hit[0]}})
					{
						$setToSetSwitched->{$setId->{$hit[1]}}->{$setId->{$hit[0]}} .= ",alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$hit[1]-$hit[0].tbl";
			
					}
					else
					{
						$setToSetSwitched->{$setId->{$hit[1]}}->{$setId->{$hit[0]}} = "alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$hit[1]-$hit[0].tbl";
					}
				}
				#write to alignment
				open (ALN,">>$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$hit[0]-$hit[1].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$hit[0]-$hit[1].tbl";
				print ALN join "\t", @hit;
				print ALN "\n";
				close(ALN);

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
				#write to alignment
				open (ALN,">>$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$hit[0]-$hit[1].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$hit[0]-$hit[1].tbl";
				print ALN join "\t", @hit;
				print ALN "\n";
				close(ALN);
			}
			else
			{
				if ($seqToSeq->{$hit[0]}->{$hit[1]} == 1) #check if this is the first hit
				{
					#prepare sequences for rerun alignment
					unless(-e "$commoncfg->{TMPDIR}/$hit[0].$$.seq")
					{
						my $getSequenceA = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
						$getSequenceA->execute($hit[0]);
						my @getSequenceA =  $getSequenceA->fetchrow_array();
						open (SEQA,">$commoncfg->{TMPDIR}/$hit[0].$$.seq") or die "can't open file: $commoncfg->{TMPDIR}/$hit[0].$$.seq";
						my $sequenceDetailsA = decode_json $getSequenceA[8];
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
					my @alignmentsSwitched;
					my $goodOverlap = ($checkGood) ? 0 : 1;
					open (CMDA,"$alignEngineList->{$alignEngine} -query $commoncfg->{TMPDIR}/$hit[0].$$.seq -subject $commoncfg->{TMPDIR}/$hit[1].$$.seq -dust no -evalue 1e-200 -perc_identity $identityAlignment -outfmt 6 |") or die "can't open CMD: $!";
					while(<CMDA>)
					{
						chop;
						/^#/ and next;
						my @detailedHit = split("\t",$_);
						$detailedHit[12] = 0; #add a hidden column
						if($detailedHit[3] >= $minOverlapAlignment)
						{
							my $alignment = join "\t",@detailedHit;
							push @alignments, $alignment;
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
							$alignment = join "\t",@detailedHit;
							push @alignmentsSwitched, $alignment;							
						}									
					}
					close(CMDA);
					if($goodOverlap)
					{
						my $lastQuery = '';
						my $lastSubject = '';
						foreach (@alignments)
						{
							my @detailedHit = split("\t",$_);
							if($lastQuery != $detailedHit[0] && $lastSubject != $detailedHit[1])
							{
								if (exists $seqToSet->{$detailedHit[0]}->{$setId->{$detailedHit[1]}})
								{
									$seqToSet->{$detailedHit[0]}->{$setId->{$detailedHit[1]}} .= ",alignments/seqToSeq$queryDir$subjectDir/$detailedHit[0]-$detailedHit[1].tbl";
								}
								else
								{
									$seqToSet->{$detailedHit[0]}->{$setId->{$detailedHit[1]}} = "alignments/seqToSeq$queryDir$subjectDir/$detailedHit[0]-$detailedHit[1].tbl";
								}
								if (exists $setToSet->{$setId->{$detailedHit[0]}}->{$setId->{$detailedHit[1]}})
								{
									$setToSet->{$setId->{$detailedHit[0]}}->{$setId->{$detailedHit[1]}} .= ",alignments/seqToSeq$queryDir$subjectDir/$detailedHit[0]-$detailedHit[1].tbl";
								}
								else
								{
									$setToSet->{$setId->{$detailedHit[0]}}->{$setId->{$detailedHit[1]}} = "alignments/seqToSeq$queryDir$subjectDir/$detailedHit[0]-$detailedHit[1].tbl";
								}
							}
							#write to alignment
							open (ALN,">>$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$detailedHit[0]-$detailedHit[1].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$detailedHit[0]-$detailedHit[1].tbl";
							print ALN join "\t", @detailedHit;
							print ALN "\n";
							close(ALN);
							$lastQuery = $detailedHit[0];
							$lastSubject = $detailedHit[1];
						}
						$lastQuery = '';
						$lastSubject = '';						
						foreach (@alignmentsSwitched)
						{
							my @detailedHit = split("\t",$_);
							if($lastQuery != $detailedHit[0] && $lastSubject != $detailedHit[1])
							{
								if (exists $seqToSetSwitched->{$detailedHit[0]}->{$setId->{$detailedHit[1]}})
								{
									$seqToSetSwitched->{$detailedHit[0]}->{$setId->{$detailedHit[1]}} .= ",alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$detailedHit[0]-$detailedHit[1].tbl";
								}
								else
								{
									$seqToSetSwitched->{$detailedHit[0]}->{$setId->{$detailedHit[1]}} = "alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$detailedHit[0]-$detailedHit[1].tbl";
								}
								if (exists $setToSetSwitched->{$setId->{$detailedHit[0]}}->{$setId->{$detailedHit[1]}})
								{
									$setToSetSwitched->{$setId->{$detailedHit[0]}}->{$setId->{$detailedHit[1]}} .= ",alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$detailedHit[0]-$detailedHit[1].tbl";
								}
								else
								{
									$setToSetSwitched->{$setId->{$detailedHit[0]}}->{$setId->{$detailedHit[1]}} = "alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$detailedHit[0]-$detailedHit[1].tbl";
								}
							}
							#write to alignment
							open (ALN,">>$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$detailedHit[0]-$detailedHit[1].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$detailedHit[0]-$detailedHit[1].tbl";
							print ALN join "\t", @detailedHit;
							print ALN "\n";
							close(ALN);
							$lastQuery = $detailedHit[0];
							$lastSubject = $detailedHit[1];
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

		foreach my $queryId (keys %$seqToSeq)
		{
			unlink("$commoncfg->{TMPDIR}/$queryId.$$.seq");
			foreach my $subjectId (keys %{$seqToSeq->{$queryId}})
			{
				unlink("$commoncfg->{TMPDIR}/$subjectId.$$.seq");
			}
		}

		foreach my $sequenceId (keys %$seqToSet)
		{
			my $queryDirLocal;
			for (my $position = 0; $position < length($sequenceId); $position += 2)
			{
				$queryDirLocal .= "/q". substr($sequenceId,$position,2);
				unless (-e "$commoncfg->{DATADIR}/alignments/seqToSet$queryDirLocal")
				{
					mkdir "$commoncfg->{DATADIR}/alignments/seqToSet$queryDirLocal";
				}
			}
			foreach my $subjectSetId (keys %{$seqToSet->{$sequenceId}})
			{
				my $subjectDirLocal;
				for (my $position = 0; $position < length($subjectSetId); $position += 2)
				{
					$subjectDirLocal .= "/s". substr($subjectSetId,$position,2);
					unless (-e "$commoncfg->{DATADIR}/alignments/seqToSet$queryDirLocal$subjectDirLocal")
					{
						mkdir "$commoncfg->{DATADIR}/alignments/seqToSet$queryDirLocal$subjectDirLocal";
					}
				}
				open (LIST,">$commoncfg->{DATADIR}/alignments/seqToSet$queryDirLocal$subjectDirLocal/$sequenceId-$subjectSetId.list") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSet$queryDirLocal$subjectDirLocal/$sequenceId-$subjectSetId.list";
				foreach (split ",", $seqToSet->{$sequenceId}->{$subjectSetId})
				{
					print LIST "$_\n";
				}
				close(LIST);
			}
		}

		foreach my $querySetId (keys %$setToSet)
		{
			my $queryDirLocal;
			for (my $position = 0; $position < length($querySetId); $position += 2)
			{
				$queryDirLocal .= "/q". substr($querySetId,$position,2);
				unless (-e "$commoncfg->{DATADIR}/alignments/setToSet$queryDirLocal")
				{
					mkdir "$commoncfg->{DATADIR}/alignments/setToSet$queryDirLocal";
				}
			}
			foreach my $subjectSetId (keys %{$setToSet->{$querySetId}})
			{
				my $subjectDirLocal;
				for (my $position = 0; $position < length($subjectSetId); $position += 2)
				{
					$subjectDirLocal .= "/s". substr($subjectSetId,$position,2);
					unless (-e "$commoncfg->{DATADIR}/alignments/setToSet$queryDirLocal$subjectDirLocal")
					{
						mkdir "$commoncfg->{DATADIR}/alignments/setToSet$queryDirLocal$subjectDirLocal";
					}
				}
				open (LIST,">$commoncfg->{DATADIR}/alignments/setToSet$queryDirLocal$subjectDirLocal/$querySetId-$subjectSetId.list") or die "can't open file: $commoncfg->{DATADIR}/alignments/setToSet$queryDirLocal$subjectDirLocal/$querySetId-$subjectSetId.list";
				foreach (split ",", $setToSet->{$querySetId}->{$subjectSetId})
				{
					print LIST "$_\n";
				}
				close(LIST);
			}
		}

		foreach my $sequenceId (keys %$seqToSetSwitched)
		{
			my $queryDirLocal;
			for (my $position = 0; $position < length($sequenceId); $position += 2)
			{
				$queryDirLocal .= "/q". substr($sequenceId,$position,2);
				unless (-e "$commoncfg->{DATADIR}/alignments/seqToSet$queryDirLocal")
				{
					mkdir "$commoncfg->{DATADIR}/alignments/seqToSet$queryDirLocal";
				}
			}
			foreach my $subjectSetId (keys %{$seqToSetSwitched->{$sequenceId}})
			{
				my $subjectDirLocal;
				for (my $position = 0; $position < length($subjectSetId); $position += 2)
				{
					$subjectDirLocal .= "/s". substr($subjectSetId,$position,2);
					unless (-e "$commoncfg->{DATADIR}/alignments/seqToSet$queryDirLocal$subjectDirLocal")
					{
						mkdir "$commoncfg->{DATADIR}/alignments/seqToSet$queryDirLocal$subjectDirLocal";
					}
				}
				open (LIST,">$commoncfg->{DATADIR}/alignments/seqToSet$queryDirLocal$subjectDirLocal/$sequenceId-$subjectSetId.list") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSet$queryDirLocal$subjectDirLocal/$sequenceId-$subjectSetId.list";
				foreach (split ",", $seqToSetSwitched->{$sequenceId}->{$subjectSetId})
				{
					print LIST "$_\n";
				}
				close(LIST);
			}
		}

		foreach my $querySetId (keys %$setToSetSwitched)
		{
			my $queryDirLocal;
			for (my $position = 0; $position < length($querySetId); $position += 2)
			{
				$queryDirLocal .= "/q". substr($querySetId,$position,2);
				unless (-e "$commoncfg->{DATADIR}/alignments/setToSet$queryDirLocal")
				{
					mkdir "$commoncfg->{DATADIR}/alignments/setToSet$queryDirLocal";
				}
			}
			foreach my $subjectSetId (keys %{$setToSetSwitched->{$querySetId}})
			{
				my $subjectDirLocal;
				for (my $position = 0; $position < length($subjectSetId); $position += 2)
				{
					$subjectDirLocal .= "/s". substr($subjectSetId,$position,2);
					unless (-e "$commoncfg->{DATADIR}/alignments/setToSet$queryDirLocal$subjectDirLocal")
					{
						mkdir "$commoncfg->{DATADIR}/alignments/setToSet$queryDirLocal$subjectDirLocal";
					}
				}
				open (LIST,">$commoncfg->{DATADIR}/alignments/setToSet$queryDirLocal$subjectDirLocal/$querySetId-$subjectSetId.list") or die "can't open file: $commoncfg->{DATADIR}/alignments/setToSet$queryDirLocal$subjectDirLocal/$querySetId-$subjectSetId.list";
				foreach (split ",", $setToSetSwitched->{$querySetId}->{$subjectSetId})
				{
					print LIST "$_\n";
				}
				close(LIST);
			}
		}

		#mark repeat regions to be updated!!!
# 		if ($markRepeatRegion)
# 		{
# 			my $converageCutoff = 0.95;
# 			my $todo = 1;
# 			do{
# 				$todo = 0;
# 				my $lastQuery = 0;
# 				my $lastQend = 0;
# 				my $lastAlignmentLength = 0;
# 				my $lastAlignmentId = 0;
# 				my $getAlignment = $dbh->prepare("SELECT * FROM alignment WHERE program = '$alignEngine\_1e-200\_$identityAlignment\_$minOverlapAlignment' AND align_length < 5000 AND hidden = 0 ORDER BY query,q_start,q_end");
# 				$getAlignment->execute();
# 				while (my @getAlignment = $getAlignment->fetchrow_array())
# 				{
# 					if($lastQuery == $getAlignment[2])
# 					{
# 						if($lastQend > $getAlignment[8])
# 						{
# 							my $converage = $lastQend - $getAlignment[8];
# 							if ( $converage/$lastAlignmentLength >= $converageCutoff)
# 							{
# 								my $hideLastAlignment=$dbh->do("UPDATE alignment SET hidden = 1 WHERE id = $lastAlignmentId");
# 								$todo = 1;
# 							}
# 							if ( $converage/$getAlignment[5] >= $converageCutoff)
# 							{
# 								my $hideAlignment=$dbh->do("UPDATE alignment SET hidden = 1 WHERE id = $getAlignment[0]");
# 								$todo = 1;
# 							}
# 						}
# 					}
# 					$lastQuery = $getAlignment[2];
# 					$lastQend = $getAlignment[9];
# 					$lastAlignmentLength = $getAlignment[5];
# 					$lastAlignmentId = $getAlignment[0];
# 				}
# 			} while($todo);
# 		}
		
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
		exit 0;
	}
	else
	{
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