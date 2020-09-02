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
use userConfig;
use userCookie;

my $userCookie = new userCookie;
my $userId = (cookie('cid')) ? $userCookie->checkCookie(cookie('cid')) : 0;
exit if (!$userId);

my $user = new user;
my $userDetail = $user->getAllFieldsWithUserId($userId);
my $userName = $userDetail->{"userName"};

my $commoncfg = readConfig("main.conf");
my $userConfig = new userConfig;

my $alignEngineList;
$alignEngineList->{'blastn'} = "blast+/bin/blastn";
$alignEngineList->{'BLAT'} = "blat";
my $windowmasker = 'blast+/bin/windowmasker';
my $makeblastdb = 'blast+/bin/makeblastdb';

my $identityBlast = param ('identityBlast') || $userConfig->getFieldValueWithUserIdAndFieldName($userId,"SEQTOSEQIDENTITY");
my $minOverlapBlast = param ('minOverlapBlast') || $userConfig->getFieldValueWithUserIdAndFieldName($userId,"SEQTOSEQMINOVERLAP");
my $task = param ('megablast') || 'blastn';

my $speedyMode = param ('speedyMode') || '0';
my $checkGood = param ('checkGood') || '0';
my $seqOne = param ('seqOne') || '';
my $seqTwo = param ('seqTwo') || '';
my $assemblyId = param ('assemblyId') || '';
my $alignmentCheckFormUrl = "alignmentCheckForm.cgi";
if($assemblyId)
{
	$alignmentCheckFormUrl .= "?assemblyId=$assemblyId";
}
my $blastn = 'blast+/bin/blastn';


print header;

if($seqOne && $seqTwo)
{
	`rm $commoncfg->{TMPDIR}/*$seqOne*.aln.html`; #delete cached files
	`rm $commoncfg->{TMPDIR}/*$seqTwo*.aln.html`; #delete cached files
	if($seqOne == $seqTwo)
	{
		print <<END;
<script>
	parent.errorPop("Please provide two different sequences!");
</script>
END
		exit;
	}
	my $pid = fork();
	if ($pid) {
		print <<END;
<script>
	parent.closeDialog();
	parent.informationPop("It's running! This processing might take a while.");
	parent.openDialog('$alignmentCheckFormUrl&seqOne=$seqOne&seqTwo=$seqTwo');
</script>	
END
	}
	elsif($pid == 0){
		close (STDOUT);
		#connect to the mysql server
		my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});
		my $getSequenceA = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		my $deleteAlignmentA = $dbh->do("DELETE FROM alignment WHERE query = $seqOne AND subject = $seqTwo");
		my $deleteAlignmentB = $dbh->do("DELETE FROM alignment WHERE query = $seqTwo AND subject = $seqOne");
		$getSequenceA->execute($seqOne);
		my @getSequenceA = $getSequenceA->fetchrow_array();
		open (SEQA,">$commoncfg->{TMPDIR}/$getSequenceA[0].$$.seq") or die "can't open file: $commoncfg->{TMPDIR}/$getSequenceA[0].$$.seq";
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
		my $getSequenceB = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$getSequenceB->execute($seqTwo);
		my @getSequenceB =  $getSequenceB->fetchrow_array();
		open (SEQB,">$commoncfg->{TMPDIR}/$getSequenceB[0].$$.seq") or die "can't open file: $commoncfg->{TMPDIR}/$getSequenceB[0].$$.seq";
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

		my @alignments;
		my $goodOverlap = ($checkGood) ? 0 : 1;
		open (CMD,"$blastn -query $commoncfg->{TMPDIR}/$getSequenceA[0].$$.seq -subject $commoncfg->{TMPDIR}/$getSequenceB[0].$$.seq -dust no -evalue 1e-200 -perc_identity $identityBlast -outfmt 6 |") or die "can't open CMD: $!";
		while(<CMD>)
		{
			/^#/ and next;
			my @hit = split("\t",$_);
			if($hit[3] >= $minOverlapBlast)
			{
				push @alignments, $_;
				if($hit[6] == 1 || $hit[7] == $getSequenceA[5])
				{
					$goodOverlap = 1;
				}
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
				if($hit[6] == 1 || $hit[7] == $getSequenceB[5])
				{
					$goodOverlap = 1;
				}
				my $reverseBlast = join "\t",@hit;
				push @alignments, $reverseBlast;							
			}									
		}
		close(CMD);
		if($goodOverlap)
		{
			foreach (@alignments)
			{
				my @hit = split("\t",$_);
				#write to alignment
				my $insertAlignment=$dbh->prepare("INSERT INTO alignment VALUES ('', 'SEQtoSEQ_1e-200_$identityBlast\_$minOverlapBlast', ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 0)");
				$insertAlignment->execute(@hit);
			}
		}
		unlink("$commoncfg->{TMPDIR}/$getSequenceB[0].$$.seq");
		unlink("$commoncfg->{TMPDIR}/$getSequenceA[0].$$.seq");
	}
	else{
		die "couldn't fork: $!\n";
	} 
}
else
{
	if($seqOne || $seqTwo)
	{
		my $querySeq = ($seqOne) ? $seqOne : $seqTwo;
		`rm $commoncfg->{TMPDIR}/*$querySeq*.aln.html`; #delete cached files
		if (!$assemblyId)
		{
			#right now, this function can only be used for an assembly.
			print <<END;
<script>
	parent.errorPop("Please provide extra information!");
</script>	
END
			exit;
		}
		
		my $pid = fork();
		if ($pid) {
			print <<END;
	<script>
		parent.closeDialog();
		parent.informationPop("It's running! This processing might take a while.");
		parent.openDialog('$alignmentCheckFormUrl&seqOne=$querySeq');
	</script>	
END
		}
		elsif($pid == 0){
			close (STDOUT);
			#connect to the mysql server
			my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});
			my $assembly=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
			$assembly->execute($assemblyId);
			my @assembly = $assembly->fetchrow_array();
			my $target=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
			$target->execute($assembly[4]);
			my @target = $target->fetchrow_array();
			my $deleteAlignment = $dbh->do("DELETE FROM alignment WHERE (query = $querySeq OR subject = $querySeq) AND program LIKE 'SEQtoSEQ%'");

			my $assemblySequenceLength;
			open (SEQALL,">$commoncfg->{TMPDIR}/$assembly[4].$$.seq") or die "can't open file: $commoncfg->{TMPDIR}/$assembly[4].$$.seq";
			open (SEQNEW,">$commoncfg->{TMPDIR}/$assembly[4].$querySeq.$$.seq") or die "can't open file: $commoncfg->{TMPDIR}/$assembly[4].$querySeq.$$.seq";
			if($target[1] eq 'library')
			{
				my $getClones = $dbh->prepare("SELECT * FROM clones WHERE sequenced > 0 AND libraryId = ?");
				$getClones->execute($assembly[4]);
				while(my @getClones = $getClones->fetchrow_array())
				{
					my $getSequences = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o < 50 AND name LIKE ?");
					$getSequences->execute($getClones[1]);
					while(my @getSequences = $getSequences->fetchrow_array())
					{
						$assemblySequenceLength->{$getSequences[0]} = $getSequences[5];
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
						print SEQALL ">$getSequences[0]\n$sequence\n" if ($getSequences[0] ne $querySeq);
						print SEQNEW ">$getSequences[0]\n$sequence\n" if ($getSequences[0] eq $querySeq);
					}
				}
			}
			if($target[1] eq 'genome')
			{
				my $getSequences = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 99 AND x = ?");
				$getSequences->execute($assembly[4]);
				while(my @getSequences = $getSequences->fetchrow_array())
				{
					$assemblySequenceLength->{$getSequences[0]} = $getSequences[5];
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
					print SEQALL ">$getSequences[0]\n$sequence\n" if ($getSequences[0] ne $querySeq);
					print SEQNEW ">$getSequences[0]\n$sequence\n" if ($getSequences[0] eq $querySeq);
				}
			}
			close(SEQALL);
			close(SEQNEW);

			system( "$makeblastdb -in $commoncfg->{TMPDIR}/$assembly[4].$$.seq -dbtype nucl" );
			my $goodSequenceId;
			open (CMD,"$alignEngineList->{'blastn'} -query $commoncfg->{TMPDIR}/$assembly[4].$querySeq.$$.seq -task $task -db $commoncfg->{TMPDIR}/$assembly[4].$$.seq -dust no -evalue 1e-200 -perc_identity $identityBlast -max_target_seqs 10 -num_threads 8 -outfmt 6 |") or die "can't open CMD: $!";
			while(<CMD>)
			{
				/^#/ and next;
				my @hit = split("\t",$_);
				next if($hit[0] eq $hit[1]);
				next if($hit[3] < $minOverlapBlast);
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
					my $insertAlignmentA=$dbh->prepare("INSERT INTO alignment VALUES ('', 'SEQtoSEQ\_1e-200\_$identityBlast\_$minOverlapBlast', ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 0)");
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

					my $insertAlignmentB=$dbh->prepare("INSERT INTO alignment VALUES ('', 'SEQtoSEQ\_1e-200\_$identityBlast\_$minOverlapBlast', ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 0)");
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
						open (CMDA,"$alignEngineList->{'blastn'} -query $commoncfg->{TMPDIR}/$hit[0].$$.seq -subject $commoncfg->{TMPDIR}/$hit[1].$$.seq -dust no -evalue 1e-200 -perc_identity $identityBlast -outfmt 6 |") or die "can't open CMD: $!";
						while(<CMDA>)
						{
							/^#/ and next;
							my @detailedHit = split("\t",$_);
							if($detailedHit[3] >= $minOverlapBlast)
							{
								push @alignments, $_;
								if($detailedHit[6] == 1 || $detailedHit[7] == $assemblySequenceLength->{$detailedHit[0]})
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

								if($detailedHit[6] == 1 || $detailedHit[7] == $assemblySequenceLength->{$detailedHit[0]})
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
								my $insertAlignment=$dbh->prepare("INSERT INTO alignment VALUES ('', 'SEQtoSEQ\_1e-200\_$identityBlast\_$minOverlapBlast', ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 0)");
								$insertAlignment->execute(@detailedHit);
							}
						}
					}
				}
			}
			close(CMD);
			unlink("$commoncfg->{TMPDIR}/$assembly[4].$$.seq");
			unlink("$commoncfg->{TMPDIR}/$assembly[4].$querySeq.$$.seq");
			unlink("$commoncfg->{TMPDIR}/$assembly[4].$$.seq.nhr");
			unlink("$commoncfg->{TMPDIR}/$assembly[4].$$.seq.nin");
			unlink("$commoncfg->{TMPDIR}/$assembly[4].$$.seq.nsq");
			foreach my $queryId (keys %$goodSequenceId)
			{
				unlink("$commoncfg->{TMPDIR}/$queryId.$$.seq");
				foreach my $subjectId (keys %{$goodSequenceId->{$queryId}})
				{
					unlink("$commoncfg->{TMPDIR}/$subjectId.$$.seq");
				}
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
	parent.errorPop("Please give at least one seqeunce!");
</script>	
END
	}
}