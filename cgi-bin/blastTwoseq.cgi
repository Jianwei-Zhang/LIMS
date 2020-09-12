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
## sequences will be saved to $commoncfg->{DATADIR}/alignments
until (-e "$commoncfg->{DATADIR}/alignments")
{
	mkdir "$commoncfg->{DATADIR}/alignments";
}
until (-e "$commoncfg->{DATADIR}/alignments/seqToSeq")
{
	mkdir "$commoncfg->{DATADIR}/alignments/seqToSeq";
}
until (-e "$commoncfg->{DATADIR}/alignments/seqToSet")
{
	mkdir "$commoncfg->{DATADIR}/alignments/seqToSet";
}
until (-e "$commoncfg->{DATADIR}/alignments/setToSet")
{
	mkdir "$commoncfg->{DATADIR}/alignments/setToSet";
}

my $alignEngineList;
$alignEngineList->{'blastn'} = "blast+/bin/blastn";
$alignEngineList->{'BLAT'} = "blat";
my $windowmasker = 'blast+/bin/windowmasker';
my $makeblastdb = 'blast+/bin/makeblastdb';

my $identityBlast = param ('identityBlast') || $userConfig->getFieldValueWithUserIdAndFieldName($userId,"SEQTOSEQIDENTITY");
my $minOverlapBlast = param ('minOverlapBlast') || $userConfig->getFieldValueWithUserIdAndFieldName($userId,"SEQTOSEQMINOVERLAP");
my $alignEngine = param ('alignEngine') || 'blastn';
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
		my $getSequenceOne = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$getSequenceOne->execute($seqOne);
		my @getSequenceOne = $getSequenceOne->fetchrow_array();
		open (SEQA,">$commoncfg->{TMPDIR}/$getSequenceOne[0].$$.seq") or die "can't open file: $commoncfg->{TMPDIR}/$getSequenceOne[0].$$.seq";
		my $sequenceOneDetails = decode_json $getSequenceOne[8];
		my $sequenceOne = 'ERROR: NO SEQUENCE FOUND! PLEASE CONTACT YOUR ADMINISTRATOR.';
		my $in = Bio::SeqIO->new(-file => "$commoncfg->{DATADIR}/$sequenceOneDetails->{'sequence'}",
								-format => 'Fasta');
		while ( my $seq = $in->next_seq() )
		{
			$sequenceOne = $seq->seq;
		}
		print SEQA ">$getSequenceOne[0]\n$sequenceOne\n";
		close(SEQA);
		my $getSequenceTwo = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$getSequenceTwo->execute($seqTwo);
		my @getSequenceTwo =  $getSequenceTwo->fetchrow_array();
		open (SEQB,">$commoncfg->{TMPDIR}/$getSequenceTwo[0].$$.seq") or die "can't open file: $commoncfg->{TMPDIR}/$getSequenceTwo[0].$$.seq";
		my $sequenceTwoDetails = decode_json $getSequenceTwo[8];
		my $sequenceTwo = 'ERROR: NO SEQUENCE FOUND! PLEASE CONTACT YOUR ADMINISTRATOR.';
		my $in = Bio::SeqIO->new(-file => "$commoncfg->{DATADIR}/$sequenceTwoDetails->{'sequence'}",
								-format => 'Fasta');
		while ( my $seq = $in->next_seq() )
		{
			$sequenceTwo = $seq->seq;
		}
		print SEQB ">$getSequenceTwo[0]\n$sequenceTwo\n";
		close(SEQB);

		my $seqToSeq;
		my $queryDir;
		my $subjectDir;
		my $queryDirSwitched;
		my $subjectDirSwitched;

		unless (exists $seqToSeq->{$seqOne}->{$seqTwo}) # clean old data first
		{

			for (my $position = 0; $position < length($seqOne); $position += 2)
			{
				$queryDir .= "/q". substr($seqOne,$position,2);
				until (-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir")
				{
					mkdir "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir";
				}
			}

			for (my $position = 0; $position < length($seqTwo); $position += 2)
			{
				$subjectDir .= "/s". substr($seqTwo,$position,2);
				until (-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir")
				{
					mkdir "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir";
				}
			}

			unlink "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$seqOne-$seqTwo.tbl"; #delete old alignments
			until (-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$seqOne-$seqTwo.tbl")
			{
				open (ALN,">$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$seqOne-$seqTwo.tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$seqOne-$seqTwo.tbl";
				print ALN "#SEQtoSEQ_1e-200_$identityBlast\_$minOverlapBlast\n";
				print ALN "#query\tsubject\tperc_indentity\talign_length\tmismatches\tgaps\tq_start\tq_end\ts_start\ts_end\te_val\tbit_score\thidden\n";
				close(ALN);

			}

			for (my $position = 0; $position < length($seqTwo); $position += 2)
			{
				$queryDirSwitched .= "/q". substr($seqTwo,$position,2);
				until (-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched")
				{
					mkdir "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched";
				}
			}

			for (my $position = 0; $position < length($seqOne); $position += 2)
			{
				$subjectDirSwitched .= "/s". substr($seqOne,$position,2);
				until (-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched")
				{
					mkdir "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched";
				}
			}

			unlink "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$seqTwo-$seqOne.tbl"; #delete old alignments
			until (-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$seqTwo-$seqOne.tbl")
			{
				open (ALN,">$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$seqTwo-$seqOne.tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$seqTwo-$seqOne.tbl";
				print ALN "#SEQtoSEQ_1e-200_$identityBlast\_$minOverlapBlast\n";
				print ALN "#query\tsubject\tperc_indentity\talign_length\tmismatches\tgaps\tq_start\tq_end\ts_start\ts_end\te_val\tbit_score\thidden\n";
				close(ALN);
			}
			$seqToSeq->{$seqOne}->{$seqTwo} = 1;
		}

		my @alignments;
		my @alignmentsSwitched;
		my $goodOverlap = ($checkGood) ? 0 : 1;
		open (CMD,"$alignEngineList->{$alignEngine} -query $commoncfg->{TMPDIR}/$getSequenceOne[0].$$.seq -subject $commoncfg->{TMPDIR}/$getSequenceTwo[0].$$.seq -dust no -evalue 1e-200 -perc_identity $identityBlast -outfmt 6 |") or die "can't open CMD: $!";
		while(<CMD>)
		{
			chop;
			/^#/ and next;
			my @hit = split("\t",$_);
			$hit[12] = 0; #add a hidden column
			if($hit[3] >= $minOverlapBlast)
			{
				my $alignment = join "\t",@hit;
				push @alignments, $alignment;
				if($hit[6] == 1 || $hit[7] == $getSequenceOne[5])
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
					$exchange = $seqTwo;
					$seqTwo = $hit[0];
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
					$exchange = $seqTwo;
					$seqTwo = $hit[0];
					$hit[0] = $exchange;
				}
				if($hit[6] == 1 || $hit[7] == $getSequenceTwo[5])
				{
					$goodOverlap = 1;
				}
				$alignment = join "\t",@detailedHit;
				push @alignmentsSwitched, $alignment;							
			}									
		}
		close(CMD);
		if($goodOverlap)
		{
			foreach (@alignments)
			{
				my @detailedHit = split("\t",$_);
				#write to alignment
				open (ALN,">>$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$detailedHit[0]-$detailedHit[1].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$detailedHit[0]-$detailedHit[1].tbl";
				print ALN join "\t", @detailedHit;
				print ALN "\n";
				close(ALN);
			}
			foreach (@alignmentsSwitched)
			{
				my @detailedHit = split("\t",$_);
				#write to alignment
				open (ALN,">>$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$detailedHit[0]-$detailedHit[1].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$detailedHit[0]-$detailedHit[1].tbl";
				print ALN join "\t", @detailedHit;
				print ALN "\n";
				close(ALN);
			}
		}
		unlink("$commoncfg->{TMPDIR}/$getSequenceTwo[0].$$.seq");
		unlink("$commoncfg->{TMPDIR}/$getSequenceOne[0].$$.seq");
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
			my $seqToSeq;
			my $queryDir;
			my $subjectDir;
			my $queryDirSwitched;
			my $subjectDirSwitched;
			open (CMD,"$alignEngineList->{$alignEngine} -query $commoncfg->{TMPDIR}/$assembly[4].$querySeq.$$.seq -task $task -db $commoncfg->{TMPDIR}/$assembly[4].$$.seq -dust no -evalue 1e-200 -perc_identity $identityBlast -max_target_seqs 10 -num_threads 8 -outfmt 6 |") or die "can't open CMD: $!";
			while(<CMD>)
			{
				chop;
				/^#/ and next;
				my @hit = split("\t",$_);
				$hit[12] = 0; #add a hidden column
				next if($hit[0] eq $hit[1]);
				next if($hit[3] < $minOverlapBlast);

				unless (exists $seqToSeq->{$hit[0]}->{$hit[1]}) # clean old data first
				{
					for (my $position = 0; $position < length($hit[0]); $position += 2)
					{
						$queryDir .= "/q". substr($hit[0],$position,2);
						until (-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir")
						{
							mkdir "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir";
						}
					}

					for (my $position = 0; $position < length($hit[1]); $position += 2)
					{
						$subjectDir .= "/s". substr($hit[1],$position,2);
						until (-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir")
						{
							mkdir "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir";
						}
					}

					unlink "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$hit[0]-$hit[1].tbl"; #delete old alignments
					until (-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$hit[0]-$hit[1].tbl")
					{
						open (ALN,">$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$hit[0]-$hit[1].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$hit[0]-$hit[1].tbl";
						print ALN "#SEQtoSEQ_1e-200_$identityBlast\_$minOverlapBlast\n";
						print ALN "#query\tsubject\tperc_indentity\talign_length\tmismatches\tgaps\tq_start\tq_end\ts_start\ts_end\te_val\tbit_score\thidden\n";
						close(ALN);

					}

					for (my $position = 0; $position < length($hit[1]); $position += 2)
					{
						$queryDirSwitched .= "/q". substr($hit[1],$position,2);
						until (-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched")
						{
							mkdir "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched";
						}
					}

					for (my $position = 0; $position < length($hit[0]); $position += 2)
					{
						$subjectDirSwitched .= "/s". substr($hit[0],$position,2);
						until (-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched")
						{
							mkdir "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched";
						}
					}

					unlink "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$hit[1]-$hit[0].tbl"; #delete old alignments
					until (-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$hit[1]-$hit[0].tbl")
					{
						open (ALN,">$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$hit[1]-$hit[0].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$hit[1]-$hit[0].tbl";
						print ALN "#SEQtoSEQ_1e-200_$identityBlast\_$minOverlapBlast\n";
						print ALN "#query\tsubject\tperc_indentity\talign_length\tmismatches\tgaps\tq_start\tq_end\ts_start\ts_end\te_val\tbit_score\thidden\n";
						close(ALN);
					}
					$seqToSeq->{$hit[0]}->{$hit[1]} = 1;
				}

				if($speedyMode)
				{
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
					if ($seqToSeq->{$hit[1]}->{$hit[0]} == 1) #check if this is the first hit
					{
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
							foreach (@alignments)
							{
								my @detailedHit = split("\t",$_);
								#write to alignment
								open (ALN,">>$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$detailedHit[0]-$detailedHit[1].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$detailedHit[0]-$detailedHit[1].tbl";
								print ALN join "\t", @detailedHit;
								print ALN "\n";
								close(ALN);
							}
							foreach (@alignmentsSwitched)
							{
								my @detailedHit = split("\t",$_);
								#write to alignment
								open (ALN,">>$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$detailedHit[0]-$detailedHit[1].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$detailedHit[0]-$detailedHit[1].tbl";
								print ALN join "\t", @detailedHit;
								print ALN "\n";
								close(ALN);
							}
						}
					}
				}
				$seqToSeq->{$hit[1]}->{$hit[0]}++;
			}
			close(CMD);
			unlink("$commoncfg->{TMPDIR}/$assembly[4].$$.seq");
			unlink("$commoncfg->{TMPDIR}/$assembly[4].$querySeq.$$.seq");
			unlink("$commoncfg->{TMPDIR}/$assembly[4].$$.seq.nhr");
			unlink("$commoncfg->{TMPDIR}/$assembly[4].$$.seq.nin");
			unlink("$commoncfg->{TMPDIR}/$assembly[4].$$.seq.nsq");
			foreach my $queryId (keys %$seqToSeq)
			{
				unlink("$commoncfg->{TMPDIR}/$queryId.$$.seq");
				foreach my $subjectId (keys %{$seqToSeq->{$queryId}})
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