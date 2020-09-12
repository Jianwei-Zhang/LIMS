#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use JSON::XS; #JSON::XS is recommended to be installed for handling JSON string of big size 
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
my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});

my $queryGenomeId = param ('queryGenomeId') || '';
my $subjectGenomeId = param ('subjectGenomeId') || '';
my $alignEngine = param ('alignEngine') || 'blastn';
my $alignmentFile = upload ('alignmentFile');
my $alignmentFilePath = param ('alignmentFilePath') || '';
my $alignmentInfile = "$commoncfg->{TMPDIR}/$$.tbl";

print header;

if($queryGenomeId && $subjectGenomeId && ($alignmentFile || $alignmentFilePath))
{
	my $pid = fork();
	if ($pid) {
		print <<END;
<script>
	parent.closeDialog();
	parent.errorPop("It's being loaded! This processing might take a while.");
</script>	
END
	}
	elsif($pid == 0){
		if($alignmentFilePath)
		{
			#$alignmentInfile = $alignmentFilePath;
			`cp $alignmentFilePath $alignmentInfile`;
		}
		else
		{
			open (FILE, ">$alignmentInfile");
			while (read ($alignmentFile, my $Buffer, 1024)) {
				print FILE $Buffer;
			}
			close FILE;
		}		
		`perl -p -i -e 's/\r/\n/g' $alignmentInfile`; #convert CR line terminators (MAC OS) into LF line terminators (Unix)
		close (STDOUT);
		#connect to the mysql server
		my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});

		my $setId;
		my $queryId;
		my $queryGenome=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$queryGenome->execute($queryGenomeId);
		my @queryGenome = $queryGenome->fetchrow_array();
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
					$queryId->{$getSequences[2]} = $getSequences[0];
					$setId->{$getSequences[0]} = $getSequences[4];

				}
			}
		}
		if($queryGenome[1] eq 'genome')
		{
			my $getSequences = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 99 AND x = ?");
			$getSequences->execute($queryGenomeId);
			while(my @getSequences = $getSequences->fetchrow_array())
			{
				$queryId->{$getSequences[2]} = $getSequences[0];
				$setId->{$getSequences[0]} = $getSequences[4];
			}
		}

		my $subjectId;
		if ($queryGenomeId eq $subjectGenomeId)
		{
			$subjectId = $queryId;
		}
		else
		{
			my $subjectGenome=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
			$subjectGenome->execute($subjectGenomeId);
			my @subjectGenome = $subjectGenome->fetchrow_array();
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
						$subjectId->{$getSequences[2]} = $getSequences[0];
						$setId->{$getSequences[0]} = $getSequences[4];
					}
				}
			}
			if($subjectGenome[1] eq 'genome')
			{
				my $getSequences = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 99 AND x = ?");
				$getSequences->execute($subjectGenomeId);
				while(my @getSequences = $getSequences->fetchrow_array())
				{
					$subjectId->{$getSequences[2]} = $getSequences[0];
					$setId->{$getSequences[0]} = $getSequences[4];
				}
			}
		}

		my $seqToSeq;
		my $seqToSet;
		my $setToSet;
		my $queryDir;
		my $subjectDir;
		my $seqToSetSwitched;
		my $setToSetSwitched;
		my $queryDirSwitched;
		my $subjectDirSwitched;

		open(ALN, "$alignmentInfile") or die "cannot open file $alignmentInfile";
		while(<ALN>)
		{
			chop;
			/^#/ and next;
			my @hit = split("\t",$_);
			$hit[12] = 0; #add a hidden column
			next if($hit[0] eq $hit[1]);
			if (exists $queryId->{$hit[0]} && exists $subjectId->{$hit[1]})
			{
				$hit[0] = $queryId->{$hit[0]};
				$hit[1] = $subjectId->{$hit[1]};

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
						print ALN "#$alignEngine\_1e-200\_$identityAlignment\_$minOverlapAlignment\n";
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
						print ALN "#$alignEngine\_manual_loading\n";
						print ALN "#query\tsubject\tperc_indentity\talign_length\tmismatches\tgaps\tq_start\tq_end\ts_start\ts_end\te_val\tbit_score\thidden\n";
						close(ALN);
					}
					$seqToSeq->{$hit[0]}->{$hit[1]} = 1;
				}

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

				if ($seqToSeq->{$hit[1]}->{$hit[0]} == 1) #check if this is the first hit
				{
					if (exists $seqToSetSwitched->{$hit[0]}->{$setId->{$hit[1]}})
					{
						$seqToSetSwitched->{$hit[0]}->{$setId->{$hit[1]}} .= ",alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$hit[0]-$hit[1].tbl";
			
					}
					else
					{
						$seqToSetSwitched->{$hit[0]}->{$setId->{$hit[1]}} = "alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$hit[0]-$hit[1].tbl";
					}
					if (exists $setToSetSwitched->{$setId->{$hit[0]}}->{$setId->{$hit[1]}})
					{
						$setToSetSwitched->{$setId->{$hit[0]}}->{$setId->{$hit[1]}} .= ",alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$hit[0]-$hit[1].tbl";
			
					}
					else
					{
						$setToSetSwitched->{$setId->{$hit[0]}}->{$setId->{$hit[1]}} = "alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$hit[0]-$hit[1].tbl";
					}
				}
				#write to alignment
				open (ALN,">>$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$hit[0]-$hit[1].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$hit[0]-$hit[1].tbl";
				print ALN join "\t", @hit;
				print ALN "\n";
				close(ALN);
				$seqToSeq->{$hit[0]}->{$hit[1]}++;
			}
			else
			{
				next;
			}
		}
		close(ALN);
		unlink ($alignmentInfile);

		foreach my $sequenceId (keys %$seqToSet)
		{
			my $queryDirLocal;
			for (my $position = 0; $position < length($sequenceId); $position += 2)
			{
				$queryDirLocal .= "/q". substr($sequenceId,$position,2);
				until (-e "$commoncfg->{DATADIR}/alignments/seqToSet$queryDirLocal")
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
					until (-e "$commoncfg->{DATADIR}/alignments/seqToSet$queryDirLocal$subjectDirLocal")
					{
						mkdir "$commoncfg->{DATADIR}/alignments/seqToSet$queryDirLocal$subjectDirLocal";
					}
				}
		
				until (-e "$commoncfg->{DATADIR}/alignments/seqToSet$queryDirLocal$subjectDirLocal/$sequenceId-$subjectSetId.list")
				{
					open (LIST,">$commoncfg->{DATADIR}/alignments/seqToSet$queryDirLocal$subjectDirLocal/$sequenceId-$subjectSetId.list") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSet$queryDirLocal$subjectDirLocal/$sequenceId-$subjectSetId.list";
					foreach (split ",", $seqToSet->{$sequenceId}->{$subjectSetId})
					{
						print LIST "$_\n";
					}
					close(LIST);
				}
			}
		}

		foreach my $querySetId (keys %$setToSet)
		{
			my $queryDirLocal;
			for (my $position = 0; $position < length($querySetId); $position += 2)
			{
				$queryDirLocal .= "/q". substr($querySetId,$position,2);
				until (-e "$commoncfg->{DATADIR}/alignments/setToSet$queryDirLocal")
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
					until (-e "$commoncfg->{DATADIR}/alignments/setToSet$queryDirLocal$subjectDirLocal")
					{
						mkdir "$commoncfg->{DATADIR}/alignments/setToSet$queryDirLocal$subjectDirLocal";
					}
				}
		
				until (-e "$commoncfg->{DATADIR}/alignments/setToSet$queryDirLocal$subjectDirLocal/$querySetId-$subjectSetId.list")
				{
					open (LIST,">$commoncfg->{DATADIR}/alignments/setToSet$queryDirLocal$subjectDirLocal/$querySetId-$subjectSetId.list") or die "can't open file: $commoncfg->{DATADIR}/alignments/setToSet$queryDirLocal$subjectDirLocal/$querySetId-$subjectSetId.list";
					foreach (split ",", $setToSet->{$querySetId}->{$subjectSetId})
					{
						print LIST "$_\n";
					}
					close(LIST);
				}
			}
		}

		foreach my $sequenceId (keys %$seqToSetSwitched)
		{
			my $queryDirLocal;
			for (my $position = 0; $position < length($sequenceId); $position += 2)
			{
				$queryDirLocal .= "/q". substr($sequenceId,$position,2);
				until (-e "$commoncfg->{DATADIR}/alignments/seqToSet$queryDirLocal")
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
					until (-e "$commoncfg->{DATADIR}/alignments/seqToSet$queryDirLocal$subjectDirLocal")
					{
						mkdir "$commoncfg->{DATADIR}/alignments/seqToSet$queryDirLocal$subjectDirLocal";
					}
				}
		
				until (-e "$commoncfg->{DATADIR}/alignments/seqToSet$queryDirLocal$subjectDirLocal/$sequenceId-$subjectSetId.list")
				{
					open (LIST,">$commoncfg->{DATADIR}/alignments/seqToSet$queryDirLocal$subjectDirLocal/$sequenceId-$subjectSetId.list") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSet$queryDirLocal$subjectDirLocal/$sequenceId-$subjectSetId.list";
					foreach (split ",", $seqToSetSwitched->{$sequenceId}->{$subjectSetId})
					{
						print LIST "$_\n";
					}
					close(LIST);
				}
			}
		}

		foreach my $querySetId (keys %$setToSetSwitched)
		{
			my $queryDirLocal;
			for (my $position = 0; $position < length($querySetId); $position += 2)
			{
				$queryDirLocal .= "/q". substr($querySetId,$position,2);
				until (-e "$commoncfg->{DATADIR}/alignments/setToSet$queryDirLocal")
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
					until (-e "$commoncfg->{DATADIR}/alignments/setToSet$queryDirLocal$subjectDirLocal")
					{
						mkdir "$commoncfg->{DATADIR}/alignments/setToSet$queryDirLocal$subjectDirLocal";
					}
				}
		
				until (-e "$commoncfg->{DATADIR}/alignments/setToSet$queryDirLocal$subjectDirLocal/$querySetId-$subjectSetId.list")
				{
					open (LIST,">$commoncfg->{DATADIR}/alignments/setToSet$queryDirLocal$subjectDirLocal/$querySetId-$subjectSetId.list") or die "can't open file: $commoncfg->{DATADIR}/alignments/setToSet$queryDirLocal$subjectDirLocal/$querySetId-$subjectSetId.list";
					foreach (split ",", $setToSetSwitched->{$querySetId}->{$subjectSetId})
					{
						print LIST "$_\n";
					}
					close(LIST);
				}
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
	parent.errorPop("Please provide required infomation or file!");
</script>	
END
}