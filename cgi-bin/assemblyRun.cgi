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

my $assemblyId = param ('assemblyId') || '';
my $replace = param ('replace') || '0';
my $withAllSequences = param ('withAllSequences') || '0';
my $fpcOrAgpId = param ('fpcOrAgpId') || '0';
my $refGenomeId = param ('refGenomeId') || '0';
my $assignChr = param ('assignChr') || '0';
my $reposition = param ('reposition') || '0';
my $alignmentBlockSize = param ('alignmentBlockSize') || $userConfig->getFieldValueWithUserIdAndFieldName($userId,"SEQTOGNMMINOVERLAP"); #1 kb
my $alignmentCoverage = param ('alignmentCoverage') || '25'; #25%
my $orientContigs = param ('orientContigs') || '0';
my $assemblySeqMinLength = param ('assemblySeqMinLength') || '0';
my $endToEnd = param ('endToEnd') || '0';
my $identityEndToEnd = param ('identityEndToEnd') || $userConfig->getFieldValueWithUserIdAndFieldName($userId,"ENDTOENDIDENTITY");
my $minOverlapEndToEnd = param ('minOverlapEndToEnd') || $userConfig->getFieldValueWithUserIdAndFieldName($userId,"ENDTOENDMINOVERLAP");
my $redundancyFilterSeq = param ('redundancyFilterSeq') || '0';
my $redundancySeq = 0;
my $redundancyFilterOverlap = param ('redundancyFilterOverlap') || '0';
my $orientSeqs = param ('orientSeqs') || '0';
my $renumber = param ('renumber') || '0';
my $emailNotification = param ('emailNotification') || '0';
my $gapLength = 100;
print header;

if($assemblyId)
{
	my $pid = fork();
	if ($pid) {
		print <<END;
<script>
	parent.closeDialog();
</script>	
END
	}
	elsif($pid == 0){
		#connect to the mysql server
		my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});

		my $checkRunningAssembly=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assembly' AND barcode < 0 ");
		$checkRunningAssembly->execute();
		if ($checkRunningAssembly->rows > 0)
		{
			print <<END;
<script>
	parent.errorPop("Assembly stopped due to another running assembly. Please try again later!");
</script>	
END
			exit;
		}

		my $assembly=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$assembly->execute($assemblyId);
		my @assembly = $assembly->fetchrow_array();
		if ($assembly[7] > 1 || $assembly[7] < 0) # exit if not for assembly
		{
			print <<END;
<script>
	parent.errorPop("This assembly is frozen or running. Nothing changed.");
</script>	
END
			exit;
		}
		else
		{
		print <<END;
<script>
	parent.refresh("menu");
	parent.errorPop("It's running! This processing might take a while.");
</script>	
END
		}
		close (STDOUT);

		my $target=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$target->execute($assembly[4]);
		my @target = $target->fetchrow_array();

		my $inAssemblySequenceId;
		my $sequenceLength;
		my $sequenceName;
		my $sequenceIdOfAssemblySeq;
		my $updateAssemblyToRunning=$dbh->do("UPDATE matrix SET barcode = '-1' WHERE id = $assemblyId");
		if($replace) #delete existing information
		{
			my $assemblySeqs = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblySeq' AND o = ?");
			$assemblySeqs->execute($assemblyId);
			while(my @assemblySeqs = $assemblySeqs->fetchrow_array())
			{
				$inAssemblySequenceId->{$assemblySeqs[5]} = 1;
			}
			my $deleteAssemblyCtg=$dbh->do("DELETE FROM matrix WHERE container LIKE 'assemblyCtg' AND o = $assemblyId");
			my $deleteAssemblySeq=$dbh->do("DELETE FROM matrix WHERE container LIKE 'assemblySeq' AND o = $assemblyId");

			#initialize assembly
			my $assemblyCtgNumber = 0;
			if($target[1] eq 'library')
			{
				my $getClones = $dbh->prepare("SELECT * FROM clones WHERE libraryId = ? AND sequenced > 0");
				$getClones->execute($assembly[4]);
				while(my @getClones= $getClones->fetchrow_array())
				{
					my @assemblySeqList;
					my $getSequences = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o < 50 AND y >= ? AND name LIKE ? ORDER BY y DESC");
					$getSequences->execute($assemblySeqMinLength,$getClones[1]);
					while(my @getSequences = $getSequences->fetchrow_array())
					{
						$sequenceLength->{$getSequences[0]} = $getSequences[5];
						$sequenceName->{$getSequences[0]} = $getSequences[2];
						my $insertAssemblySeq=$dbh->prepare("INSERT INTO matrix VALUES ('', 'assemblySeq', ?, ?, 0, ?, ?, 1, ?, ?, NOW())");
						$insertAssemblySeq->execute($getSequences[2],$assemblyId,$getSequences[0],$getSequences[5],"1,$getSequences[5]",$userName);
						my $assemblySeqId = $dbh->{mysql_insertid};
						$sequenceIdOfAssemblySeq->{$assemblySeqId} = $getSequences[0];
						if (@assemblySeqList)
						{
							push @assemblySeqList,"-($assemblySeqId)";
						}
						else
						{
							push @assemblySeqList,"($assemblySeqId)";
						}
					}
					my $assemblySeqList = join ",", @assemblySeqList;
					$assemblyCtgNumber++;
					my $insertAssemblyCtg=$dbh->prepare("INSERT INTO matrix VALUES ('', 'assemblyCtg', ?, ?, 0, 0, 0, 0, ?, ?, NOW())");
					$insertAssemblyCtg->execute($assemblyCtgNumber,$assemblyId,$assemblySeqList,$userName);
				}
			}
			if($target[1] eq 'genome')
			{
				my $getSequences = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 99 AND y >= ? AND x = ? ORDER BY id");
				$getSequences->execute($assemblySeqMinLength,$assembly[4]);
				while(my @getSequences = $getSequences->fetchrow_array())
				{
					$sequenceLength->{$getSequences[0]} = $getSequences[5];
					$sequenceName->{$getSequences[0]} = $getSequences[2];
					my $insertAssemblySeq=$dbh->prepare("INSERT INTO matrix VALUES ('', 'assemblySeq', ?, ?, 0, ?, ?, 1, ?, ?, NOW())");
					$insertAssemblySeq->execute($getSequences[2],$assemblyId,$getSequences[0],$getSequences[5],"1,$getSequences[5]",$userName);
					my $assemblySeqId = $dbh->{mysql_insertid};
					$sequenceIdOfAssemblySeq->{$assemblySeqId} = $getSequences[0];
					$assemblyCtgNumber++;
					my $insertAssemblyCtg=$dbh->prepare("INSERT INTO matrix VALUES ('', 'assemblyCtg', ?, ?, 0, 0, 0, 0, ?, ?, NOW())");
					$insertAssemblyCtg->execute($assemblyCtgNumber,$assemblyId,"($assemblySeqId)",$userName);
				}
			}
		}
		else
		{
			my $assemblyCtgNumber = 0;
			my $cloneInAssemblyCtg;
			my $assemblySeqs = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblySeq' AND o = ?");
			$assemblySeqs->execute($assemblyId);
			while(my @assemblySeqs = $assemblySeqs->fetchrow_array())
			{
				$inAssemblySequenceId->{$assemblySeqs[5]} = 1;
				my $assemblyCtgOfSeq = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblyCtg' AND MATCH (note) AGAINST (?)");
				$assemblyCtgOfSeq->execute($assemblySeqs[0]);
				my @assemblyCtgOfSeq = $assemblyCtgOfSeq->fetchrow_array();
				$sequenceLength->{$assemblySeqs[5]} = $assemblySeqs[6];
				$sequenceName->{$assemblySeqs[5]} = $assemblySeqs[2];
				$sequenceIdOfAssemblySeq->{$assemblySeqs[0]} = $assemblySeqs[5];
				$cloneInAssemblyCtg->{$assemblySeqs[2]} = $assemblyCtgOfSeq[0];
				$assemblyCtgNumber = $assemblyCtgOfSeq[2] if($assemblyCtgOfSeq[2] > $assemblyCtgNumber);
			}
			if($target[1] eq 'library')
			{
				my $getClones = $dbh->prepare("SELECT * FROM clones WHERE libraryId = ? AND sequenced > 0");
				$getClones->execute($assembly[4]);
				while(my @getClones= $getClones->fetchrow_array())
				{
					my $getSequences = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o < 50 AND y >= ? AND name LIKE ? ORDER BY y DESC");
					$getSequences->execute($assemblySeqMinLength,$getClones[1]);
					while(my @getSequences = $getSequences->fetchrow_array())
					{
						next if (exists $sequenceLength->{$getSequences[0]} || $withAllSequences < 1); #skip if seq has been assembled or not needed.
						$sequenceLength->{$getSequences[0]} = $getSequences[5];
						$sequenceName->{$getSequences[0]} = $getSequences[2];
						my $insertAssemblySeq=$dbh->prepare("INSERT INTO matrix VALUES ('', 'assemblySeq', ?, ?, 0, ?, ?, 1, ?, ?, NOW())");
						$insertAssemblySeq->execute($getSequences[2],$assemblyId,$getSequences[0],$getSequences[5],"1,$getSequences[5]",$userName);
						my $assemblySeqId = $dbh->{mysql_insertid};
						$sequenceIdOfAssemblySeq->{$assemblySeqId} = $getSequences[0];
						if(exists $cloneInAssemblyCtg->{$getSequences[2]})
						{
							my $assemblyCtg=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
							$assemblyCtg->execute($cloneInAssemblyCtg->{$getSequences[2]});
							my @assemblyCtg = $assemblyCtg->fetchrow_array();
							my @newAssemblyCtgSeqList;
							my $flag = 1;
							foreach (split ",", $assemblyCtg[8])
							{
								next unless ($_);
								push @newAssemblyCtgSeqList, $_;
								if($flag)
								{
									$_ =~ s/[^a-zA-Z0-9]//g;
									my $assemblySeqList=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
									$assemblySeqList->execute($_);
									my @assemblySeqList = $assemblySeqList->fetchrow_array();
									if($assemblySeqList[2] eq $getSequences[2])
									{
										push @newAssemblyCtgSeqList, "-($assemblySeqId)";
										$flag = 0;
									}
								}
							}								
							$assemblyCtg[8] = join ",",@newAssemblyCtgSeqList;
							my $updateAssemblyCtg=$dbh->prepare("UPDATE matrix SET note = ? WHERE id = ?");
							$updateAssemblyCtg->execute($assemblyCtg[8],$assemblyCtg[0]);
						}
						else
						{
							$assemblyCtgNumber++;
							my $insertAssemblyCtg=$dbh->prepare("INSERT INTO matrix VALUES ('', 'assemblyCtg', ?, ?, 0, 0, 0, 0, ?, ?, NOW())");
							$insertAssemblyCtg->execute($assemblyCtgNumber,$assemblyId,"($assemblySeqId)",$userName);
						}
					}
				}
			}
			if($target[1] eq 'genome')
			{
				my $getSequences = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 99 AND y >= ? AND x = ? ORDER BY id");
				$getSequences->execute($assemblySeqMinLength,$assembly[4]);
				while(my @getSequences = $getSequences->fetchrow_array())
				{
					next if (exists $sequenceLength->{$getSequences[0]} || $withAllSequences < 1); #skip if seq has been assembled or not needed.
					$assemblyCtgNumber++;
					$sequenceLength->{$getSequences[0]} = $getSequences[5];
					$sequenceName->{$getSequences[0]} = $getSequences[2];
					my $insertAssemblySeq=$dbh->prepare("INSERT INTO matrix VALUES ('', 'assemblySeq', ?, ?, 0, ?, ?, 1, ?, ?, NOW())");
					$insertAssemblySeq->execute($getSequences[2],$assemblyId,$getSequences[0],$getSequences[5],"1,$getSequences[5]",$userName);
					my $assemblySeqId = $dbh->{mysql_insertid};
					$sequenceIdOfAssemblySeq->{$assemblySeqId} = $getSequences[0];
					my $insertAssemblyCtg=$dbh->prepare("INSERT INTO matrix VALUES ('', 'assemblyCtg', ?, ?, 0, 0, 0, 0, ?, ?, NOW())");
					$insertAssemblyCtg->execute($assemblyCtgNumber,$assemblyId,"($assemblySeqId)",$userName);
				}
			}
		}

		if($fpcOrAgpId) #merge contigs based on fpc
		{
			my $updateAssemblyToRunningFpc=$dbh->do("UPDATE matrix SET barcode = '-3' WHERE id = $assemblyId");
			my $updateAssemblyFpcId=$dbh->do("UPDATE matrix SET z = $fpcOrAgpId WHERE id = $assemblyId");
			my $assemblySeqCtg;
			my $assemblySeqCtgId;
			my $assemblyCtgList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblyCtg' AND o = ?");
			$assemblyCtgList->execute($assemblyId);
			while (my @assemblyCtgList = $assemblyCtgList->fetchrow_array())
			{
				foreach (split ",", $assemblyCtgList[8])
				{
					next unless ($_);
					$_ =~ s/[^a-zA-Z0-9]//g;
					my $assemblySeq = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
					$assemblySeq->execute($_);
					my @assemblySeq = $assemblySeq->fetchrow_array();
					$assemblySeqCtg->{$assemblySeq[5]} = $assemblyCtgList[2];
					$assemblySeqCtgId->{$assemblySeq[5]} = $assemblyCtgList[0];
				}
			}

			my $fpcOrAgp=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
			$fpcOrAgp->execute($fpcOrAgpId);
			my @fpcOrAgp = $fpcOrAgp->fetchrow_array();
			
			if($fpcOrAgp[1] eq 'fpc')
			{
				my $fpcCtg;
				my $fpcCtgList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'fpcCtg' AND o = ?");
				$fpcCtgList->execute($fpcOrAgpId);
				while (my @fpcCtgList = $fpcCtgList->fetchrow_array())
				{
					$fpcCtg->{$fpcCtgList[2]} = $fpcCtgList[8];
				}
				my @fpcCtg = map  { $_->[0] }
						 sort { $a->[1] <=> $b->[1] }
						 map  { [$_, $_=~/(\d+)/] }
						 keys %$fpcCtg;
				for my $eachFpcCtg (@fpcCtg)
				{
					next if ($eachFpcCtg eq 'Ctg0');
					my $fpcSeq; #fpcSeq
					my $fpcCloneLeftEnd;
					my $fpcCloneRightEnd;
					my $fpcCloneMaxEnd = 0;
					my $fpcCloneList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'fpcClone' AND x > 0 AND o = ? AND MATCH (note) AGAINST (?)");
																							# fpcClone Sequenced(x > 0), no matter MTP marked(y > 0) or not.
					$fpcCloneList->execute($fpcOrAgpId,$eachFpcCtg);
					while (my @fpcCloneList = $fpcCloneList->fetchrow_array())
					{
						my $getSequences = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o < 50 AND name LIKE ?");
						$getSequences->execute($fpcCloneList[2]);
						while(my @getSequences = $getSequences->fetchrow_array())
						{
							$fpcSeq->{$getSequences[0]} = $fpcCloneList[8];
							$fpcCloneList[8] =~ /Ends Left (\d*)/;
							$fpcCloneLeftEnd->{$getSequences[0]} = $1;
							$fpcCloneList[8] =~ /Ends Right (\d*)/;
							$fpcCloneRightEnd->{$getSequences[0]} = $1;
							$fpcCloneMaxEnd = $1 if ($1 > $fpcCloneMaxEnd);
						}
					}
					my @fpcSeq = sort { $fpcCloneLeftEnd->{$a} <=> $fpcCloneLeftEnd->{$b} } keys %$fpcCloneLeftEnd;
					if(@fpcSeq > 1)
					{
						my $preSeq = shift @fpcSeq;
						for my $eachFpcSeq (@fpcSeq)
						{
							#1 read preSeq contig;
							#2 read currentSeq contig;
							#3 compare contigs, if the same, do nothing; if not, go to 4
							if($assemblySeqCtg->{$preSeq} ne $assemblySeqCtg->{$eachFpcSeq})
							{
								#4 check overlap, if not overlapped, do nothing; if yes, go to 5
								my $queryDir;
								for (my $position = 0; $position < length($preSeq); $position += 2)
								{
									$queryDir .= "/q". substr($preSeq,$position,2);
								}
								my $subjectDir;
								for (my $position = 0; $position < length($eachFpcSeq); $position += 2)
								{
									$subjectDir .= "/s". substr($eachFpcSeq,$position,2);
								}
								if(-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$getSequenceA[0]$preSeq-$eachFpcSeq.tbl")
								{
									#5 merge two contigs;
									my $assemblyPreCtg=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
									$assemblyPreCtg->execute($assemblySeqCtgId->{$preSeq});
									my @assemblyPreCtg = $assemblyPreCtg->fetchrow_array();
									my $assemblyCurrentCtg=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
									$assemblyCurrentCtg->execute($assemblySeqCtgId->{$eachFpcSeq});
									my @assemblyCurrentCtg = $assemblyCurrentCtg->fetchrow_array();
									my $mergedName = ($assemblyPreCtg[2] < $assemblyCurrentCtg[2]) ? $assemblyPreCtg[2] : $assemblyCurrentCtg[2];
									my $mergedAssemblySeq = "$assemblyPreCtg[8],$assemblyCurrentCtg[8]";
									foreach (split ",", $mergedAssemblySeq)
									{
										next unless ($_);
										$_ =~ s/[^a-zA-Z0-9]//g;
										my $assemblySeq = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
										$assemblySeq->execute($_);
										my @assemblySeq = $assemblySeq->fetchrow_array();
										$assemblySeqCtg->{$assemblySeq[5]} = $mergedName;
										$assemblySeqCtgId->{$assemblySeq[5]} = $assemblyPreCtg[0];
									}
									my $updatedAssemblyCtg = $dbh->prepare("UPDATE matrix SET name = ?, note = ? WHERE id = ?");
									$updatedAssemblyCtg->execute($mergedName,$mergedAssemblySeq,$assemblyPreCtg[0]);
									my $deleteAssemblyCtg=$dbh->do("DELETE FROM matrix WHERE id = $assemblyCurrentCtg[0]");
								}
							}
							$preSeq = $eachFpcSeq;
						}
					}
				}
			}
			elsif($fpcOrAgp[1] eq 'agp')
			{
				my $component;
				my $partNumber = -1;
				my $preAssemblyCtgId = 0;
				foreach (split "\n", $fpcOrAgp[8])
				{
					/^#/ and next;
					my @agpLine = split/\t/, $_;
				    if ($agpLine[4] ne 'N' && $agpLine[4] ne 'U' )
					{
						if (exists $component->{$agpLine[5]})
						{							
							#a component can't be used twice or more.
							$partNumber = $agpLine[3];
							next;
						}
						else
						{
							#set assemblySeq orientation based on $agpLine[8] and range based on "$agpLine[6]&$agpLine[7]"
							my $agpOrientation = ($agpLine[8] ne '-') ?  "1" : "-1";
							my $assemblySeqByName = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblySeq' AND o = ? AND name LIKE ?");
							$assemblySeqByName->execute($assemblyId,$agpLine[5]);
							my @assemblySeqByName = $assemblySeqByName->fetchrow_array();
							$agpLine[6] = $assemblySeqByName[6] if ($agpLine[6] > $assemblySeqByName[0]);
							$agpLine[7] = $assemblySeqByName[6] if ($agpLine[7] > $assemblySeqByName[0]);
							my $updateAssemblySeq=$dbh->prepare("UPDATE matrix SET barcode = ?, note = ? WHERE id = ?");
							$updateAssemblySeq->execute($agpOrientation,"$agpLine[6],$agpLine[7]",$assemblySeqByName[0]);
							if ($agpLine[3] - $partNumber == 1) #connected
							{
								#merge pre and current assemblyCtg
								unless($assemblySeqCtgId->{$assemblySeqByName[5]} == $preAssemblyCtgId)
								{
									my $assemblyPreCtg=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
									$assemblyPreCtg->execute($preAssemblyCtgId);
									my @assemblyPreCtg = $assemblyPreCtg->fetchrow_array();
									my $assemblyCurrentCtg=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
									$assemblyCurrentCtg->execute($assemblySeqCtgId->{$assemblySeqByName[5]});
									my @assemblyCurrentCtg = $assemblyCurrentCtg->fetchrow_array();
									my $mergedName = ($assemblyPreCtg[2] < $assemblyCurrentCtg[2]) ? $assemblyPreCtg[2] : $assemblyCurrentCtg[2];
									my $mergedAssemblySeq = "$assemblyPreCtg[8],$assemblyCurrentCtg[8]";
									foreach (split ",", $mergedAssemblySeq)
									{
										next unless ($_);
										$_ =~ s/[^a-zA-Z0-9]//g;
										my $assemblySeq = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
										$assemblySeq->execute($_);
										my @assemblySeq = $assemblySeq->fetchrow_array();
										$assemblySeqCtg->{$assemblySeq[5]} = $mergedName;
										$assemblySeqCtgId->{$assemblySeq[5]} = $assemblyPreCtg[0];
									}
									my $updatedAssemblyCtg = $dbh->prepare("UPDATE matrix SET name = ?, note = ? WHERE id = ?");
									$updatedAssemblyCtg->execute($mergedName,$mergedAssemblySeq,$assemblyPreCtg[0]);
									my $deleteAssemblyCtg=$dbh->do("DELETE FROM matrix WHERE id = $assemblyCurrentCtg[0]");
								}
							}
							else
							{
								unless($fpcOrAgp[5] != 1) # Chr-Seq
								{
									#set assemblyCtg position based on $agpLine[1] and chrNumber based on $agpLine[0]
									($agpLine[0] =~ /^ChrUN/) and next;
									($agpLine[0] =~ /^Ctg/) and next;
									my $chrNumber = 0;
									if ($agpLine[0] =~ /Chloroplast/)
									{
										$chrNumber = 98;
									}
									elsif ($agpLine[0] =~ /Mitochondrion/)
									{
										$chrNumber = 99;
									}
									elsif ($agpLine[0] =~ /Contamination/)
									{
										$chrNumber = 100;
									}
									else
									{
										$agpLine[0] =~ s/\D//g;
										$chrNumber = $agpLine[0];
									}
									my $updatedAssemblyCtg = $dbh->prepare("UPDATE matrix SET x = ?, z = ? WHERE id = ?");
									$updatedAssemblyCtg->execute($chrNumber,$agpLine[1],$assemblySeqCtgId->{$assemblySeqByName[5]});
								}
								# Ctg-Seq or unknown, no chrNumber and position assigned
								#do nothing
								$preAssemblyCtgId = $assemblySeqCtgId->{$assemblySeqByName[5]};
							}
							$partNumber = $agpLine[3];
							$component->{$agpLine[5]} = 1;
						}
					}
				}
			}
			my $updateAssemblyToRunningStatus=$dbh->do("UPDATE matrix SET barcode = '-1' WHERE id = $assemblyId");
		}

		if($redundancyFilterSeq)
		{
			my $updateAssemblyToFilteringSeq=$dbh->do("UPDATE matrix SET barcode = '-6' WHERE id = $assemblyId");
			my $assemblyMultiCtgList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblyCtg' AND o = ? AND note LIKE '%,%'");
			$assemblyMultiCtgList->execute($assemblyId);
			while (my @assemblyMultiCtgList = $assemblyMultiCtgList->fetchrow_array())
			{
				#buried seq
				my $todo = 1;
				do{
					$todo = 0;
					my @seqInCtg;
					foreach (split ",", $assemblyMultiCtgList[8])
					{
						next unless ($_);
						/^-/ and next;
						$_ =~ s/[^a-zA-Z0-9]//g;
						push @seqInCtg,$_;
					}
					my $totalSeqInCtg = @seqInCtg;
					for (my $i = 1; $i < $totalSeqInCtg; $i++)
					{
						my $startFound = 0;
						my $endFound = 0;
						my $coveredLength = 0;
						if($sequenceLength->{$sequenceIdOfAssemblySeq->{$seqInCtg[$i-1]}} <= $sequenceLength->{$sequenceIdOfAssemblySeq->{$seqInCtg[$i]}})
						{
							my $queryDir;
							for (my $position = 0; $position < length($sequenceIdOfAssemblySeq->{$seqInCtg[$i-1]}); $position += 2)
							{
								$queryDir .= "/q". substr($sequenceIdOfAssemblySeq->{$seqInCtg[$i-1]},$position,2);
							}
							my $subjectDir;
							for (my $position = 0; $position < length($sequenceIdOfAssemblySeq->{$seqInCtg[$i]}); $position += 2)
							{
								$subjectDir .= "/s". substr($sequenceIdOfAssemblySeq->{$seqInCtg[$i]},$position,2);
							}
							if(-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$sequenceIdOfAssemblySeq->{$seqInCtg[$i-1]}-$sequenceIdOfAssemblySeq->{$seqInCtg[$i]}.tbl")
							{
								open (TBL,"$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$sequenceIdOfAssemblySeq->{$seqInCtg[$i-1]}-$sequenceIdOfAssemblySeq->{$seqInCtg[$i]}.tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$sequenceIdOfAssemblySeq->{$seqInCtg[$i-1]}-$sequenceIdOfAssemblySeq->{$seqInCtg[$i]}.tbl";
								while(<TBL>)
								{
									chop;
									/^#/ and next;
									my @getAlignment = split("\t",$_);
									$startFound = 1 if ($getAlignment[6] == 1);
									$endFound = 1 if ($getAlignment[7] == $sequenceLength->{$sequenceIdOfAssemblySeq->{$seqInCtg[$i-1]}});
									$coveredLength += $getAlignment[7] - $getAlignment[6] + 1;
								}
								close(TBL);
							}

							if ($startFound > 0 && $endFound > 0 && $coveredLength / $sequenceLength->{$sequenceIdOfAssemblySeq->{$seqInCtg[$i-1]}} >= 0.99)
							{
								$assemblyMultiCtgList[8] =~ s/\($seqInCtg[$i-1]\)/-($seqInCtg[$i-1])/g;
								my $updateAssemblyCtgNote=$dbh->do("UPDATE matrix SET note = '$assemblyMultiCtgList[8]' WHERE id = $assemblyMultiCtgList[0]");
								$todo = 1;
								$redundancySeq++;
							}
						}
						else
						{
							my $queryDir;
							for (my $position = 0; $position < length($sequenceIdOfAssemblySeq->{$seqInCtg[$i]}); $position += 2)
							{
								$queryDir .= "/q". substr($sequenceIdOfAssemblySeq->{$seqInCtg[$i]},$position,2);
							}
							my $subjectDir;
							for (my $position = 0; $position < length($sequenceIdOfAssemblySeq->{$seqInCtg[$i-1]}); $position += 2)
							{
								$subjectDir .= "/s". substr($sequenceIdOfAssemblySeq->{$seqInCtg[$i-1]},$position,2);
							}
							if(-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$sequenceIdOfAssemblySeq->{$seqInCtg[$i]}-$sequenceIdOfAssemblySeq->{$seqInCtg[$i-1]}.tbl")
							{
								open (TBL,"$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$sequenceIdOfAssemblySeq->{$seqInCtg[$i]}-$sequenceIdOfAssemblySeq->{$seqInCtg[$i-1]}.tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$sequenceIdOfAssemblySeq->{$seqInCtg[$i]}-$sequenceIdOfAssemblySeq->{$seqInCtg[$i-1]}.tbl";
								while(<TBL>)
								{
									chop;
									/^#/ and next;
									my @getAlignment = split("\t",$_);
									$startFound = 1 if ($getAlignment[6] == 1);
									$endFound = 1 if ($getAlignment[7] == $sequenceLength->{$sequenceIdOfAssemblySeq->{$seqInCtg[$i]}});
									$coveredLength += $getAlignment[7] - $getAlignment[6] + 1;
								}
								close(TBL);
							}
							if ($startFound > 0 && $endFound > 0 && $coveredLength / $sequenceLength->{$sequenceIdOfAssemblySeq->{$seqInCtg[$i]}} >= 0.99)
							{
								$assemblyMultiCtgList[8] =~ s/\($seqInCtg[$i]\)/-($seqInCtg[$i])/g;
								my $updateAssemblyCtgNote=$dbh->do("UPDATE matrix SET note = '$assemblyMultiCtgList[8]' WHERE id = $assemblyMultiCtgList[0]");
								$todo = 1;
								$redundancySeq++;
							}
						}
					}
				} while ($todo);

				#redundant seq (covered by pre and next seqs)
				$todo = 1;
				do{
					$todo = 0;
					my @seqInCtg;
					foreach (split ",", $assemblyMultiCtgList[8])
					{
						next unless ($_);
						/^-/ and next;
						$_ =~ s/[^a-zA-Z0-9]//g;
						push @seqInCtg,$_;
					}
					my $totalSeqInCtg = @seqInCtg;
					for (my $i = 1; $i < $totalSeqInCtg - 1; $i++)
					{
						my $queryDir;
						for (my $position = 0; $position < length($sequenceIdOfAssemblySeq->{$seqInCtg[$i-1]}); $position += 2)
						{
							$queryDir .= "/q". substr($sequenceIdOfAssemblySeq->{$seqInCtg[$i-1]},$position,2);
						}
						my $subjectDir;
						for (my $position = 0; $position < length($sequenceIdOfAssemblySeq->{$seqInCtg[$i+1]}); $position += 2)
						{
							$subjectDir .= "/s". substr($sequenceIdOfAssemblySeq->{$seqInCtg[$i+1]},$position,2);
						}
						if (-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$sequenceIdOfAssemblySeq->{$seqInCtg[$i-1]}-$sequenceIdOfAssemblySeq->{$seqInCtg[$i+1]}.tbl")
						{
							my $goodEndFound = 0;
							open (TBL,"$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$sequenceIdOfAssemblySeq->{$seqInCtg[$i-1]}-$sequenceIdOfAssemblySeq->{$seqInCtg[$i+1]}.tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$sequenceIdOfAssemblySeq->{$seqInCtg[$i-1]}-$sequenceIdOfAssemblySeq->{$seqInCtg[$i+1]}.tbl";
							while(<TBL>)
							{
								chop;
								/^#/ and next;
								my @getAlignment = split("\t",$_);
								if($getAlignment[8] < $getAlignment[9])
								{
									$goodEndFound = 1 if ($getAlignment[6] == 1 || $getAlignment[7] == $sequenceLength->{$sequenceIdOfAssemblySeq->{$seqInCtg[$i-1]}} || $getAlignment[8] == 1 || $getAlignment[9] == $sequenceLength->{$sequenceIdOfAssemblySeq->{$seqInCtg[$i+1]}});
								}
								else
								{
									$goodEndFound = 1 if ($getAlignment[6] == 1 || $getAlignment[7] == $sequenceLength->{$sequenceIdOfAssemblySeq->{$seqInCtg[$i-1]}} || $getAlignment[8] == $sequenceLength->{$sequenceIdOfAssemblySeq->{$seqInCtg[$i+1]}} || $getAlignment[9] == 1);
								}
							}
							close(TBL);
							if($goodEndFound)
							{
								my $startFound = 0;
								my $endFound = 0;
								my $coveredLength = 0;

								my $queryDirPre;
								for (my $position = 0; $position < length($sequenceIdOfAssemblySeq->{$seqInCtg[$i]}); $position += 2)
								{
									$queryDirPre .= "/q". substr($sequenceIdOfAssemblySeq->{$seqInCtg[$i]},$position,2);
								}
								my $subjectDirPre;
								for (my $position = 0; $position < length($sequenceIdOfAssemblySeq->{$seqInCtg[$i-1]}); $position += 2)
								{
									$subjectDirPre .= "/s". substr($sequenceIdOfAssemblySeq->{$seqInCtg[$i-1]},$position,2);
								}
								if(-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirPre$subjectDirPre/$sequenceIdOfAssemblySeq->{$seqInCtg[$i]}-$sequenceIdOfAssemblySeq->{$seqInCtg[$i-1]}.tbl")
								{
									open (TBL,"$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirPre$subjectDirPre/$sequenceIdOfAssemblySeq->{$seqInCtg[$i]}-$sequenceIdOfAssemblySeq->{$seqInCtg[$i-1]}.tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDirPre$subjectDirPre/$sequenceIdOfAssemblySeq->{$seqInCtg[$i]}-$sequenceIdOfAssemblySeq->{$seqInCtg[$i-1]}.tbl";
									while(<TBL>)
									{
										chop;
										/^#/ and next;
										my @getAlignmentPre = split("\t",$_);
										$startFound = 1 if ($getAlignmentPre[6] == 1);
										$endFound = 1 if ($getAlignmentPre[7] == $sequenceLength->{$sequenceIdOfAssemblySeq->{$seqInCtg[$i]}});
										$coveredLength += $getAlignmentPre[7] - $getAlignmentPre[6] + 1;
									}
									close(TBL);
								}

								my $queryDirNext;
								for (my $position = 0; $position < length($sequenceIdOfAssemblySeq->{$seqInCtg[$i]}); $position += 2)
								{
									$queryDirNext .= "/q". substr($sequenceIdOfAssemblySeq->{$seqInCtg[$i]},$position,2);
								}
								my $subjectDirNext;
								for (my $position = 0; $position < length($sequenceIdOfAssemblySeq->{$seqInCtg[$i+1]}); $position += 2)
								{
									$subjectDirNext .= "/s". substr($sequenceIdOfAssemblySeq->{$seqInCtg[$i+1]},$position,2);
								}
								if(-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirNext$subjectDirNext/$sequenceIdOfAssemblySeq->{$seqInCtg[$i]}-$sequenceIdOfAssemblySeq->{$seqInCtg[$i+1]}.tbl")
								{
									open (TBL,"$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirNext$subjectDirNext/$sequenceIdOfAssemblySeq->{$seqInCtg[$i]}-$sequenceIdOfAssemblySeq->{$seqInCtg[$i+1]}.tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDirNext$subjectDirNext/$sequenceIdOfAssemblySeq->{$seqInCtg[$i]}-$sequenceIdOfAssemblySeq->{$seqInCtg[$i+1]}.tbl";
									while(<TBL>)
									{
										chop;
										/^#/ and next;
										my @getAlignmentNext = split("\t",$_);
										$startFound = 1 if ($getAlignmentNext[6] == 1);
										$endFound = 1 if ($getAlignmentNext[7] == $sequenceLength->{$sequenceIdOfAssemblySeq->{$seqInCtg[$i]}});
										$coveredLength += $getAlignmentNext[7] - $getAlignmentNext[6] + 1;
									}
									close(TBL);
								}
								if ($startFound > 0 && $endFound > 0 && $coveredLength / $sequenceLength->{$sequenceIdOfAssemblySeq->{$seqInCtg[$i]}} > 1)
								{
									$assemblyMultiCtgList[8] =~ s/\($seqInCtg[$i]\)/-($seqInCtg[$i])/g;
									my $updateAssemblyCtgNote=$dbh->do("UPDATE matrix SET note = '$assemblyMultiCtgList[8]' WHERE id = $assemblyMultiCtgList[0]");
									$todo = 1;
									$redundancySeq++;
								}
							}
						}
						$i = $totalSeqInCtg if ($todo);
					}
				} while ($todo);
			}
			my $updateAssemblyToRunningStatus=$dbh->do("UPDATE matrix SET barcode = '-1' WHERE id = $assemblyId");
		}

		if($orientSeqs)
		{
			my $updateAssemblyToAutoOrientation=$dbh->do("UPDATE matrix SET barcode = '-7' WHERE id = $assemblyId");
			my $assemblyMultiCtgList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblyCtg' AND o = ? AND note LIKE '%,%'");
			$assemblyMultiCtgList->execute($assemblyId);
			while (my @assemblyMultiCtgList = $assemblyMultiCtgList->fetchrow_array())
			{
				my @seqInCtg;
				foreach (split ",", $assemblyMultiCtgList[8])
				{
					next unless ($_);
					/^-/ and next;
					$_ =~ s/[^a-zA-Z0-9]//g;						
					push @seqInCtg,$_;
				}
				next if (@seqInCtg < 2);
				#determine orientation based on ends
				my $orientation;
				my $preAssemblySeqId = shift @seqInCtg;
				$orientation->{$preAssemblySeqId} = 1;
				my $firstAssemblySeqFlag = 1;
				for my $currentAssemblySeqId (@seqInCtg)
				{
					#determine orientation based on ends
					my $preAssemblySeq=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
					$preAssemblySeq->execute($preAssemblySeqId);
					my @preAssemblySeq = $preAssemblySeq->fetchrow_array();
					my $currentAssemblySeq=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
					$currentAssemblySeq->execute($currentAssemblySeqId);
					my @currentAssemblySeq = $currentAssemblySeq->fetchrow_array();

					my $twoSeqDirection;
					my $queryDirTwo;
					for (my $position = 0; $position < length($preAssemblySeq[5]); $position += 2)
					{
						$queryDirTwo .= "/q". substr($preAssemblySeq[5],$position,2);
					}
					my $subjectDirTwo;
					for (my $position = 0; $position < length($currentAssemblySeq[5]); $position += 2)
					{
						$subjectDirTwo .= "/s". substr($currentAssemblySeq[5],$position,2);
					}
					if(-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirTwo$subjectDirTwo/$preAssemblySeq[5]-$currentAssemblySeq[5].tbl")
					{
						open (TBL,"$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirTwo$subjectDirTwo/$preAssemblySeq[5]-$currentAssemblySeq[5].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDirTwo$subjectDirTwo/$preAssemblySeq[5]-$currentAssemblySeq[5].tbl";
						while(<TBL>)
						{
							chop;
							/^#/ and next;
							my @getAlignmentTwo = split("\t",$_);

							unless(exists $twoSeqDirection->{$preAssemblySeq[5]}->{$currentAssemblySeq[5]})
							{
								$twoSeqDirection->{$preAssemblySeq[5]}->{$currentAssemblySeq[5]} = ($getAlignmentTwo[8] < $getAlignmentTwo[9]) ? 1 : -1;
							}
							if($getAlignmentTwo[6] == 1)
							{
								$orientation->{$preAssemblySeqId} = -1;
								last;
							}
							if($getAlignmentTwo[7] == $preAssemblySeq[6])
							{
								$orientation->{$preAssemblySeqId} = 1;
								last;
							}
						}
						close(TBL);
					}
					if ($firstAssemblySeqFlag)
					{
						my $updatePreAssemblySeq=$dbh->do("UPDATE matrix SET barcode = $orientation->{$preAssemblySeqId} WHERE id = $preAssemblySeqId");
						$firstAssemblySeqFlag = 0;
					}
					$twoSeqDirection->{$preAssemblySeq[5]}->{$currentAssemblySeq[5]} = 1 unless (exists $twoSeqDirection->{$preAssemblySeq[5]}->{$currentAssemblySeq[5]});
					$orientation->{$currentAssemblySeqId} = $orientation->{$preAssemblySeqId} * $twoSeqDirection->{$preAssemblySeq[5]}->{$currentAssemblySeq[5]};
					my $updateAssemblySeq=$dbh->do("UPDATE matrix SET barcode = $orientation->{$currentAssemblySeqId} WHERE id = $currentAssemblySeqId");
					$preAssemblySeqId = $currentAssemblySeqId;
				}
			}			
			my $updateAssemblyToRunningStatus=$dbh->do("UPDATE matrix SET barcode = '-1' WHERE id = $assemblyId");
		}

		# end-to-end merge
		if($endToEnd)
		{
			my $updateAssemblyToRunningEndToEnd=$dbh->do("UPDATE matrix SET barcode = '-9' WHERE id = $assemblyId");
			if($refGenomeId)
			{
				my $ctgSeq;
				my $assemblyCtgIdFromCtg;
				my $assemblyChr;
				my $preCtg = '';
				my $ctgToChr=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblyCtg' AND o = ?  ORDER BY x, y");
				$ctgToChr->execute($assemblyId);
				while (my @ctgToChr = $ctgToChr->fetchrow_array())
				{
					$ctgSeq->{$ctgToChr[2]}=$ctgToChr[8];
					$assemblyChr->{$ctgToChr[2]}=$ctgToChr[4];
					$assemblyCtgIdFromCtg->{$ctgToChr[2]}=$ctgToChr[0];
					my $setPreCtg = 1;
					if ($preCtg && $assemblyChr->{$preCtg}==$assemblyChr->{$ctgToChr[2]})
					{
						my @preSeqList=split ",", $ctgSeq->{$preCtg};
						my $preCtgSeqId = '';
						for (@preSeqList)
						{
							next unless ($_);
							next if ($_ =~ /^-/);
							$_ =~ s/[^a-zA-Z0-9]//g;
							$preCtgSeqId = $_;
						}
						my $preCtgSeq=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
						$preCtgSeq->execute($preCtgSeqId);
						my @preCtgSeq = $preCtgSeq->fetchrow_array();
					  
						my @currentSeqList=split ",", $ctgSeq->{$ctgToChr[2]};
						my $currentCtgSeqId = '';
						for (@currentSeqList)
						{
							next unless ($_);
							next if ($_ =~ /^-/);
							$_ =~ s/[^a-zA-Z0-9]//g;
						 	$currentCtgSeqId = $_;
						 	last if ($currentCtgSeqId);
						}
						my $currentCtgSeq=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
						$currentCtgSeq->execute($currentCtgSeqId);
						my @currentCtgSeq = $currentCtgSeq->fetchrow_array();
						my $queryDir;
						for (my $position = 0; $position < length($preCtgSeq[5]); $position += 2)
						{
							$queryDir .= "/q". substr($preCtgSeq[5],$position,2);
						}
						my $subjectDir;
						for (my $position = 0; $position < length($currentCtgSeq[5]); $position += 2)
						{
							$subjectDir .= "/s". substr($currentCtgSeq[5],$position,2);
						}
						if (-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$preCtgSeq[5]-$currentCtgSeq[5].tbl")
						{
							my $goodEndFoundA = 0;
							my $goodEndFoundB = 0;
							open (TBL,"$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$preCtgSeq[5]-$currentCtgSeq[5].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$preCtgSeq[5]-$currentCtgSeq[5].tbl";
							while(<TBL>)
							{
								chop;
								/^#/ and next;
								my @getAlignment = split("\t",$_);
								$goodEndFoundA = 1 if ($getAlignment[6] == 1 || $getAlignment[7] == $preCtgSeq[6]);
								if($getAlignment[8] < $getAlignment[9])
								{
									$goodEndFoundB = 1 if ($getAlignment[8] == 1 || $getAlignment[9] == $currentCtgSeq[6]);
								}
								else
								{	
									 $goodEndFoundB = 1 if ($getAlignment[8] == $currentCtgSeq[6] || $getAlignment[9] == 1);
								}
							}
							close(TBL);
							if($goodEndFoundA && $goodEndFoundB)
							{
								if($preCtg < $ctgToChr[2])
								{
									$ctgSeq->{$preCtg} = "$ctgSeq->{$preCtg},$ctgSeq->{$ctgToChr[2]}";
									my $updatedAssemblyCtg = $dbh->prepare("UPDATE matrix SET name = ?, note = ? WHERE id = ?");
									$updatedAssemblyCtg->execute($preCtg,$ctgSeq->{$preCtg},$assemblyCtgIdFromCtg->{$preCtg});
									my $deleteAssemblyCtg=$dbh->do("DELETE FROM matrix WHERE id = $assemblyCtgIdFromCtg->{$ctgToChr[2]}");
									delete $assemblyCtgIdFromCtg->{$ctgToChr[2]};       
									delete $ctgSeq->{$ctgToChr[2]};
									$setPreCtg = 0;
								 }
							   else
							   {
									$ctgSeq->{$ctgToChr[2]} = "$ctgSeq->{$preCtg},$ctgSeq->{$ctgToChr[2]}";
									my $updatedAssemblyCtg = $dbh->prepare("UPDATE matrix SET name = ?, note = ? WHERE id = ?");
									$updatedAssemblyCtg->execute($ctgToChr[2],$ctgSeq->{$ctgToChr[2]},$assemblyCtgIdFromCtg->{$ctgToChr[2]});
									my $deleteAssemblyCtg=$dbh->do("DELETE FROM matrix WHERE id = $assemblyCtgIdFromCtg->{$preCtg}");
									delete $assemblyCtgIdFromCtg->{$preCtg};        
									delete $ctgSeq->{$preCtg};
							   }
							}
						}
					}
					$preCtg=$ctgToChr[2] if ($setPreCtg);
				}
			}
			else
			{
				#merge singletons first
				my $singletonCtg;
				my $singletonCtgId;
				my $assemblySingletonCtgList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblyCtg' AND o = ? AND note NOT LIKE '%,%'");
				$assemblySingletonCtgList->execute($assemblyId);
				while (my @assemblySingletonCtgList = $assemblySingletonCtgList->fetchrow_array())
				{
					next if ($assemblySingletonCtgList[8] =~ /^-/);
					$singletonCtg->{$assemblySingletonCtgList[2]} = $assemblySingletonCtgList[8];
					$singletonCtgId->{$assemblySingletonCtgList[2]} = $assemblySingletonCtgList[0];
				}
				my $merged;
				for my $eachSingletonCtgA (sort {$a <=> $b} keys %$singletonCtg)
				{
					next if (exists $merged->{$eachSingletonCtgA});
					$merged->{$eachSingletonCtgA} = 1;
					for my $eachSingletonCtgB (sort {$a <=> $b} keys %$singletonCtg)
					{
						next if (exists $merged->{$eachSingletonCtgB});
						my $singletonCtgSeqIdA = $singletonCtg->{$eachSingletonCtgA};
						$singletonCtgSeqIdA =~ s/[^a-zA-Z0-9]//g;
						my $eachSingletonCtgSeqA=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
						$eachSingletonCtgSeqA->execute($singletonCtgSeqIdA);
						my @eachSingletonCtgSeqA = $eachSingletonCtgSeqA->fetchrow_array();
						
						my $singletonCtgSeqIdB = $singletonCtg->{$eachSingletonCtgB};
						$singletonCtgSeqIdB =~ s/[^a-zA-Z0-9]//g;
						my $eachSingletonCtgSeqB=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
						$eachSingletonCtgSeqB->execute($singletonCtgSeqIdB);
						my @eachSingletonCtgSeqB = $eachSingletonCtgSeqB->fetchrow_array();

						my $queryDir;
						for (my $position = 0; $position < length($eachSingletonCtgSeqA[5]); $position += 2)
						{
							$queryDir .= "/q". substr($eachSingletonCtgSeqA[5],$position,2);
						}
						my $subjectDir;
						for (my $position = 0; $position < length($eachSingletonCtgSeqB[5]); $position += 2)
						{
							$subjectDir .= "/s". substr($eachSingletonCtgSeqB[5],$position,2);
						}
						if (-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$eachSingletonCtgSeqA[5]-$eachSingletonCtgSeqB[5].tbl")
						{
							my $goodEndFoundA = 0;
							my $goodEndFoundB = 0;
							open (TBL,"$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$eachSingletonCtgSeqA[5]-$eachSingletonCtgSeqB[5].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$preCtgSeq[5]-$eachSingletonCtgSeqB[5].tbl";
							while(<TBL>)
							{
								chop;
								/^#/ and next;
								my @getAlignment = split("\t",$_);
								$goodEndFoundA = 1 if ($getAlignment[6] == 1 || $getAlignment[7] == $eachSingletonCtgSeqA[6]);
								if($getAlignment[8] < $getAlignment[9])
								{
									$goodEndFoundB = 1 if ($getAlignment[8] == 1 || $getAlignment[9] == $eachSingletonCtgSeqB[6]);
								}
								else
								{
									$goodEndFoundB = 1 if ($getAlignment[8] == $eachSingletonCtgSeqB[6] || $getAlignment[9] == 1);
								}
							}
							close(TBL);
							if($goodEndFoundA && $goodEndFoundB)
							{
								my $mergedName = ($eachSingletonCtgA < $eachSingletonCtgB) ? $eachSingletonCtgA : $eachSingletonCtgB;
								my $updatedAssemblyCtg = $dbh->prepare("UPDATE matrix SET name = ?, note = ? WHERE id = ?");
								$updatedAssemblyCtg->execute($mergedName,"$singletonCtg->{$eachSingletonCtgA},$singletonCtg->{$eachSingletonCtgB}",$singletonCtgId->{$eachSingletonCtgA});
								my $deleteAssemblyCtg=$dbh->do("DELETE FROM matrix WHERE id = $singletonCtgId->{$eachSingletonCtgB}");
								$merged->{$eachSingletonCtgB} = 1;
								delete $singletonCtg->{$eachSingletonCtgA};
								delete $singletonCtg->{$eachSingletonCtgB};
								delete $singletonCtgId->{$eachSingletonCtgA};
								delete $singletonCtgId->{$eachSingletonCtgB};
								last;
							}						
						}
					}
				}

				my $multiCtg;
				my $multiCtgId;
				my $assemblyMultiCtgList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblyCtg' AND o = ? AND note LIKE '%,%'");
				$assemblyMultiCtgList->execute($assemblyId);
				while (my @assemblyMultiCtgList = $assemblyMultiCtgList->fetchrow_array())
				{
					$multiCtg->{$assemblyMultiCtgList[2]} = $assemblyMultiCtgList[8];
					$multiCtgId->{$assemblyMultiCtgList[2]} = $assemblyMultiCtgList[0];
				}

				#merge multies and singletons
				for my $eachSingletonCtg (sort {$a <=> $b} keys %$singletonCtg)
				{
					my $singletonCtgSeqId = $singletonCtg->{$eachSingletonCtg};
					$singletonCtgSeqId =~ s/[^a-zA-Z0-9]//g;
					my $eachSingletonCtgSeq=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
					$eachSingletonCtgSeq->execute($singletonCtgSeqId);
					my @eachSingletonCtgSeq = $eachSingletonCtgSeq->fetchrow_array();

					my $mergeFlag = 0;
					for my $eachMultiCtg (sort {$a <=> $b} keys %$multiCtg)
					{
						my $headSeqId = '';
						my $tailSeqId = '';
						for ( split ",", $multiCtg->{$eachMultiCtg})
						{
							next unless ($_);
							/^-/ and next;
							$_ =~ s/[^a-zA-Z0-9]//g;
							$headSeqId = $_ if ($headSeqId eq '');
							$tailSeqId = $_;
						}

						my $headCtgSeq=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
						$headCtgSeq->execute($headSeqId);
						my @headCtgSeq = $headCtgSeq->fetchrow_array();
						my $tailCtgSeq=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
						$tailCtgSeq->execute($tailSeqId);
						my @tailCtgSeq = $tailCtgSeq->fetchrow_array();

						my $queryDirHead;
						for (my $position = 0; $position < length($eachSingletonCtgSeq[5]); $position += 2)
						{
							$queryDirHead .= "/q". substr($eachSingletonCtgSeq[5],$position,2);
						}
						my $subjectDirHead;
						for (my $position = 0; $position < length($headCtgSeq[5]); $position += 2)
						{
							$subjectDirHead .= "/s". substr($headCtgSeq[5],$position,2);
						}
						if(-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirHead$subjectDirHead/$eachSingletonCtgSeq[5]-$headCtgSeq[5].tbl")
						{
							my $goodEndFoundA = 0;
							my $goodEndFoundB = 0;
							open (TBL,"$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirHead$subjectDirHead/$eachSingletonCtgSeq[5]-$headCtgSeq[5].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDirHead$subjectDirHead/$eachSingletonCtgSeq[5]-$headCtgSeq[5].tbl";
							while(<TBL>)
							{
								chop;
								/^#/ and next;
								my @getAlignmentHead = split("\t",$_);
								$goodEndFoundA = 1 if ($getAlignmentHead[6] == 1 || $getAlignmentHead[7] == $eachSingletonCtgSeq[6]);
								if($getAlignmentHead[8] < $getAlignmentHead[9])
								{
									$goodEndFoundB = 1 if ($getAlignmentHead[8] == 1 || $getAlignmentHead[9] == $headCtgSeq[6]);
								}
								else
								{
									$goodEndFoundB = 1 if ($getAlignmentHead[8] == $headCtgSeq[6] || $getAlignmentHead[9] == 1);
								}
							}
							close(TBL);
							if($goodEndFoundA && $goodEndFoundB)
							{
								my $mergedName = ($eachSingletonCtg < $eachMultiCtg) ? $eachSingletonCtg : $eachMultiCtg;
								$multiCtg->{$mergedName} = "$singletonCtg->{$eachSingletonCtg},$multiCtg->{$eachMultiCtg}";
								$multiCtgId->{$mergedName} = $multiCtgId->{$eachMultiCtg};
								my $updatedAssemblyCtg = $dbh->prepare("UPDATE matrix SET name = ?, note = ? WHERE id = ?");
								$updatedAssemblyCtg->execute($mergedName,$multiCtg->{$mergedName},$multiCtgId->{$mergedName});
								my $deleteAssemblyCtg=$dbh->do("DELETE FROM matrix WHERE id = $singletonCtgId->{$eachSingletonCtg}");
								if ($eachSingletonCtg < $eachMultiCtg)
								{
									delete $multiCtg->{$eachMultiCtg};
									delete $multiCtgId->{$eachMultiCtg};
								}
								delete $singletonCtg->{$eachSingletonCtg};
								delete $singletonCtgId->{$eachSingletonCtg};
								$mergeFlag = 1;
								last;
							}
						}
						last if ($mergeFlag);
						next if ($tailSeqId == $headSeqId);

						my $queryDirTail;
						for (my $position = 0; $position < length($eachSingletonCtgSeq[5]); $position += 2)
						{
							$queryDirTail .= "/q". substr($eachSingletonCtgSeq[5],$position,2);
						}
						my $subjectDirTail;
						for (my $position = 0; $position < length($tailCtgSeq[5]); $position += 2)
						{
							$subjectDirTail .= "/s". substr($tailCtgSeq[5],$position,2);
						}
						if(-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirTail$subjectDirTail/$eachSingletonCtgSeq[5]-$tailCtgSeq[5].tbl")
						{
							my $goodEndFoundA = 0;
							my $goodEndFoundB = 0;
							open (TBL,"$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirTail$subjectDirTail/$eachSingletonCtgSeq[5]-$tailCtgSeq[5].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDirTail$subjectDirTail/$eachSingletonCtgSeq[5]-$tailCtgSeq[5].tbl";
							while(<TBL>)
							{
								chop;
								/^#/ and next;
								my @getAlignmentTail = split("\t",$_);
								$goodEndFoundA = 1 if ($getAlignmentTail[6] == 1 || $getAlignmentTail[7] == $eachSingletonCtgSeq[6]);
								if($getAlignmentTail[8] < $getAlignmentTail[9])
								{
									$goodEndFoundB = 1 if ($getAlignmentTail[8] == 1 || $getAlignmentTail[9] == $tailCtgSeq[6]);
								}
								else
								{
									$goodEndFoundB = 1 if ($getAlignmentTail[8] == $tailCtgSeq[6] || $getAlignmentTail[9] == 1);
								}
							}
							close(TBL);
							if($goodEndFoundA && $goodEndFoundB)
							{
								my $mergedName = ($eachSingletonCtg < $eachMultiCtg) ? $eachSingletonCtg : $eachMultiCtg;
								$multiCtg->{$mergedName} = "$multiCtg->{$eachMultiCtg},$singletonCtg->{$eachSingletonCtg}";
								$multiCtgId->{$mergedName} = $multiCtgId->{$eachMultiCtg};
								my $updatedAssemblyCtg = $dbh->prepare("UPDATE matrix SET name = ?, note = ? WHERE id = ?");
								$updatedAssemblyCtg->execute($mergedName,$multiCtg->{$mergedName},$multiCtgId->{$mergedName});
								my $deleteAssemblyCtg=$dbh->do("DELETE FROM matrix WHERE id = $singletonCtgId->{$eachSingletonCtg}");
								if($eachSingletonCtg < $eachMultiCtg)
								{
									delete $multiCtgId->{$eachMultiCtg};								
									delete $multiCtg->{$eachMultiCtg};
								}
								delete $singletonCtg->{$eachSingletonCtg};
								delete $singletonCtgId->{$eachSingletonCtg};
								$mergeFlag = 1;									
								last;
							}
						}
						last if ($mergeFlag);
					}
				}
			
				#merge multies
				my $mergeMulties = 1;
				my $checkedCtg;
				do{
					my $mergeFlag = 0;
					for my $eachMultiCtgA (sort {$a <=> $b} keys %$multiCtg)
					{
						next if (exists $checkedCtg->{$eachMultiCtgA});
						$checkedCtg->{$eachMultiCtgA} = 1;

						my $headSeqIdA = '';
						my $tailSeqIdA = '';
						for ( split ",", $multiCtg->{$eachMultiCtgA})
						{
							next unless ($_);
							/^-/ and next;
							$_ =~ s/[^a-zA-Z0-9]//g;
							$headSeqIdA = $_ if ($headSeqIdA eq '');
							$tailSeqIdA = $_;
						}

						my $headCtgSeqA=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
						$headCtgSeqA->execute($headSeqIdA);
						my @headCtgSeqA = $headCtgSeqA->fetchrow_array();

						my $tailCtgSeqA=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
						$tailCtgSeqA->execute($tailSeqIdA);
						my @tailCtgSeqA = $tailCtgSeqA->fetchrow_array();

						for my $eachMultiCtgB (sort {$a <=> $b} keys %$multiCtg)
						{
							next if (exists $checkedCtg->{$eachMultiCtgB});

							my $headSeqIdB = '';
							my $tailSeqIdB = '';
							for ( split ",", $multiCtg->{$eachMultiCtgB})
							{
								next unless ($_);
								/^-/ and next;
								$_ =~ s/[^a-zA-Z0-9]//g;
								$headSeqIdB = $_ if ($headSeqIdB eq '');
								$tailSeqIdB = $_;
							}

							my $headCtgSeqB=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
							$headCtgSeqB->execute($headSeqIdB);
							my @headCtgSeqB = $headCtgSeqB->fetchrow_array();

							my $tailCtgSeqB=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
							$tailCtgSeqB->execute($tailSeqIdB);
							my @tailCtgSeqB = $tailCtgSeqB->fetchrow_array();

							my $queryDirAHBH;
							for (my $position = 0; $position < length($headCtgSeqA[5]); $position += 2)
							{
								$queryDirAHBH .= "/q". substr($headCtgSeqA[5],$position,2);
							}
							my $subjectDirAHBH;
							for (my $position = 0; $position < length($headCtgSeqB[5]); $position += 2)
							{
								$subjectDirAHBH .= "/s". substr($headCtgSeqB[5],$position,2);
							}
							if(-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirAHBH$subjectDirAHBH/$headCtgSeqA[5]-$headCtgSeqB[5].tbl")
							{
								my $goodEndFoundA = 0;
								my $goodEndFoundB = 0;
								open (TBL,"$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirAHBH$subjectDirAHBH/$headCtgSeqA[5]-$headCtgSeqB[5].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDirAHBH$subjectDirAHBH/$headCtgSeqA[5]-$headCtgSeqB[5].tbl";
								while(<TBL>)
								{
									chop;
									/^#/ and next;
									my @getAlignmentAHBH = split("\t",$_);
									$goodEndFoundA = 1 if ($getAlignmentAHBH[6] == 1 || $getAlignmentAHBH[7] == $headCtgSeqA[6]);
									if($getAlignmentAHBH[8] < $getAlignmentAHBH[9])
									{
										$goodEndFoundB = 1 if ($getAlignmentAHBH[8] == 1 || $getAlignmentAHBH[9] == $headCtgSeqB[6]);
									}
									else
									{
										$goodEndFoundB = 1 if ($getAlignmentAHBH[8] == $headCtgSeqB[6] || $getAlignmentAHBH[9] == 1);
									}
								}
								close(TBL);
								if($goodEndFoundA && $goodEndFoundB)
								{
									$multiCtg->{$eachMultiCtgA} = join ",", (reverse split ",", $multiCtg->{$eachMultiCtgA}); #flip ctgA
									foreach (split ",", $multiCtg->{$eachMultiCtgA})
									{
										next unless ($_);
										$_ =~ s/[^a-zA-Z0-9]//g;
										my $updateAssemblySeq=$dbh->do("UPDATE matrix SET barcode = barcode * (-1) WHERE id = $_");
									}
									if($eachMultiCtgA < $eachMultiCtgB)
									{
										$multiCtg->{$eachMultiCtgA} = "$multiCtg->{$eachMultiCtgA},$multiCtg->{$eachMultiCtgB}";
										my $updatedAssemblyCtg = $dbh->prepare("UPDATE matrix SET note = ? WHERE id = ?");
										$updatedAssemblyCtg->execute($multiCtg->{$eachMultiCtgA},$multiCtgId->{$eachMultiCtgA});
										my $deleteAssemblyCtg=$dbh->do("DELETE FROM matrix WHERE id = $multiCtgId->{$eachMultiCtgB}");
										delete $multiCtgId->{$eachMultiCtgB};								
										delete $multiCtg->{$eachMultiCtgB};
									}
									else
									{
										$multiCtg->{$eachMultiCtgB} = "$multiCtg->{$eachMultiCtgA},$multiCtg->{$eachMultiCtgB}";
										my $updatedAssemblyCtg = $dbh->prepare("UPDATE matrix SET note = ? WHERE id = ?");
										$updatedAssemblyCtg->execute($multiCtg->{$eachMultiCtgB},$multiCtgId->{$eachMultiCtgB});
										my $deleteAssemblyCtg=$dbh->do("DELETE FROM matrix WHERE id = $multiCtgId->{$eachMultiCtgA}");
										delete $multiCtgId->{$eachMultiCtgA};
										delete $multiCtg->{$eachMultiCtgA};
									}
									$mergeFlag = 1;
									last;
								}
							}
							last if($mergeFlag);

							my $queryDirAHBT;
							for (my $position = 0; $position < length($headCtgSeqA[5]); $position += 2)
							{
								$queryDirAHBT .= "/q". substr($headCtgSeqA[5],$position,2);
							}
							my $subjectDirAHBT;
							for (my $position = 0; $position < length($tailCtgSeqB[5]); $position += 2)
							{
								$subjectDirAHBT .= "/s". substr($tailCtgSeqB[5],$position,2);
							}
							if(-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirAHBT$subjectDirAHBT/$headCtgSeqA[5]-$tailCtgSeqB[5].tbl")
							{
								my $goodEndFoundA = 0;
								my $goodEndFoundB = 0;
								open (TBL,"$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirAHBT$subjectDirAHBT/$headCtgSeqA[5]-$tailCtgSeqB[5].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDirAHBT$subjectDirAHBT/$headCtgSeqA[5]-$tailCtgSeqB[5].tbl";
								while(<TBL>)
								{
									chop;
									/^#/ and next;
									my @getAlignmentAHBT = split("\t",$_);
									$goodEndFoundA = 1 if ($getAlignmentAHBT[6] == 1 || $getAlignmentAHBT[7] == $headCtgSeqA[6]);
									if($getAlignmentAHBT[10] < $getAlignmentAHBT[11])
									{
										$goodEndFoundB = 1 if ($getAlignmentAHBT[8] == 1 || $getAlignmentAHBT[9] == $tailCtgSeqB[6]);
									}
									else
									{
										$goodEndFoundB = 1 if ($getAlignmentAHBT[8] == $tailCtgSeqB[6] || $getAlignmentAHBT[9] == 1);
									}
								}
								close(TBL);
								if($goodEndFoundA && $goodEndFoundB)
								{
									if($eachMultiCtgA < $eachMultiCtgB)
									{
										$multiCtg->{$eachMultiCtgA} = "$multiCtg->{$eachMultiCtgB},$multiCtg->{$eachMultiCtgA}";
										my $updatedAssemblyCtg = $dbh->prepare("UPDATE matrix SET note = ? WHERE id = ?");
										$updatedAssemblyCtg->execute($multiCtg->{$eachMultiCtgA},$multiCtgId->{$eachMultiCtgA});
										my $deleteAssemblyCtg=$dbh->do("DELETE FROM matrix WHERE id = $multiCtgId->{$eachMultiCtgB}");
										delete $multiCtgId->{$eachMultiCtgB};								
										delete $multiCtg->{$eachMultiCtgB};
									}
									else
									{
										$multiCtg->{$eachMultiCtgB} = "$multiCtg->{$eachMultiCtgB},$multiCtg->{$eachMultiCtgA}";
										my $updatedAssemblyCtg = $dbh->prepare("UPDATE matrix SET note = ? WHERE id = ?");
										$updatedAssemblyCtg->execute($multiCtg->{$eachMultiCtgB},$multiCtgId->{$eachMultiCtgB});
										my $deleteAssemblyCtg=$dbh->do("DELETE FROM matrix WHERE id = $multiCtgId->{$eachMultiCtgA}");
										delete $multiCtgId->{$eachMultiCtgA};
										delete $multiCtg->{$eachMultiCtgA};
									}
									$mergeFlag = 1;
									last;
								}
							}
							last if($mergeFlag);

							my $queryDirATBH;
							for (my $position = 0; $position < length($tailCtgSeqA[5]); $position += 2)
							{
								$queryDirATBH .= "/q". substr($tailCtgSeqA[5],$position,2);
							}
							my $subjectDirATBH;
							for (my $position = 0; $position < length($headCtgSeqB[5]); $position += 2)
							{
								$subjectDirATBH .= "/s". substr($headCtgSeqB[5],$position,2);
							}
							if(-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirATBH$subjectDirATBH/$tailCtgSeqA[5]-$headCtgSeqB[5].tbl")
							{
								my $goodEndFoundA = 0;
								my $goodEndFoundB = 0;
								open (TBL,"$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirATBH$subjectDirATBH/$tailCtgSeqA[5]-$headCtgSeqB[5].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDirATBH$subjectDirATBH/$tailCtgSeqA[5]-$headCtgSeqB[5].tbl";
								while(<TBL>)
								{
									chop;
									/^#/ and next;
									my @getAlignmentATBH = split("\t",$_);
									$goodEndFoundA = 1 if ($getAlignmentATBH[6] == 1 || $getAlignmentATBH[7] == $tailCtgSeqA[6]);
									if($getAlignmentATBH[10] < $getAlignmentATBH[11])
									{
										$goodEndFoundB = 1 if ($getAlignmentATBH[8] == 1 || $getAlignmentATBH[9] == $headCtgSeqB[6]);
									}
									else
									{
										$goodEndFoundB = 1 if ($getAlignmentATBH[8] == $headCtgSeqB[6] || $getAlignmentATBH[9] == 1);
									}
								}
								close(TBL);
								if($goodEndFoundA && $goodEndFoundB)
								{
									if($eachMultiCtgA < $eachMultiCtgB)
									{
										$multiCtg->{$eachMultiCtgA} = "$multiCtg->{$eachMultiCtgA},$multiCtg->{$eachMultiCtgB}";
										my $updatedAssemblyCtg = $dbh->prepare("UPDATE matrix SET note = ? WHERE id = ?");
										$updatedAssemblyCtg->execute($multiCtg->{$eachMultiCtgA},$multiCtgId->{$eachMultiCtgA});
										my $deleteAssemblyCtg=$dbh->do("DELETE FROM matrix WHERE id = $multiCtgId->{$eachMultiCtgB}");
										delete $multiCtgId->{$eachMultiCtgB};								
										delete $multiCtg->{$eachMultiCtgB};
									}
									else
									{
										$multiCtg->{$eachMultiCtgB} = "$multiCtg->{$eachMultiCtgA},$multiCtg->{$eachMultiCtgB}";
										my $updatedAssemblyCtg = $dbh->prepare("UPDATE matrix SET note = ? WHERE id = ?");
										$updatedAssemblyCtg->execute($multiCtg->{$eachMultiCtgB},$multiCtgId->{$eachMultiCtgB});
										my $deleteAssemblyCtg=$dbh->do("DELETE FROM matrix WHERE id = $multiCtgId->{$eachMultiCtgA}");
										delete $multiCtgId->{$eachMultiCtgA};
										delete $multiCtg->{$eachMultiCtgA};
									}
									$mergeFlag = 1;
									last;
								}
							}
							last if($mergeFlag);

							my $queryDirATBT;
							for (my $position = 0; $position < length($tailCtgSeqA[5]); $position += 2)
							{
								$queryDirATBT .= "/q". substr($tailCtgSeqA[5],$position,2);
							}
							my $subjectDirATBT;
							for (my $position = 0; $position < length($tailCtgSeqB[5]); $position += 2)
							{
								$subjectDirATBT .= "/s". substr($tailCtgSeqB[5],$position,2);
							}
							if(-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirATBT$subjectDirATBT/$tailCtgSeqA[5]-$tailCtgSeqB[5].tbl")
							{
								my $goodEndFoundA = 0;
								my $goodEndFoundB = 0;
								open (TBL,"$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirATBT$subjectDirATBT/$tailCtgSeqA[5]-$tailCtgSeqB[5].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDirATBT$subjectDirATBT/$tailCtgSeqA[5]-$tailCtgSeqB[5].tbl";
								while(<TBL>)
								{
									chop;
									/^#/ and next;
									my @getAlignmentATBT = split("\t",$_);
									$goodEndFoundA = 1 if ($getAlignmentATBT[6] == 1 || $getAlignmentATBT[7] == $tailCtgSeqA[6]);
									if($getAlignmentATBT[10] < $getAlignmentATBT[11])
									{
										$goodEndFoundB = 1 if ($getAlignmentATBT[8] == 1 || $getAlignmentATBT[9] == $tailCtgSeqB[6]);
									}
									else
									{
										$goodEndFoundB = 1 if ($getAlignmentATBT[8] == $tailCtgSeqB[6] || $getAlignmentATBT[9] == 1);
									}
								}
								close(TBL);
								if($goodEndFoundA && $goodEndFoundB)
								{
									$multiCtg->{$eachMultiCtgB} = join ",", (reverse split ",", $multiCtg->{$eachMultiCtgB}); #flip ctgB
									foreach (split ",", $multiCtg->{$eachMultiCtgB})
									{
										next unless ($_);
										$_ =~ s/[^a-zA-Z0-9]//g;
										my $updateAssemblySeq=$dbh->do("UPDATE matrix SET barcode = barcode * (-1) WHERE id = $_");
									}
									if($eachMultiCtgA < $eachMultiCtgB)
									{
										$multiCtg->{$eachMultiCtgA} = "$multiCtg->{$eachMultiCtgA},$multiCtg->{$eachMultiCtgB}";
										my $updatedAssemblyCtg = $dbh->prepare("UPDATE matrix SET note = ? WHERE id = ?");
										$updatedAssemblyCtg->execute($multiCtg->{$eachMultiCtgA},$multiCtgId->{$eachMultiCtgA});
										my $deleteAssemblyCtg=$dbh->do("DELETE FROM matrix WHERE id = $multiCtgId->{$eachMultiCtgB}");
										delete $multiCtgId->{$eachMultiCtgB};								
										delete $multiCtg->{$eachMultiCtgB};
									}
									else
									{
										$multiCtg->{$eachMultiCtgB} = "$multiCtg->{$eachMultiCtgA},$multiCtg->{$eachMultiCtgB}";
										my $updatedAssemblyCtg = $dbh->prepare("UPDATE matrix SET note = ? WHERE id = ?");
										$updatedAssemblyCtg->execute($multiCtg->{$eachMultiCtgB},$multiCtgId->{$eachMultiCtgB});
										my $deleteAssemblyCtg=$dbh->do("DELETE FROM matrix WHERE id = $multiCtgId->{$eachMultiCtgA}");
										delete $multiCtgId->{$eachMultiCtgA};
										delete $multiCtg->{$eachMultiCtgA};
									}							
									$mergeFlag = 1;
									delete $checkedCtg->{$eachMultiCtgA};
									last;
								}
							}
							last if($mergeFlag);
						}
						delete $checkedCtg->{$eachMultiCtgA} if($mergeFlag);
						last if($mergeFlag);
					}
					unless($mergeFlag)
					{
						$mergeMulties = 0;
					}
				} while ($mergeMulties);

			}
			my $updateAssemblyToRunningStatus=$dbh->do("UPDATE matrix SET barcode = '-1' WHERE id = $assemblyId");
		}

		if($redundancyFilterOverlap)
		{
			my $updateAssemblyToFilteringOverlap=$dbh->do("UPDATE matrix SET barcode = '-10' WHERE id = $assemblyId");
			my $assemblyMultiCtgList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblyCtg' AND o = ? AND note LIKE '%,%'");
			$assemblyMultiCtgList->execute($assemblyId);
			while (my @assemblyMultiCtgList = $assemblyMultiCtgList->fetchrow_array())
			{
				my @seqInCtg;
				foreach (split ",", $assemblyMultiCtgList[8])
				{
					next unless ($_);
					/^-/ and next;
					$_ =~ s/[^a-zA-Z0-9]//g;
					push @seqInCtg,$_;
				}
				my $totalSeqInCtg = @seqInCtg;
				for (@seqInCtg)
				{
					my $assemblySeqResetLength=$dbh->do("UPDATE matrix SET note = CONCAT('1,', z) WHERE id = $_");
				}
				for (my $i = 1; $i < $totalSeqInCtg; $i++)
				{
					my $preAssemblySeq=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
					$preAssemblySeq->execute($seqInCtg[$i-1]);
					my @preAssemblySeq = $preAssemblySeq->fetchrow_array();
					my $preSeqStart;
					my $preSeqEnd;
					if($preAssemblySeq[8])
					{
						($preSeqStart,$preSeqEnd) = split ",",$preAssemblySeq[8];
					}
					else
					{
						$preSeqStart = 1;
						$preSeqEnd = $preAssemblySeq[6];
					}
					my $preSequence = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
					$preSequence->execute($preAssemblySeq[5]);
					my @preSequence = $preSequence->fetchrow_array();

					my $nextAssemblySeq=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
					$nextAssemblySeq->execute($seqInCtg[$i]);
					my @nextAssemblySeq = $nextAssemblySeq->fetchrow_array();
					my $nextSeqStart;
					my $nextSeqEnd;
					if($nextAssemblySeq[8])
					{
						($nextSeqStart,$nextSeqEnd) = split ",",$nextAssemblySeq[8];
					}
					else
					{
						$nextSeqStart = 1;
						$nextSeqEnd = $nextAssemblySeq[6];
					}
					my $nextSequence = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
					$nextSequence->execute($nextAssemblySeq[5]);
					my @nextSequence = $nextSequence->fetchrow_array();

					my $queryDir;
					for (my $position = 0; $position < length($preAssemblySeq[5]); $position += 2)
					{
						$queryDir .= "/q". substr($preAssemblySeq[5],$position,2);
					}
					my $subjectDir;
					for (my $position = 0; $position < length($nextAssemblySeq[5]); $position += 2)
					{
						$subjectDir .= "/s". substr($nextAssemblySeq[5],$position,2);
					}
					
					if(-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$preAssemblySeq[5]-$nextAssemblySeq[5].tbl")
					{
						my @getAlignment;
						my $alignmentCount = 0;
						open (TBL,"$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$preAssemblySeq[5]-$nextAssemblySeq[5].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$preAssemblySeq[5]-$nextAssemblySeq[5].tbl";
						while(<TBL>)
						{
							chop;
							/^#/ and next;
							$alignmentCount++;
							@getAlignment = split("\t",$_);
							last if ($alignmentCount > 0);
						}
						close(TBL);

						my $preSeqStartCandidate = $preSeqStart;
						my $preSeqEndCandidate = $preSeqEnd;
						my $nextSeqStartCandidate = $nextSeqStart;
						my $nextSeqEndCandidate = $nextSeqEnd;

						if($preSequence[3] < 4) #if preSeq has non-gapped sequence
						{
							#keep preSeq sequence
							if($preAssemblySeq[7] > 0)
							{
								if ($getAlignment[8] < $getAlignment[9])
								{
									$preSeqEndCandidate = $getAlignment[7];
									$nextSeqStartCandidate = $getAlignment[9] + 1;
								}
								else
								{
									$preSeqEndCandidate = $getAlignment[7];
									$nextSeqEndCandidate = $getAlignment[9] - 1;
								}
							}
							else
							{
								if ($getAlignment[8] < $getAlignment[9])
								{
									$preSeqStartCandidate = $getAlignment[6];
									$nextSeqEndCandidate = $getAlignment[8] - 1;
								}
								else
								{
									$preSeqStartCandidate = $getAlignment[6];
									$nextSeqStartCandidate = $getAlignment[8] + 1;
								}
							}

							if ($preSeqStartCandidate >= 1 && $preSeqEndCandidate <= $preAssemblySeq[6] && $nextSeqStartCandidate >= 1 && $nextSeqEndCandidate <= $nextAssemblySeq[6] && $preSeqStartCandidate <= $preSeqEndCandidate && $nextSeqStartCandidate <= $nextSeqEndCandidate)
							{
								my $updateAssemblySeqPre=$dbh->do("UPDATE matrix SET note = '$preSeqStartCandidate,$preSeqEndCandidate' WHERE id = $preAssemblySeq[0]");
								my $updateAssemblySeqNext=$dbh->do("UPDATE matrix SET note = '$nextSeqStartCandidate,$nextSeqEndCandidate' WHERE id = $nextAssemblySeq[0]");
							}
							else
							{
								#keep nextSeq sequence
								$preSeqStartCandidate = $preSeqStart;
								$preSeqEndCandidate = $preSeqEnd;
								$nextSeqStartCandidate = $nextSeqStart;
								$nextSeqEndCandidate = $nextSeqEnd;
								if($preAssemblySeq[7] > 0)
								{
									if ($getAlignment[8] < $getAlignment[9])
									{
										$preSeqEndCandidate = $getAlignment[6] - 1;
										$nextSeqStartCandidate = $getAlignment[8];
									}
									else
									{
										$preSeqEndCandidate = $getAlignment[6] - 1;
										$nextSeqEndCandidate = $getAlignment[8];
									}
								}
								else
								{
									if ($getAlignment[8] < $getAlignment[9])
									{
										$preSeqStartCandidate = $getAlignment[7] + 1;
										$nextSeqEndCandidate = $getAlignment[9];
									}
									else
									{
										$preSeqStartCandidate = $getAlignment[7] + 1;
										$nextSeqStartCandidate = $getAlignment[9];
									}
								}

								if ($preSeqStartCandidate >= 1 && $preSeqEndCandidate <= $preAssemblySeq[6] && $nextSeqStartCandidate >= 1 && $nextSeqEndCandidate <= $nextAssemblySeq[6] && $preSeqStartCandidate <= $preSeqEndCandidate && $nextSeqStartCandidate <= $nextSeqEndCandidate)
								{
									my $updateAssemblySeqPre=$dbh->do("UPDATE matrix SET note = '$preSeqStartCandidate,$preSeqEndCandidate' WHERE id = $preAssemblySeq[0]");
									my $updateAssemblySeqNext=$dbh->do("UPDATE matrix SET note = '$nextSeqStartCandidate,$nextSeqEndCandidate' WHERE id = $nextAssemblySeq[0]");
								}
								else
								{
									#do nothing if conflicts found
								}
							}
						}
						else
						{
							#keep nextSeq sequence
							if($preAssemblySeq[7] > 0)
							{
								if ($getAlignment[8] < $getAlignment[9])
								{
									$preSeqEndCandidate = $getAlignment[6] - 1;
									$nextSeqStartCandidate = $getAlignment[8];
								}
								else
								{
									$preSeqEndCandidate = $getAlignment[6] - 1;
									$nextSeqEndCandidate = $getAlignment[8];
								}
							}
							else
							{
								if ($getAlignment[8] < $getAlignment[9])
								{
									$preSeqStartCandidate = $getAlignment[7] + 1;
									$nextSeqEndCandidate = $getAlignment[9];
								}
								else
								{
									$preSeqStartCandidate = $getAlignment[7] + 1;
									$nextSeqStartCandidate = $getAlignment[9];
								}
							}

							if ($preSeqStartCandidate >= 1 && $preSeqEndCandidate <= $preAssemblySeq[6] && $nextSeqStartCandidate >= 1 && $nextSeqEndCandidate <= $nextAssemblySeq[6] && $preSeqStartCandidate <= $preSeqEndCandidate && $nextSeqStartCandidate <= $nextSeqEndCandidate)
							{
								my $updateAssemblySeqPre=$dbh->do("UPDATE matrix SET note = '$preSeqStartCandidate,$preSeqEndCandidate' WHERE id = $preAssemblySeq[0]");
								my $updateAssemblySeqNext=$dbh->do("UPDATE matrix SET note = '$nextSeqStartCandidate,$nextSeqEndCandidate' WHERE id = $nextAssemblySeq[0]");
							}
							else
							{
								$preSeqStartCandidate = $preSeqStart;
								$preSeqEndCandidate = $preSeqEnd;
								$nextSeqStartCandidate = $nextSeqStart;
								$nextSeqEndCandidate = $nextSeqEnd;
								#keep preSeq sequence
								if($preAssemblySeq[7] > 0)
								{
									if ($getAlignment[8] < $getAlignment[9])
									{
										$preSeqEndCandidate = $getAlignment[7];
										$nextSeqStartCandidate = $getAlignment[9] + 1;
									}
									else
									{
										$preSeqEndCandidate = $getAlignment[7];
										$nextSeqEndCandidate = $getAlignment[9] - 1;
									}
								}
								else
								{
									if ($getAlignment[8] < $getAlignment[9])
									{
										$preSeqStartCandidate = $getAlignment[6];
										$nextSeqEndCandidate = $getAlignment[8] - 1;
									}
									else
									{
										$preSeqStartCandidate = $getAlignment[6];
										$nextSeqStartCandidate = $getAlignment[8] + 1;
									}
								}

								if ($preSeqStartCandidate >= 1 && $preSeqEndCandidate <= $preAssemblySeq[6] && $nextSeqStartCandidate >= 1 && $nextSeqEndCandidate <= $nextAssemblySeq[6] && $preSeqStartCandidate <= $preSeqEndCandidate && $nextSeqStartCandidate <= $nextSeqEndCandidate)
								{
									my $updateAssemblySeqPre=$dbh->do("UPDATE matrix SET note = '$preSeqStartCandidate,$preSeqEndCandidate' WHERE id = $preAssemblySeq[0]");
									my $updateAssemblySeqNext=$dbh->do("UPDATE matrix SET note = '$nextSeqStartCandidate,$nextSeqEndCandidate' WHERE id = $nextAssemblySeq[0]");
								}
								else
								{
									#do nothing if conflicts found
								}
							}
						}
					}
				}
			}
			my $updateAssemblyToRunningStatus=$dbh->do("UPDATE matrix SET barcode = '-1' WHERE id = $assemblyId");
		}

		if($refGenomeId)
		{
			my $updateAssemblyGenomeId=$dbh->do("UPDATE matrix SET y = $refGenomeId WHERE id = $assemblyId");
			my $sequenceInRefGenome;
			my $getGenomeSequence = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 99 AND x = ?");
			$getGenomeSequence->execute($refGenomeId);
			while(my @getGenomeSequence = $getGenomeSequence->fetchrow_array())
			{
				$sequenceInRefGenome->{$getGenomeSequence[0]} = $getGenomeSequence[6];
			}
			
			my $updateAssemblyToAssignChr=$dbh->do("UPDATE matrix SET barcode = '-5' WHERE id = $assemblyId");
			#assign ctg to genome chr
			my $assemblyAllCtgList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblyCtg' AND o = ?");
			$assemblyAllCtgList->execute($assemblyId);
			while (my @assemblyAllCtgList = $assemblyAllCtgList->fetchrow_array())
			{
				my $assemblyCtgLength = 0;
				my $assemblySeqList = '';
				my $chrNumber;
				my $chrPosition;

				my @assemblySeqListAll;
				my $firstAssemblySeq = "";
				my $lastAssemblySeq = "";
			
				foreach (split ",", $assemblyAllCtgList[8])
				{
					next unless ($_);
					$assemblySeqList .= ($assemblySeqList ne '') ? ",$_": $_;
					/^-/ and next;
					$_ =~ s/[^a-zA-Z0-9]//g;
					push @assemblySeqListAll, $_;
					$firstAssemblySeq = $_ unless ($firstAssemblySeq);
					$lastAssemblySeq = $_;
				}
				my $lastComponentType = '';
				for my $assemblySeqIndividual (@assemblySeqListAll)
				{
					my $queryDir;
					for (my $position = 0; $position < length($sequenceIdOfAssemblySeq->{$assemblySeqIndividual}); $position += 2)
					{
						$queryDir .= "/q". substr($sequenceIdOfAssemblySeq->{$assemblySeqIndividual},$position,2);
					}

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
									open (TBL, "$dir/$file") or die "can't open file: $dir/$file";
									while(<TBL>)
									{
										chop;
										/^#/ and next;
										my @getAlignment = split("\t",$_);
										next if ($getAlignment[12] > 0);
										next unless (exists $sequenceInRefGenome->{$getAlignment[1]}); #check if subject is in refGenome
										next unless ($getAlignment[3] >= $alignmentBlockSize || $getAlignment[3]*100/$sequenceLength->{$sequenceIdOfAssemblySeq->{$assemblySeqIndividual}} >= $alignmentCoverage);
										$chrNumber->{$sequenceInRefGenome->{$getAlignment[1]}} = 0 unless (exists $chrNumber->{$sequenceInRefGenome->{$getAlignment[1]}});
										$chrNumber->{$sequenceInRefGenome->{$getAlignment[1]}} += $getAlignment[3];
										my $estimatedPosition = ($getAlignment[8] < $getAlignment[9]) ? $getAlignment[8] - $getAlignment[6] - $assemblyCtgLength : $getAlignment[9] - $sequenceLength->{$sequenceIdOfAssemblySeq->{$assemblySeqIndividual}} + $getAlignment[7] - $assemblyCtgLength;
										if(exists $chrPosition->{$sequenceInRefGenome->{$getAlignment[1]}})
										{
											$chrPosition->{$sequenceInRefGenome->{$getAlignment[1]}} .= ",$estimatedPosition";
										}
										else
										{
											$chrPosition->{$sequenceInRefGenome->{$getAlignment[1]}} = $estimatedPosition;
										}
									}
									close(TBL);
								}
							}
							elsif(-d "$dir/$file")
							{
								unshift @dirs, "$dir/$file";
							}
						}
					}while (@dirs);


					my $assemblySeq=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
					$assemblySeq->execute($assemblySeqIndividual);
					my @assemblySeq = $assemblySeq->fetchrow_array();
					my $assemblySeqStart;
					my $assemblySeqEnd;
					if($assemblySeq[8])
					{
						($assemblySeqStart,$assemblySeqEnd) = split ",",$assemblySeq[8];
					}
					else
					{
						$assemblySeqStart = 1;
						$assemblySeqEnd = $assemblySeq[6];
						my $updateAssemblySeq=$dbh->do("UPDATE matrix SET note = '$assemblySeqStart,$assemblySeqEnd' WHERE id = $assemblySeq[0]");
					}
					my $filterLength = 0;
					my $sequence = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
					$sequence->execute($assemblySeq[5]);
					my @sequence =  $sequence->fetchrow_array();
					my $sequenceDetails = decode_json $sequence[8];
					$sequenceDetails->{'filter'} = '' unless (exists $sequenceDetails->{'filter'});
					if ($sequenceDetails->{'filter'}) 
					{
						foreach (split ",", $sequenceDetails->{'filter'} )
						{
							my ($filterStart,$filterEnd) = split "-", $_;
							next if ($assemblySeqStart > $filterEnd);
							next if ($assemblySeqEnd < $filterStart);
							if ($assemblySeqStart >= $filterStart && $assemblySeqEnd <= $filterEnd)
							{
								$filterLength += $assemblySeqEnd - $assemblySeqStart + 1;
							}
							elsif ($assemblySeqStart >= $filterStart && $assemblySeqStart <= $filterEnd)
							{
								$filterLength += $filterEnd - $assemblySeqStart + 1;
							}
							elsif ($assemblySeqEnd >= $filterStart && $assemblySeqEnd <= $filterEnd)
							{
								$filterLength += $assemblySeqEnd - $filterStart + 1;
							}
							else
							{
								$filterLength += $filterEnd - $filterStart + 1;
							}
						}
					}
					#add non-end gaps to assemblyCtg length
					if($assemblySeqIndividual ne $firstAssemblySeq && $lastComponentType ne 'U')
					{
						if ($assemblySeq[4] eq 1 || $assemblySeq[4] eq 3 || $assemblySeq[4] eq 4 || $assemblySeq[4] eq 6 || $assemblySeq[4] eq 7 || $assemblySeq[4] eq 8) # add 5' 100 Ns
						{
							$assemblyCtgLength += $gapLength;
							$lastComponentType = 'U';
						}
					}
					$assemblyCtgLength += $assemblySeqEnd - $assemblySeqStart + 1 - $filterLength;
					$lastComponentType = 'D';
					if($assemblySeqIndividual ne $lastAssemblySeq && $lastComponentType ne 'U')
					{
						if ($assemblySeq[4] eq 2 || $assemblySeq[4] eq 3 || $assemblySeq[4] eq 5 || $assemblySeq[4] eq 6 || $assemblySeq[4] eq 7 || $assemblySeq[4] eq 8) # add 5' 100 Ns
						{
							$assemblyCtgLength += $gapLength;				
							$lastComponentType = 'U';
						}
					}
				}
				my @assignedChr = sort {$chrNumber->{$b} <=> $chrNumber->{$a}} keys %$chrNumber;
				my $assignedChr = ($assemblyAllCtgList[4]) ? $assemblyAllCtgList[4] : (@assignedChr) ? shift @assignedChr : 0;
				my $assignedPosition = 0;
					
				if($assignedChr)
				{
					my @assignedPositionCandidates = split ",",$chrPosition->{$assignedChr};
					sort {$a <=> $b} @assignedPositionCandidates;
					my $estimatedMedian = int ($#assignedPositionCandidates/2);
					$assignedPosition = $assignedPositionCandidates[$estimatedMedian];
				}
				if($assignChr && $assemblyAllCtgList[4] == 0) #assign chrNumber for unplaced contigs
				{
					my $updateAssemblyCtg=$dbh->do("UPDATE matrix SET x = $assignedChr, z = $assignedPosition WHERE id = $assemblyAllCtgList[0]");
				}
				if ($reposition && $assemblyAllCtgList[4]) #relocate position
				{
					my $updateAssemblyCtg=$dbh->do("UPDATE matrix SET x = $assignedChr, z = $assignedPosition WHERE id = $assemblyAllCtgList[0]");
				}
				my $updateAssemblyCtgLength=$dbh->do("UPDATE matrix SET barcode = $assemblyCtgLength, note = '$assemblySeqList' WHERE id = $assemblyAllCtgList[0]");
			}
			my $updateAssemblyToRunningStatus=$dbh->do("UPDATE matrix SET barcode = '-1' WHERE id = $assemblyId");

			if($orientContigs)
			{
				my $updateAssemblyToOrientingContigs=$dbh->do("UPDATE matrix SET barcode = '-8' WHERE id = $assemblyId");
				my $assemblyAllCtgList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblyCtg' AND o = ? ORDER BY x,z");
				$assemblyAllCtgList->execute($assemblyId);
				while (my @assemblyAllCtgList = $assemblyAllCtgList->fetchrow_array())
				{
					my $toBeFlipped = 0;
					for (split ",", $assemblyAllCtgList[8])
					{
						next unless ($_);
						/^-/ and next;
						$_ =~ s/[^a-zA-Z0-9]//g;						
						my $assemblySeq=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
						$assemblySeq->execute($_);
						my @assemblySeq = $assemblySeq->fetchrow_array();

						my $queryDir;
						for (my $position = 0; $position < length($assemblySeq[5]); $position += 2)
						{
							$queryDir .= "/q". substr($assemblySeq[5],$position,2);
						}

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
										open (TBL, "$dir/$file") or die "can't open file: $dir/$file";
										while(<TBL>)
										{
											chop;
											/^#/ and next;
											my @getAlignment = split("\t",$_);
											next if ($getAlignment[12] > 0);
											next unless (exists $sequenceInRefGenome->{$getAlignment[1]});
											if($assemblySeq[7] > 0)
											{
												if ($getAlignment[8] > $getAlignment[9])
												{
													$toBeFlipped += $getAlignment[3];
												}
												else
												{
													$toBeFlipped -= $getAlignment[3];
												}
											}
											else
											{
												if ($getAlignment[8] < $getAlignment[9])
												{
													$toBeFlipped += $getAlignment[3];
												}
												else
												{
													$toBeFlipped -= $getAlignment[3];
												}
											}
										}
										close(TBL);
									}
								}
								elsif(-d "$dir/$file")
								{
									unshift @dirs, "$dir/$file";
								}
							}
						}while (@dirs);

					}
					if ($toBeFlipped > 0)
					{
						$assemblyAllCtgList[8] = join ",", (reverse split ",", $assemblyAllCtgList[8]);
						foreach (split ",", $assemblyAllCtgList[8])
						{
							next unless ($_);
							$_ =~ s/[^a-zA-Z0-9]//g;
							my $updateAssemblySeq=$dbh->do("UPDATE matrix SET barcode = barcode * (-1) WHERE id = $_");
						}
						my $updateAssemblyCtg=$dbh->prepare("UPDATE matrix SET note = ? WHERE id = ?");
						$updateAssemblyCtg->execute($assemblyAllCtgList[8],$assemblyAllCtgList[0]);
					}
				}
				my $updateAssemblyToRunningStatus=$dbh->do("UPDATE matrix SET barcode = '-1' WHERE id = $assemblyId");
 	 		}

			my $updateAssemblyToEstimatingLength=$dbh->do("UPDATE matrix SET barcode = '-11' WHERE id = $assemblyId");
			if($renumber)
			{
				my $assemblyAllCtgList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblyCtg' AND o = ? ORDER BY x,z,barcode");
				$assemblyAllCtgList->execute($assemblyId);
				my $totalCtgNumber = $assemblyAllCtgList->rows;
				my $assignedChrOrder = 1;
				while (my @assemblyAllCtgList = $assemblyAllCtgList->fetchrow_array())
				{
					if($assemblyAllCtgList[4] > 0)
					{
						my $updateAssemblyCtgLength=$dbh->do("UPDATE matrix SET name = $assignedChrOrder, y = $assignedChrOrder WHERE id = $assemblyAllCtgList[0]");
						$assignedChrOrder++;
					}
					else
					{
						my $updateAssemblyCtgLength=$dbh->do("UPDATE matrix SET name = $totalCtgNumber, y = $totalCtgNumber WHERE id = $assemblyAllCtgList[0]");
						$totalCtgNumber--;
					}
				}
			}
		}
		else
		{
			my $updateAssemblyToEstimatingLength=$dbh->do("UPDATE matrix SET barcode = '-11' WHERE id = $assemblyId");
			my $assemblyAllCtgList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblyCtg' AND o = ? ORDER BY x,z,barcode");
			$assemblyAllCtgList->execute($assemblyId);
			my $totalCtgNumber = $assemblyAllCtgList->rows;
			my $assignedChrOrder = 1;
			while (my @assemblyAllCtgList = $assemblyAllCtgList->fetchrow_array())
			{
				my $assemblyCtgLength = 0;
				my $assemblySeqList = '';

				my @assemblySeqListAll;
				my $firstAssemblySeq = "";
				my $lastAssemblySeq = "";
			
				foreach (split ",", $assemblyAllCtgList[8])
				{
					next unless ($_);
					$assemblySeqList .= ($assemblySeqList ne '') ? ",$_": $_;
					/^-/ and next;
					$_ =~ s/[^a-zA-Z0-9]//g;
					push @assemblySeqListAll, $_;
					$firstAssemblySeq = $_ unless ($firstAssemblySeq);
					$lastAssemblySeq = $_;
				}

				my $lastComponentType = '';
				for my $assemblySeqIndividual (@assemblySeqListAll)
				{
					my $assemblySeq=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
					$assemblySeq->execute($assemblySeqIndividual);
					my @assemblySeq = $assemblySeq->fetchrow_array();
					my $assemblySeqStart;
					my $assemblySeqEnd;
					if($assemblySeq[8])
					{
						($assemblySeqStart,$assemblySeqEnd) = split ",",$assemblySeq[8];
					}
					else
					{
						$assemblySeqStart = 1;
						$assemblySeqEnd = $assemblySeq[6];
						my $updateAssemblySeq=$dbh->do("UPDATE matrix SET note = '$assemblySeqStart,$assemblySeqEnd' WHERE id = $assemblySeq[0]");
					}
					my $filterLength = 0;
					my $sequence = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
					$sequence->execute($assemblySeq[5]);
					my @sequence =  $sequence->fetchrow_array();
					my $sequenceDetails = decode_json $sequence[8];
					$sequenceDetails->{'filter'} = '' unless (exists $sequenceDetails->{'filter'});
					if ($sequenceDetails->{'filter'}) 
					{
						foreach (split ",", $sequenceDetails->{'filter'} )
						{
							my ($filterStart,$filterEnd) = split "-", $_;
							next if ($assemblySeqStart > $filterEnd);
							next if ($assemblySeqEnd < $filterStart);
							if ($assemblySeqStart >= $filterStart && $assemblySeqEnd <= $filterEnd)
							{
								$filterLength += $assemblySeqEnd - $assemblySeqStart + 1;
							}
							elsif ($assemblySeqStart >= $filterStart && $assemblySeqStart <= $filterEnd)
							{
								$filterLength += $filterEnd - $assemblySeqStart + 1;
							}
							elsif ($assemblySeqEnd >= $filterStart && $assemblySeqEnd <= $filterEnd)
							{
								$filterLength += $assemblySeqEnd - $filterStart + 1;
							}
							else
							{
								$filterLength += $filterEnd - $filterStart + 1;
							}
						}
					}
					#add non-end gaps to assemblyCtg length
					if($assemblySeqIndividual ne $firstAssemblySeq && $lastComponentType ne 'U')
					{
						if ($assemblySeq[4] eq 1 || $assemblySeq[4] eq 3 || $assemblySeq[4] eq 4 || $assemblySeq[4] eq 6 || $assemblySeq[4] eq 7 || $assemblySeq[4] eq 8) # add 5' 100 Ns
						{
							$assemblyCtgLength += $gapLength;
							$lastComponentType = 'U';
						}
					}
					$assemblyCtgLength += $assemblySeqEnd - $assemblySeqStart + 1 - $filterLength;
					$lastComponentType = 'D';
					if($assemblySeqIndividual ne $lastAssemblySeq && $lastComponentType ne 'U')
					{
						if ($assemblySeq[4] eq 2 || $assemblySeq[4] eq 3 || $assemblySeq[4] eq 5 || $assemblySeq[4] eq 6 || $assemblySeq[4] eq 7 || $assemblySeq[4] eq 8) # add 5' 100 Ns
						{
							$assemblyCtgLength += $gapLength;				
							$lastComponentType = 'U';
						}
					}
				}

				if($renumber)
				{
					if($assemblyAllCtgList[4] > 0)
					{
						my $updateAssemblyCtgLength=$dbh->do("UPDATE matrix SET name = $assignedChrOrder, y = $assignedChrOrder, barcode = $assemblyCtgLength, note = '$assemblySeqList' WHERE id = $assemblyAllCtgList[0]");
						$assignedChrOrder++;
					}
					else
					{
						my $updateAssemblyCtgLength=$dbh->do("UPDATE matrix SET name = $totalCtgNumber, y = $totalCtgNumber, barcode = $assemblyCtgLength, note = '$assemblySeqList' WHERE id = $assemblyAllCtgList[0]");
						$totalCtgNumber--;
					}
				}
				else
				{
					my $updateAssemblyCtgLength=$dbh->do("UPDATE matrix SET y = $assignedChrOrder, barcode = $assemblyCtgLength, note = '$assemblySeqList' WHERE id = $assemblyAllCtgList[0]");
					$assignedChrOrder++;
				}
			}
		}
		my $assemblyDetails = decode_json $assembly[8];
		$assemblyDetails->{'log'} = '' if (!exists $assemblyDetails->{'log'});
		$assemblyDetails->{'log'} .= "\n" if ($assemblyDetails->{'log'});
		#add log info
		$assemblyDetails->{'log'} .= "==== ". localtime() . " ====\n"
			."New Assembly: $replace;\n"
			."With All Sequences: $withAllSequences;\n"
			."Physical Reference: $fpcOrAgpId;\n"
			."Reference Genome: $refGenomeId; alignmentBlockSize $alignmentBlockSize; alignmentCoverage $alignmentCoverage; Assign Chr: $assignChr; Orient Contigs: $orientContigs;\n"
			."End-to-End: $endToEnd, Identity:$identityEndToEnd, Overlap:$minOverlapEndToEnd;\n"
			."Filter Redundancy: Seq - $redundancyFilterSeq ($redundancySeq Hidden), Overlap - $redundancyFilterOverlap;\n"
			."Auto-Orient Seqs: $orientSeqs;\n"
			."Renumber Contigs: $renumber."
			;
		my $json = JSON::XS->new->allow_nonref;
		$assemblyDetails = $json->encode($assemblyDetails);
		my $updateAssemblyToWork=$dbh->prepare("UPDATE matrix SET barcode = '1', note = ? WHERE id = ?");
		$updateAssemblyToWork->execute($assemblyDetails,$assemblyId);
		if($emailNotification)
		{
			#email to user after alignment finishes.
			open(MAIL,"|/usr/sbin/sendmail -t -oi");
			print MAIL "To: $userEmail\n";
			print MAIL "From: $author\n";
			print MAIL "Subject: Assembly Successfully Completed\n\n";
			print MAIL <<eof;
Dear $userFullName ($userName),

Your assembly '$assembly[2]' job has successfully completed.

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
	parent.errorPop("Please give an assembly id!");
</script>	
END
}
