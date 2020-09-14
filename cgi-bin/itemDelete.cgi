#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use JSON::XS;
use DBI;
use lib "lib/";
use lib "lib/pangu";
use pangu;
use user;
use config;
use userCookie;

my $userCookie = new userCookie;
my $userId = (cookie('cid')) ? $userCookie->checkCookie(cookie('cid')) : 0;
exit if (!$userId);

my $user = new user;
my $userDetail = $user->getAllFieldsWithUserId($userId);
my $userName = $userDetail->{"userName"};

my $commoncfg = readConfig("main.conf");
my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});

my %smrtrunStatus = (
	0=>'Initialized',
	1=>'Started',
	2=>'Completed'
	);
my @items = param ('items');
my $option =  param ('option') || '';
my $config = new config;
my $userPermission;
my $userInGroup=$dbh->prepare("SELECT * FROM link WHERE type LIKE 'group' AND child = ?");
$userInGroup->execute($userId);
while (my @userInGroup = $userInGroup->fetchrow_array())
{
	my $group = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$group->execute($userInGroup[0]);
	my @group=$group->fetchrow_array();
	my $groupDetails = decode_json $group[8];
	$groupDetails->{'permission'} = '' unless (exists $groupDetails->{'permission'});
	for(split (",",$groupDetails->{'permission'}))
	{
		my $permission = $config->getFieldName($_);
		$userPermission->{$userId}->{$permission} = 1;		
	}
}
print header;

my $libraryId = 0;
for(@items)
{
	my $itemId = $_;
	if ($itemId eq 0)
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
			if ($option eq 'delTempFiles')
			{
				#delete temp files
				opendir(DIR, $commoncfg->{TMPDIR}) or die "can't opendir $commoncfg->{TMPDIR}: $!";
				my @files = readdir(DIR);
				closedir DIR;
				foreach my $file (sort @files)
				{
					next if ($file =~ /^\./);
					unlink "$commoncfg->{TMPDIR}/$file" if (-f "$commoncfg->{TMPDIR}/$file"); #delete files
				}
			}
			else
			{
				#connect to the mysql server
				my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});
				
				#background check for orphan stuffs to be deleted.
				# orphan alignments
				if(-e "$commoncfg->{DATADIR}/alignments/seqToSeq")
				{
					my @dirs;
					push @dirs, "$commoncfg->{DATADIR}/alignments/seqToSeq";

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
									my $queryId = $1;
									my $subjectId = $2;
									my $toBeDeleted = 0;
									my $querySequence = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
									$querySequence->execute($queryId);
									unless ($querySequence->rows > 0)
									{
										$toBeDeleted++;
									}
									my $subjectSequence = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
									$subjectSequence->execute($subjectId);
									unless ($subjectSequence->rows > 0)
									{
										$toBeDeleted++;
									}
									if($toBeDeleted)
									{
										unlink ("$dir/$file"); #delete old alignments
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
											 unlink ("$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$subjectId-$queryId.tbl"); #delete old alignments
										}
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

				#
				#delete comment
				my $commentIds;
				my $commentParents = $dbh->prepare("SELECT o FROM matrix WHERE container LIKE 'comment' GROUP BY o");
				$commentParents->execute();
				while (my @commentParents =  $commentParents->fetchrow_array())
				{
					$commentIds->{$commentParents[0]} = 1;
				}
				for (keys %{$commentIds})
				{
					my $commentParent = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
					$commentParent->execute($_);
					unless ($commentParent->rows > 0)
					{
						my $deleteComment=$dbh->do("DELETE FROM matrix WHERE container LIKE 'comment' AND o = $_");
					}
				}
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
		my $item = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$item->execute($itemId);
		my @item=$item->fetchrow_array();
		if($item[9] eq $userName || exists $userPermission->{$userId}->{'item'})
		{
			if($item[1] eq 'room')
			{
				my $freezerInRoom = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'freezer' AND o = $itemId");
				$freezerInRoom->execute();
				if ($freezerInRoom->rows > 0)
				{
					print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("Please delete freezers in the room first!");
	</script>	
END
					exit;
				}
				my $deleteRoom=$dbh->do("DELETE FROM matrix WHERE id = $itemId");
				print <<END;
	<script>
		parent.closeDialog();
		parent.refresh("menu");
	</script>	
END
			}
			elsif($item[1] eq 'freezer')
			{
				my $boxInFreezer = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'box' AND o = $itemId");
				$boxInFreezer->execute();
				if ($boxInFreezer->rows > 0)
				{
					print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("Please delete boxes in the freezer first!");
	</script>	
END
					exit;
				}
				my $deleteFreezer=$dbh->do("DELETE FROM matrix WHERE id = $itemId");
				print <<END;
	<script>
		parent.closeDialog();
		parent.refresh("menu");
	</script>	
END
			}
			elsif($item[1] eq 'box')
			{
				my $itemInBox = $dbh->prepare("SELECT * FROM link WHERE type LIKE 'box' AND parent = $itemId");
				$itemInBox->execute();
				if ($itemInBox->rows > 0)
				{
					print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("Items in the box will become orphans!");
	</script>	
END
				}
				my $deleteBox=$dbh->do("DELETE FROM matrix WHERE id = $itemId");
				my $deleteLinkA=$dbh->do("DELETE FROM link WHERE parent = $itemId AND type LIKE 'box'");
				print <<END;
	<script>
		parent.closeDialog();
		parent.refresh("menu");
	</script>	
END
			}
			elsif($item[1] eq "project")
			{
				my $serviceInProject = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'service' AND z = $itemId");
				$serviceInProject->execute();
				if($serviceInProject->rows > 0)
				{
					print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("Please delete related services first!");
	</script>	
END
					exit;
				}

				my $libraryInProject = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'library' AND z = $itemId");
				$libraryInProject->execute();
				if($libraryInProject->rows > 0)
				{
					print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("Please delete related libraries first!");
	</script>	
END
					exit;
				}
				my $deleteProject=$dbh->do("DELETE FROM matrix WHERE id = $itemId");
				print <<END;
	<script>
		parent.closeDialog();
		parent.refresh("menu");
	</script>	
END
			}
			elsif($item[1] eq 'asbProject')
			{
				my $childInAsbProject = $dbh->prepare("SELECT * FROM link WHERE parent = $itemId AND type LIKE 'asbProject'");
				$childInAsbProject->execute();
				if($childInAsbProject->rows > 0)
				{
					print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("Please unlink related libraries/genomes first!");
	</script>	
END
					exit;
				}

				my $deleteAsbProject=$dbh->do("DELETE FROM matrix WHERE id = $itemId");
				my $deleteLink=$dbh->do("DELETE FROM link WHERE parent = $itemId AND type LIKE 'asbProject'");
				print <<END;
	<script>
		parent.closeDialog();
		parent.refresh("menu");
	</script>	
END
			}
			elsif($item[1] eq 'assembly')
			{
				if ($option eq 'unplaced')
				{
					my $assemblyCtg = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblyCtg' AND o = $itemId AND x = 0");
					$assemblyCtg->execute();
					while (my @assemblyCtg=$assemblyCtg->fetchrow_array())
					{
						foreach (split ",", $assemblyCtg[8])
						{
							$_ =~ s/[^a-zA-Z0-9]//g;
							my $deleteAssemblySeq=$dbh->do("DELETE FROM matrix WHERE id = $_");
						}
						my $deleteAssemblyCtg=$dbh->do("DELETE FROM matrix WHERE id = $assemblyCtg[0]");
						my $deleteComment=$dbh->do("DELETE FROM matrix WHERE container LIKE 'comment' AND o = $assemblyCtg[0]");
					}
				}
				elsif ($option eq 'contamination')
				{
					my $assemblyCtg = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblyCtg' AND o = $itemId AND x = 100");
					$assemblyCtg->execute();
					while (my @assemblyCtg=$assemblyCtg->fetchrow_array())
					{
						foreach (split ",", $assemblyCtg[8])
						{
							$_ =~ s/[^a-zA-Z0-9]//g;
							my $deleteAssemblySeq=$dbh->do("DELETE FROM matrix WHERE id = $_");
						}
						my $deleteAssemblyCtg=$dbh->do("DELETE FROM matrix WHERE id = $assemblyCtg[0]");
						my $deleteComment=$dbh->do("DELETE FROM matrix WHERE container LIKE 'comment' AND o = $assemblyCtg[0]");
					}
				}
				else
				{
					my $assemblyCtg = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblyCtg' AND o = $itemId");
					$assemblyCtg->execute();
					while (my @assemblyCtg=$assemblyCtg->fetchrow_array())
					{
						foreach (split ",", $assemblyCtg[8])
						{
							$_ =~ s/[^a-zA-Z0-9]//g;
							my $deleteAssemblySeq=$dbh->do("DELETE FROM matrix WHERE id = $_");
						}
						my $deleteAssemblyCtg=$dbh->do("DELETE FROM matrix WHERE id = $assemblyCtg[0]");
						my $deleteComment=$dbh->do("DELETE FROM matrix WHERE container LIKE 'comment' AND o = $assemblyCtg[0]");
					}
					my $deleteAsbGenomeLink = $dbh->do("DELETE FROM link WHERE parent = $itemId AND type LIKE 'asbGenome'");
					my $deleteAssembly=$dbh->do("DELETE FROM matrix WHERE id = $itemId");
				}
				print <<END;
	<script>
		parent.closeDialog();
		parent.refresh("menu");
	</script>	
END
			}
			elsif($item[1] eq 'assemblyCtg')
			{
				foreach (split ",", $item[8])
				{
					$_ =~ s/[^a-zA-Z0-9]//g;
					my $deleteAssemblySeq=$dbh->do("DELETE FROM matrix WHERE id = $_");
				}
				my $deleteAssemblyCtg=$dbh->do("DELETE FROM matrix WHERE id = $itemId");
				my $deleteComment=$dbh->do("DELETE FROM matrix WHERE container LIKE 'comment' AND o = $itemId");

				print <<END;
	<script>
		parent.closeDialog();
		parent.closeViewer();
		parent.refresh("menu");
	</script>	
END
			}
			elsif($item[1] eq 'assemblySeq')
			{
				my $scrollLeft = param ('scrollLeft') || '0';
				my $assemblyCtgOfSeq = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblyCtg' AND MATCH (note) AGAINST (?)");
				$assemblyCtgOfSeq->execute($itemId);
				my @assemblyCtgOfSeq = $assemblyCtgOfSeq->fetchrow_array();
				my @seqList = split ",", $assemblyCtgOfSeq[8];
				my $number=@seqList;
				my @newSeqList;
				if ($number > 1)
				{
					for (@seqList)
					{
						if (/($itemId)/)
						{
							next;
						}
						else
						{
							push @newSeqList, $_;
						}
					}
					my $newAsbSeqList = join ",", @newSeqList;
					my $updateAssemblyCtg=$dbh->do("UPDATE matrix SET note = '$newAsbSeqList' WHERE id = $assemblyCtgOfSeq[0]");
					my $deleteAssemblySeq=$dbh->do("DELETE FROM matrix WHERE id = $itemId");
					print <<END;
	<script>
	parent.closeViewer();
	parent.openViewer("assemblyCtgView.cgi?assemblyCtgId=$assemblyCtgOfSeq[0]&scrollLeft=$scrollLeft");
	//parent.refresh("menu");
	</script>	
END
				}
				else
				{
					my $deleteAssemblyCtg=$dbh->do("DELETE FROM matrix WHERE id = $assemblyCtgOfSeq[0]");
					my $deleteAssemblySeq=$dbh->do("DELETE FROM matrix WHERE id = $itemId");
					my $deleteComment=$dbh->do("DELETE FROM matrix WHERE container LIKE 'comment' AND o = $assemblyCtgOfSeq[0]");
					print <<END;
	<script>
		parent.closeDialog();
		parent.closeViewer();
		parent.refresh("menu");
	</script>	
END
				}
			
			}
			elsif($item[1] eq 'sequence')
			{
				my $assemblySeq=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblySeq' AND y = ?");
				$assemblySeq->execute($itemId);
				my @assemblySeq = $assemblySeq->fetchrow_array();
				if ($assemblySeq->rows > 0)
				{
					print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("Please delete related assemblies using this sequence first!");
	</script>	
END
				}
				else
				{
					my $sequenceDetails = decode_json $item[8];
					unlink ("$commoncfg->{DATADIR}/$sequenceDetails->{'sequence'}");
					my $deleteSequence=$dbh->do("DELETE FROM matrix WHERE id = $itemId");

					my $queryDir;
					for (my $position = 0; $position < length($itemId); $position += 2)
					{
						$queryDir .= "/q". substr($itemId,$position,2);
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
										unlink ("$dir/$file"); #delete old alignments
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
											unlink ("$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$subjectId-$queryId.tbl"); #delete old alignments
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

					if ($item[3] eq '99' )
					{
						my $genome=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
						$genome->execute($item[4]);
						my @genome = $genome->fetchrow_array();
						$genome[3] = $genome[3]-1;
						my $updateGenome=$dbh->do("UPDATE matrix SET o = '$genome[3]' WHERE id = $item[4]");
					}
					print <<END;
	<script>
		parent.closeDialog();
		parent.closeViewer();
		parent.informationPop("Successfully deleted!");
		parent.refresh("menu");
	</script>	
END
				}
			
			}
			elsif($item[1] eq 'genome')
			{
				my $genomeAsReference = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assembly' AND y = $itemId");
				$genomeAsReference->execute();
				if($genomeAsReference->rows > 0)
				{
					print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("Please delete related assemblies using this genome as a reference first!");
	</script>	
END
					exit;
				}

				my $genomeForAssembly = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assembly' AND x = $itemId");
				$genomeForAssembly->execute();
				if ($genomeForAssembly->rows > 0)
				{
					print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("Please delete related assemblies first!");
	</script>	
END
					exit;
				}

				my $asbGenomeForAssembly = $dbh->prepare("SELECT * FROM link WHERE child = $itemId AND type LIKE 'asbGenome'");
				$asbGenomeForAssembly->execute();
				if($asbGenomeForAssembly->rows > 0)
				{
					print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("Please unlink related assemblies first!");
	</script>	
END
					exit;
				}

				my $genomeForAsbProject = $dbh->prepare("SELECT * FROM link WHERE child = $itemId AND type LIKE 'asbProject");
				$genomeForAsbProject->execute();
				if($genomeForAsbProject->rows > 0)
				{
					print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("Please unlink from all assembly projects first!");
	</script>	
END
					exit;
				}

				my $genomeSequence = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND (o = 99 OR o = 97) AND x = $itemId");
				$genomeSequence->execute();
				while(my @genomeSequence = $genomeSequence->fetchrow_array())
				{
					my $sequenceDetails = decode_json $genomeSequence[8];
					unlink ("$commoncfg->{DATADIR}/$sequenceDetails->{'sequence'}");
					my $queryDir;
					for (my $position = 0; $position < length($genomeSequence[0]); $position += 2)
					{
						$queryDir .= "/q". substr($genomeSequence[0],$position,2);
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
										unlink ("$dir/$file"); #delete old alignments
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
											unlink ("$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirSwitched$subjectDirSwitched/$subjectId-$queryId.tbl"); #delete old alignments
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
				}
				my $deleteGenomeSequence = $dbh->do("DELETE FROM matrix WHERE container LIKE 'sequence' AND (o = 99 OR o = 97) AND x = $itemId"); # both sequence and piece
				my $deleteAgp=$dbh->do("DELETE FROM matrix WHERE container LIKE 'agp' AND x = $itemId");
				my $deleteGenome=$dbh->do("DELETE FROM matrix WHERE id = $itemId");
				print <<END;
	<script>
		parent.closeDialog();
		parent.refresh("general");
	</script>	
END
			}
			elsif($item[1] eq 'dart')
			{
		
				my $deleteDartSNP = $dbh->do("DELETE FROM matrix WHERE container LIKE 'dartSNP' AND o = $itemId");
				my $deleteDartGenotype = $dbh->do("DELETE FROM matrix WHERE container LIKE 'dartGenotype' AND o = $itemId");
				my $deleteDart=$dbh->do("DELETE FROM matrix WHERE id = $itemId");
				print <<END;
	<script>
		parent.closeDialog();
		parent.refresh("general");
	</script>	
END
			}
			elsif($item[1] eq 'dataset')
			{
		
				my $deleteRecord = $dbh->do("DELETE FROM matrix WHERE container LIKE 'record' AND z = $itemId");
				my $deleteDart=$dbh->do("DELETE FROM matrix WHERE id = $itemId");
				print <<END;
	<script>
		parent.closeDialog();
		parent.refresh("general");
	</script>	
END
			}
			elsif($item[1] eq 'record')
			{
				my $deleteRecord=$dbh->do("DELETE FROM matrix WHERE id = $itemId");
				print <<END;
	<script>
		parent.closeDialog();
		parent.refresh("menu");
	</script>	
END
			}
			elsif($item[1] eq 'fpc')
			{
		
				my $fpcAsReference = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assembly' AND z = $itemId");
				$fpcAsReference->execute();
				if($fpcAsReference->rows > 0)
				{
					print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("Please delete related assemblies using this FPC as a physical reference first!");
	</script>	
END
					exit;
				}

				my $deleteFpcClone=$dbh->do("DELETE FROM matrix WHERE container LIKE 'fpcClone' AND o = $itemId");
				my $deleteFpcCtg=$dbh->do("DELETE FROM matrix WHERE container LIKE 'fpcCtg' AND o = $itemId");
				my $deleteFpc=$dbh->do("DELETE FROM matrix WHERE id = $itemId");
				print <<END;
	<script>
		parent.closeDialog();
		parent.refresh("menu");
	</script>	
END
			}
			elsif($item[1] eq 'agp')
			{
		
				my $agpAsReference = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assembly' AND z = $itemId");
				$agpAsReference->execute();
				if($agpAsReference->rows > 0)
				{
					print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("Please delete related assemblies using this AGP as a physical reference first!");
	</script>	
END
					exit;
				}
				my $deleteAgp=$dbh->do("DELETE FROM matrix WHERE id = $itemId");
				print <<END;
	<script>
		parent.closeDialog();
		parent.refresh("general");
	</script>	
END
			}
			elsif($item[1] eq 'plate')
			{
				if ($item[3] > 0)
				{
					unless($item[9] eq $userName)
					{
						print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("You can NOT delete a label-generated plate unless your created the plate!");
	</script>	
END
						exit;
					}
				}

				my $deletePlate=$dbh->do("DELETE FROM matrix WHERE id = $itemId");
				my $deleteLink=$dbh->do("DELETE FROM link WHERE child = $itemId AND type LIKE 'box'");
				my $plateId;
				my $plateInLibrary=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'plate' AND z = $item[6] ORDER BY o");
				$plateInLibrary->execute();
				while (my @plateInLibrary = $plateInLibrary->fetchrow_array())
				{
					$plateId->{$plateInLibrary[2]}->{$plateInLibrary[4]} = $plateInLibrary[0];
				}
				my $nofPlates = keys %{$plateId};
				my $updateLibrary=$dbh->do("UPDATE matrix SET o = $nofPlates WHERE id = $item[6]");
				print <<END;
	<script>
		parent.closeDialog();
		parent.refresh("menu");
	</script>	
END
			}
			elsif($item[1] eq 'pool')
			{
				$libraryId = $item[4];
				my $deletePool=$dbh->do("DELETE FROM matrix WHERE id = $itemId");
				my $deleteLinkPoolJob=$dbh->do("DELETE FROM link WHERE parent = $itemId AND type LIKE 'poolJob'");
				my $deleteLinkPoolClone=$dbh->do("DELETE FROM link WHERE parent = $itemId AND type LIKE 'poolClone'");
				print <<END;
	<script>
		parent.closeDialog();
		parent.refresh("viewLibraryTabs$libraryId");
	</script>	
END
			}
			elsif($item[1] eq "library")
			{
				my $plateInLibrary = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'plate' AND z = $itemId");
				$plateInLibrary->execute();
				if ($plateInLibrary->rows > 0)
				{
					print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("Please delete related plates first!");
	</script>	
END
					exit;
				}

				my $rearrayingLibrary = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'library' AND y = $itemId");
				$rearrayingLibrary->execute();
				if ($rearrayingLibrary->rows > 0)
				{
					print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("Please delete related rearraying library first!");
	</script>	
END
					exit;
				}

				my $fpcForLibrary = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'fpc' AND barcode = $itemId");
				$fpcForLibrary->execute();
				if ($fpcForLibrary->rows > 0)
				{
					print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("Please delete related FPC first!");
	</script>	
END
					exit;
				}

				my $besForLibrary = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 98 AND x = $itemId");
				$besForLibrary->execute();
				if ($besForLibrary->rows > 0)
				{
					print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("Please delete related BES first!");
	</script>	
END
					exit;
				}

				my $poolForLibrary = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'pool' AND x = $itemId");
				$poolForLibrary->execute();
				if ($poolForLibrary->rows > 0)
				{
					print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("Please delete related pools first!");
	</script>	
END
					exit;
				}

				my $tagForLibrary = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'tag' AND x = $itemId");
				$tagForLibrary->execute();
				if ($tagForLibrary->rows > 0)
				{
					print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("Please delete related tags first!");
	</script>	
END
					exit;
				}

				my $assemblyForLibrary = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assembly' AND x = $itemId");
				$assemblyForLibrary->execute();
				if ($assemblyForLibrary->rows > 0)
				{
					print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("Please delete related assemblies first!");
	</script>	
END
					exit;
				}

				my $deleteLibrary=$dbh->do("DELETE FROM matrix WHERE id = $itemId");					
				print <<END;
	<script>
		parent.closeDialog();
		parent.refresh("menu");
	</script>	
END
			}
			elsif($item[1] eq "service")
			{
				my $sampleInService = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sample' AND z = $itemId");
				$sampleInService->execute();
				if ($sampleInService->rows > 0)
				{
					print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("Please delete related samples first!");
	</script>	
END
					exit;
				}
				my $deleteService=$dbh->do("DELETE FROM matrix WHERE id = $itemId");
				print <<END;
	<script>
		parent.closeDialog();
		parent.refresh("menu");
		parent.refresh("general");
	</script>	
END
			}
			elsif($item[1] eq "sample")
			{
				my $paclibForSample = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'paclib' AND z = $itemId");
				$paclibForSample->execute();
				if ($paclibForSample->rows > 0)
				{
					print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("Please delete related Paclibs first!");
	</script>	
END
					exit;
				}
				my $deleteSample=$dbh->do("DELETE FROM matrix WHERE id = $itemId");
				my $deleteLink=$dbh->do("DELETE FROM link WHERE child = $itemId AND type LIKE 'box'");
				print <<END;
	<script>
		parent.closeDialog();
		parent.refresh("menu");
		parent.refresh("general");
	</script>	
END
			}
			elsif($item[1] eq "paclib")
			{
				my $smrtcellForPaclib = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'smrtcell' AND o = $itemId");
				$smrtcellForPaclib->execute();
				if ($smrtcellForPaclib->rows > 0)
				{
					print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("Please delete related SMRT cells first!");
	</script>	
END
					exit;
				}
				my $deletePaclib=$dbh->do("DELETE FROM matrix WHERE id = $itemId");
				my $deleteLink=$dbh->do("DELETE FROM link WHERE child = $itemId AND type LIKE 'box'");
				print <<END;
	<script>
		parent.closeDialog();
		parent.refresh("menu");
		parent.refresh("general");
	</script>	
END
			}
			elsif($item[1] eq "smrtrun")
			{
				if($item[7] > 0)
				{
					print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("You can NOT delete a $smrtrunStatus{$item[7]} run!");
	</script>	
END
					exit;
				}
				else
				{
					my $smrtwellInRun = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'smrtwell' AND z = $itemId");
					$smrtwellInRun->execute();
					if ($smrtwellInRun->rows > 0)
					{
						print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("Please delete related wells first!");
	</script>	
END
						exit;
					}
					my $deleteSmrtrun=$dbh->do("DELETE FROM matrix WHERE id = $itemId");
				}
				print <<END;
	<script>
		parent.closeDialog();
		parent.refresh("menu");
		parent.refresh("general");
	</script>	
END
			}
			elsif($item[1] eq "smrtwell")
			{

				my $smrtrun=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
				$smrtrun->execute($item[6]);
				my @smrtrun = $smrtrun->fetchrow_array();
				if($smrtrun[7] > 0)
				{
					print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("You can NOT delete the well from a $smrtrunStatus{$smrtrun[7]} run!");
	</script>	
END
					exit;
				}
				my $deleteSmrtwell=$dbh->do("DELETE FROM matrix WHERE id = $itemId");
				print <<END;
	<script>
		parent.closeDialog();
		parent.refresh("menu");
		parent.refresh("general");
	</script>	
END
			}
			elsif($item[1] eq "vector")
			{
				my $libraryWithVector = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE library AND x = $itemId");
				$libraryWithVector->execute();
				if ($libraryWithVector->rows > 0)
				{
					print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("Please delete related libraries first!");
	</script>	
END
						exit;
				}
				my $deleteVector=$dbh->do("DELETE FROM matrix WHERE id = $itemId");
				print <<END;
	<script>
		parent.closeDialog();
		parent.refresh("general");
	</script>	
END
			}
			elsif($item[1] eq 'group')
			{
				my $deleteGroup=$dbh->do("DELETE FROM matrix WHERE id = $itemId");
				my $deleteLink=$dbh->do("DELETE FROM link WHERE parent = $itemId AND type LIKE 'group'");
				print <<END;
	<script>
		parent.closeDialog();
		parent.refresh("menu");
	</script>	
END
			}
			elsif($item[1] eq 'comment')
			{
				my $deleteComment=$dbh->do("DELETE FROM matrix WHERE id = $itemId");
				print <<END;
	<script>
		parent.closeDialog();
		parent.refresh("menu");
	</script>	
END
				exit;
			}
			else
			{
				print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("Not a valid deletion!");
	</script>	
END
				exit;
			}
		}
		else
		{
			print <<END;
	<script>
		parent.closeDialog();
		parent.errorPop("Not a valid user!");
	</script>	
END
			exit;
		}
	}
}