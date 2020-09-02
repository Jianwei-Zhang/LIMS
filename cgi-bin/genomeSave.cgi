#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use JSON::XS;
use Bio::SeqIO;
use DBI;
use File::Basename;
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

my $asbProjectId = param ('asbProjectId') || '';
my $genomeId = param ('genomeId') || '';
my $genomeName = param ('name') || '';
my $replace = param ('replace') || '0';
my $forAssembly = param ('forAssembly') || '0';
my $asReference = param ('asReference') || '0';
my $libraryId = param ('libraryId') || '0';

my $assignChr = param ('assignChr') || '0';
my $split = param ('split') || '0';
my $genomeDescription = param('description') || '';
my $genomeFile = upload ('genomeFile');
my $genomeFilePath = param ('genomeFilePath') || '';
my $agpFile = upload ('agpFile');
my $agpFilename = param("agpFile");

my $safeFilenameCharacters = "a-zA-Z0-9_.-";
my ( $filename, $filepath, $extension ) = fileparse ( $agpFilename, '..*' );
$agpFilename = $filename . $extension;
$agpFilename =~ tr/ /_/;
$agpFilename =~ s/[^$safeFilenameCharacters]//g;

my $agpObjectComponent = param ('agpObjectComponent') || '0';
my $genomeInfile = "$commoncfg->{TMPDIR}/$$.genome";
my $agpInfile = "$commoncfg->{TMPDIR}/$$.agp";
my $json = JSON::XS->new->allow_nonref;

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

if($genomeName)
{
	if($genomeId)
	{
		my $checkGenomeName = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'genome' AND name LIKE ? AND id != ?");
		$checkGenomeName->execute($genomeName, $genomeId);
		if($checkGenomeName->rows < 1)
		{
			my $checkCreator = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
			$checkCreator->execute($genomeId);
			my @checkCreator=$checkCreator->fetchrow_array();
			if($checkCreator[9] eq $userName || exists $userPermission->{$userId}->{'genome'})
			{
				my $updateGenome=$dbh->prepare("UPDATE matrix SET name = ?, x = ?, y = ?, z = ?, note = ? WHERE id = ?");
				$updateGenome->execute($genomeName,$forAssembly,$asReference,$libraryId,$genomeDescription,$genomeId);

				if($agpFile) # since agp file is relatively small, we upload it first
				{
					open (FILE, ">$agpInfile");
					while (read ($agpFile, my $agpBuffer, 1024)) {
						print FILE $agpBuffer;
					}
					close FILE;
					my $agpDetails;
					open(AGP, "<$agpInfile") or die "cannot open file $agpInfile";
					{
						local $/;
						$agpDetails = <AGP>;
					}
					close(AGP);

					my $version = 1;
					my $checkAgp = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'agp' AND x = ?");
					$checkAgp->execute($genomeId);
					if($checkAgp->rows < 1)
					{
						my $insertAgp=$dbh->prepare("INSERT INTO matrix VALUES ('', 'agp', ?, ?, ?, ?, 0, 0, ?, ?, NOW())");
						$insertAgp->execute($agpFilename,$version,$genomeId,$agpObjectComponent,$agpDetails,$userName);
						my $agpId = $dbh->{mysql_insertid};
					}
					else
					{
						while(my @checkAgp=$checkAgp->fetchrow_array())
						{
							$version = $checkAgp[3] if ($checkAgp[3] > $version);
						}
						$version++;
						my $insertAgp=$dbh->prepare("INSERT INTO matrix VALUES ('', 'agp', ?, ?, ?, ?, 0, 0, ?, ?, NOW())");
						$insertAgp->execute($agpFilename,$version,$genomeId,$agpObjectComponent,$agpDetails,$userName);
						my $agpId = $dbh->{mysql_insertid};
					}
					unlink ($agpInfile);
				}
				if($genomeFile || $genomeFilePath)
				{
					my $pid = fork();
					if ($pid) {
						print <<END;
<script>
parent.closeDialog();
parent.informationPop("It's loading! This processing might take a while.");
parent.refresh("general");
</script>	
END
					}
					elsif($pid == 0){
						my $updateGenomeToRunning = $dbh->do("UPDATE matrix SET o = 0, barcode = -1 WHERE id = $genomeId");			

						if($genomeFilePath)
						{
							#$genomeInfile = $genomeFilePath;
							`cp $genomeFilePath $genomeInfile`;
						}
						else
						{
							open (FILE, ">$genomeInfile");
							while (read ($genomeFile, my $Buffer, 1024)) {
								print FILE $Buffer;
							}
							close FILE;
						}		
						`perl -p -i -e 's/\r/\n/g' $genomeInfile`; #convert CR line terminators (MAC OS) into LF line terminators (Unix)

						if ($replace) #delete old sequences
						{
							my $genomeAsReference = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assembly' AND y = $genomeId");
							$genomeAsReference->execute();
							if($genomeAsReference->rows > 0)
							{
								print <<END;
<script>
parent.closeDialog();
parent.errorPop("You can NOT replace this genome since it has already been used as a reference in an assembly!");
</script>	
END

								unlink ($genomeInfile);
								exit;
							}

							my $genomeForAssembly = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assembly' AND x = $genomeId");
							$genomeForAssembly->execute();
							if ($genomeForAssembly->rows > 0)
							{
								print <<END;
<script>
parent.closeDialog();
parent.errorPop("You can NOT replace this genome since it has already been used in an assembly!");
</script>	
END
								unlink ($genomeInfile);

								exit;
							}

							my $genomeSequence = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND (o = 99 OR o = 97) AND x = $genomeId");
							$genomeSequence->execute();
							while(my @genomeSequence = $genomeSequence->fetchrow_array())
							{
								my $deleteAlignment=$dbh->do("DELETE FROM alignment WHERE query = $genomeSequence[0] OR subject = $genomeSequence[0]");
							}
							my $deleteGenomeSequence = $dbh->do("DELETE FROM matrix WHERE container LIKE 'sequence' AND (o = 99 OR o = 97) AND x = $genomeId");
						}

						close (STDOUT);

						#loading sequence
						my @allSequenceId;
						my $in = Bio::SeqIO->new(-file => $genomeInfile,
												-format => 'Fasta');
						while ( my $seq = $in->next_seq() )
						{
							my $seqDetails;
							$seqDetails->{'id'} = $seq->id;
							$seqDetails->{'description'} = $seq->desc || '';
							$seqDetails->{'sequence'} = $seq->seq;
							$seqDetails->{'sequence'} =~ tr/a-zA-Z/N/c; #replace nonword characters.
#							$seqDetails->{'sequence'} =~ s/^N+|N+$//g; #replace Ns at both ends.	
							$seqDetails->{'gapList'} = '';
							my $chrNumber = 0;
							if($assignChr)
							{
								$chrNumber = ($seqDetails->{'id'} =~ /(\d+)/) ? $1 : 0;
							}
							my $seqLength = $seq->length();
							my $seqend = 0;
							my @piece = ();
							foreach (split(/([N|n]{20,})/,$seqDetails->{'sequence'})) #at least 20 Ns to be a gap
							{
								my $seqstart=$seqend+1;
								$seqend=$seqend + length($_);
								if($_ =~ /^[N|n]+$/)
								{
									$seqDetails->{'gapList'} .= ($seqDetails->{'gapList'} ne '') ? ",$seqstart-$seqend" : "$seqstart-$seqend" ;
								}
								else
								{
									my $pieceDetails;
									$pieceDetails->{'id'} = $seq->id;
									$pieceDetails->{'description'} = $seq->desc || '';
									$pieceDetails->{'sequence'} = $_;
									$pieceDetails->{'gapList'} = '';
									$pieceDetails->{'coordinates'} = "$seqstart-$seqend";
									my $pieceDetailsEncoded = $json->encode($pieceDetails);
									push @piece, $pieceDetailsEncoded;
								}
							}
							my $seqDetailsEncoded = $json->encode($seqDetails);
							my $insertSequence=$dbh->prepare("INSERT INTO matrix VALUES ('', 'sequence', ?, 99, ?, ?, ?, 0, ?, ?, NOW())");
							$insertSequence->execute($seq->id,$genomeId,$seqLength,$chrNumber,$seqDetailsEncoded,$userName);
							my $parentSequenceId = $dbh->{mysql_insertid};
							push @allSequenceId,$parentSequenceId;

							if($split && @piece > 1)
							{
								my $pieceNumber = 1;
								for(@piece)
								{
									my $pieceDetails = decode_json $_;
									$pieceDetails->{'id'} = '' unless (exists $pieceDetails->{'id'});
									$pieceDetails->{'sequence'} = '' unless (exists $pieceDetails->{'sequence'});
									$seqLength = length ($pieceDetails->{'sequence'});
									my $insertChildSequence=$dbh->prepare("INSERT INTO matrix VALUES ('', 'sequence', ?, 97, ?, ?, ?, ?, ?, ?, NOW())");
									$insertChildSequence->execute($pieceDetails->{'id'}.'-'.$pieceNumber,$genomeId,$seqLength,$pieceNumber,$parentSequenceId,$_,$userName);
									my $childSequenceId = $dbh->{mysql_insertid};
									push @allSequenceId,$childSequenceId;
									$pieceNumber++;
								}
							}					
						}
						unlink ($genomeInfile);
						unlink ($agpInfile);
						foreach (@allSequenceId)
						{
							my $getSequences = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
							$getSequences->execute($_);
							while(my @getSequences = $getSequences->fetchrow_array())
							{
								my $sequenceDetails = decode_json $getSequences[8];
								$sequenceDetails->{'sequence'} = '' unless (exists $sequenceDetails->{'sequence'});
								my $seqDir;
								for (my $position = 0; $position < length($getSequences[4]); $position += 2)
								{
									$seqDir .= "/g". substr($getSequences[4],$position,2);
									until (-e "$commoncfg->{DATADIR}/sequences$seqDir")
									{
										mkdir "$commoncfg->{DATADIR}/sequences$seqDir";
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
											mkdir "$commoncfg->{DATADIR}/sequences$seqDir$seqFile";
										}
									}
								}

								until (-e "$commoncfg->{DATADIR}/sequences$seqDir$seqFile.fa")
								{
									open (SEQ,">$commoncfg->{DATADIR}/sequences$seqDir$seqFile.fa") or die "can't open file: $commoncfg->{DATADIR}/sequences$seqDir$seqFile.fa";
									print SEQ ">$getSequences[0]\n";
									print SEQ &multiLineSeq($sequenceDetails->{'sequence'},80);
									close(SEQ);
								}

								$sequenceDetails->{'sequence'} = "sequences$seqDir$seqFile.fa";
								my $seqDetailsEncoded = $json->encode($sequenceDetails);
								my $updateSequence = $dbh->prepare("UPDATE matrix SET note = ? WHERE id = ?");
								$updateSequence->execute($seqDetailsEncoded,$getSequences[0]);
							}
						}
						
						my $seqNumber = 0;
						my $countGenomeSequence = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 99 AND x = $genomeId");
						$countGenomeSequence->execute();
						$seqNumber = $countGenomeSequence->rows;
						my $updateGenomeToLoaded = $dbh->do("UPDATE matrix SET o = $seqNumber, barcode = 1, creationDate = NOW() WHERE id = $genomeId");			
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
parent.closeDialog();
parent.refresh("general");
</script>	
END
				}
			}
			else
			{
				print <<END;
<script>
	parent.errorPop("Not a valid user!");
</script>	
END
			}
		}
		else
		{
			print <<END;
<script>
	parent.errorPop("The given genome name is existing!");
</script>	
END
		}
	}
	else
	{
		if($genomeFile || $genomeFilePath)
		{
			my $checkGenomeName = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'genome' AND name LIKE ?");
			$checkGenomeName->execute($genomeName);
			if($checkGenomeName->rows < 1)
			{
				my $insertGenome=$dbh->prepare("INSERT INTO matrix VALUES ('', 'genome', ?, 0, ?, ?, ?, 0, ?, ?, NOW())");
				$insertGenome->execute($genomeName,$forAssembly,$asReference,$libraryId,$genomeDescription,$userName);
				$genomeId = $dbh->{mysql_insertid};
				my $insertLink = $dbh->do("INSERT INTO link VALUES ($asbProjectId,$genomeId,'asbProject')") if ($asbProjectId);
				my $pid = fork();
				if ($pid) {
					print <<END;
<script>
parent.closeDialog();
parent.refresh("general");
</script>	
END
				}
				elsif($pid == 0){
					close (STDOUT);
					my $updateGenomeToRunning = $dbh->do("UPDATE matrix SET barcode = -1 WHERE id = $genomeId");			
					if($agpFile) # since agp file is relatively small, we upload it first
					{
						open (FILE, ">$agpInfile");
						while (read ($agpFile, my $agpBuffer, 1024)) {
							print FILE $agpBuffer;
						}
						close FILE;
						my $agpDetails;
						open(AGP, "<$agpInfile") or die "cannot open file $agpInfile";
						{
							local $/;
							$agpDetails = <AGP>;
						}
						close(AGP);
						my $version = 1;
						my $insertAgp=$dbh->prepare("INSERT INTO matrix VALUES ('', 'agp', ?, ?, ?, ?, 0, 0, ?, ?, NOW())");
						$insertAgp->execute($agpFilename,$version,$genomeId,$agpObjectComponent,$agpDetails,$userName);
						my $agpId = $dbh->{mysql_insertid};
					}
					if($genomeFilePath)
					{
						#$genomeInfile = $genomeFilePath;
						`cp $genomeFilePath $genomeInfile`;
					}
					else
					{
						open (FILE, ">$genomeInfile");
						while (read ($genomeFile, my $Buffer, 1024)) {
							print FILE $Buffer;
						}
						close FILE;
					}		
					`perl -p -i -e 's/\r/\n/g' $genomeInfile`; #convert CR line terminators (MAC OS) into LF line terminators (Unix)
					#loading sequence
					my @allSequenceId;
					my $in = Bio::SeqIO->new(-file => $genomeInfile,
											-format => 'Fasta');
					while ( my $seq = $in->next_seq() )
					{
						my $seqDetails;
						$seqDetails->{'id'} = $seq->id;
						$seqDetails->{'description'} = $seq->desc || '';
						$seqDetails->{'sequence'} = $seq->seq;
						$seqDetails->{'sequence'} =~ tr/a-zA-Z/N/c; #replace nonword characters.
#						$seqDetails->{'sequence'} =~ s/^N+|N+$//g; #replace Ns at both ends.	
						$seqDetails->{'gapList'} = '';
						my $chrNumber = 0;
						if($assignChr)
						{
							$chrNumber = ($seqDetails->{'id'} =~ /(\d+)/) ? $1 : 0;
						}
						my $seqLength = $seq->length();
						my $seqend = 0;
						my @piece = ();
						foreach (split(/([N|n]{20,})/,$seqDetails->{'sequence'})) #at least 20 Ns to be a gap
						{
							my $seqstart=$seqend+1;
							$seqend=$seqend + length($_);
							if($_ =~ /^[N|n]+$/)
							{
								$seqDetails->{'gapList'} .= ($seqDetails->{'gapList'} ne '') ? ",$seqstart-$seqend" : "$seqstart-$seqend" ;
							}
							else
							{
								my $pieceDetails;
								$pieceDetails->{'id'} = $seq->id;
								$pieceDetails->{'description'} = $seq->desc || '';
								$pieceDetails->{'sequence'} = $_;
								$pieceDetails->{'gapList'} = '';
								$pieceDetails->{'coordinates'} = "$seqstart-$seqend";
								my $pieceDetailsEncoded = $json->encode($pieceDetails);
								push @piece, $pieceDetailsEncoded;
							}
						}
						$seqDetails = $json->encode($seqDetails);
						my $insertSequence=$dbh->prepare("INSERT INTO matrix VALUES ('', 'sequence', ?, 99, ?, ?, ?, 0, ?, ?, NOW())");
						$insertSequence->execute($seq->id,$genomeId,$seqLength,$chrNumber,$seqDetails,$userName);
						my $parentSequenceId = $dbh->{mysql_insertid};
						push @allSequenceId,$parentSequenceId;

						if($split && @piece > 1)
						{
							my $pieceNumber = 1;
							for(@piece)
							{
								my $pieceDetails = decode_json $_;
								$pieceDetails->{'id'} = '' unless (exists $pieceDetails->{'id'});
								$pieceDetails->{'sequence'} = '' unless (exists $pieceDetails->{'sequence'});
								$seqLength = length ($pieceDetails->{'sequence'});
								my $insertChildSequence=$dbh->prepare("INSERT INTO matrix VALUES ('', 'sequence', ?, 97, ?, ?, ?, ?, ?, ?, NOW())");
								$insertChildSequence->execute($pieceDetails->{'id'}.'-'.$pieceNumber,$genomeId,$seqLength,$pieceNumber,$parentSequenceId,$_,$userName);
								my $childSequenceId = $dbh->{mysql_insertid};
								push @allSequenceId,$childSequenceId;
								$pieceNumber++;
							}
						}
					}
					unlink ($genomeInfile);
					unlink ($agpInfile);

					foreach (@allSequenceId)
					{
						my $getSequences = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
						$getSequences->execute($_);
						while(my @getSequences = $getSequences->fetchrow_array())
						{
							my $sequenceDetails = decode_json $getSequences[8];
							$sequenceDetails->{'sequence'} = '' unless (exists $sequenceDetails->{'sequence'});
							my $seqDir;
							for (my $position = 0; $position < length($getSequences[4]); $position += 2)
							{
								$seqDir .= "/g". substr($getSequences[4],$position,2);
								until (-e "$commoncfg->{DATADIR}/sequences$seqDir")
								{
									mkdir "$commoncfg->{DATADIR}/sequences$seqDir";
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
										mkdir "$commoncfg->{DATADIR}/sequences$seqDir$seqFile";
									}
								}
							}

							until (-e "$commoncfg->{DATADIR}/sequences$seqDir$seqFile.fa")
							{
								open (SEQ,">$commoncfg->{DATADIR}/sequences$seqDir$seqFile.fa") or die "can't open file: $commoncfg->{DATADIR}/sequences$seqDir$seqFile.fa";
								print SEQ ">$getSequences[0]\n";
								print SEQ &multiLineSeq($sequenceDetails->{'sequence'},80);
								close(SEQ);
							}

							$sequenceDetails->{'sequence'} = "sequences$seqDir$seqFile.fa";
							my $seqDetailsEncoded = $json->encode($sequenceDetails);
							my $updateSequence = $dbh->prepare("UPDATE matrix SET note = ? WHERE id = ?");
							$updateSequence->execute($seqDetailsEncoded,$getSequences[0]);
						}
					}

					my $seqNumber = 0;
					my $countGenomeSequence = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 99 AND x = $genomeId");
					$countGenomeSequence->execute();
					$seqNumber = $countGenomeSequence->rows;
					my $updateGenomeToLoaded = $dbh->do("UPDATE matrix SET o = $seqNumber, barcode = 1, creationDate = NOW() WHERE id = $genomeId");			
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
	parent.errorPop("The given genome name is existing!");
</script>	
END
				exit;
			}
		}
		else
		{
			print <<END;
<script>
	parent.errorPop("No Genome file found!");
</script>
END
			exit;
		}
	}
}
else
{
	print <<END;
<script>
	parent.errorPop("Please give a genome name!");
</script>	
END
}

