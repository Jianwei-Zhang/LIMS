#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use JSON::XS;
use DBI;
use Bio::SeqIO;
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

my $dartId = param ('dartId') || '';
my $dartName = param ('name') || '';
my $genebankId = param ('genebankId') || '0';

my $dartDescription = param('description') || '';
my $dartFile = upload ('dartFile');
my $dartFilePath = param ('dartFilePath') || '';

my $safeFilenameCharacters = "a-zA-Z0-9_.-";
my $dartInfile = "$commoncfg->{TMPDIR}/$$.dart";
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

if($dartName)
{
	if($dartId)
	{
		my $checkDartName = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'dart' AND name LIKE ? AND id != ?");
		$checkDartName->execute($dartName, $dartId);
		if($checkDartName->rows < 1)
		{
			my $checkCreator = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
			$checkCreator->execute($dartId);
			my @checkCreator=$checkCreator->fetchrow_array();
			if($checkCreator[9] eq $userName || exists $userPermission->{$userId}->{'dart'})
			{
				my $updateDart=$dbh->prepare("UPDATE matrix SET name = ?, z = ?, note = ? WHERE id = ?");
				$updateDart->execute($dartName,$genebankId,$dartDescription,$dartId);
				if($dartFile || $dartFilePath)
				{
					my $pid = fork();
					if ($pid)
					{
						print <<END;
<script>
parent.closeDialog();
parent.informationPop("It's loading! This processing might take a while.");
parent.refresh("general");
</script>	
END
					}
					elsif($pid == 0)
					{
						close (STDOUT);
						my $updateDartToRunning = $dbh->do("UPDATE matrix SET barcode = -1 WHERE id = $dartId");			
						if($dartFilePath)
						{
							`cp $dartFilePath $dartInfile`;
						}
						else
						{
							open (FILE, ">$dartInfile");
							while (read ($dartFile, my $Buffer, 1024)) {
								print FILE $Buffer;
							}
							close FILE;
						}		
						#delete old data set
						my $deleteDartSNP = $dbh->do("DELETE FROM matrix WHERE container LIKE 'dartSNP' AND z = $dartId");
						my $deleteDartGenotype = $dbh->do("DELETE FROM matrix WHERE container LIKE 'dartGenotype' AND z = $dartId");
						`perl -p -i -e 's/\r/\n/g' $dartInfile`; #convert CR line terminators (MAC OS) into LF line terminators (Unix)
						my @headerRows = ("Order number","DArT plate barcode","client plate barcode","well row position","well column position","sample comments","genotype name");
						my $snpNumber = 0;
						my $genotype;
						my $genotypeNumber = 0;
						my $headerRow = 0;
						my $metadataCol = 0;
						my @metadata;
						open (INFILE, "$dartInfile");
						while (<INFILE>) {
							chop;
							my @dartLine =  split("\t",$_);
							if($dartLine[0] eq '*')
							{
								$metadataCol = 0;
								$genotypeNumber = 0;
								foreach (@dartLine)
								{
									if($_ eq '*')
									{
										$metadataCol++;
									}
									else
									{
										$genotype->{$genotypeNumber}->{$headerRows[$headerRow]} = $_;
										$genotypeNumber++;
									}
								
								}
								$headerRow++;
							}
							else
							{
								if($dartLine[0] =~ /ID/) # any words with 'ID', e.g. CloneID, SNP-ID, etc.
								{
									my $metadataLineCol = 0;
									$genotypeNumber = 0;
									foreach (@dartLine)
									{
										if($metadataLineCol < $metadataCol)
										{
											push @metadata, $_;
										}
										else
										{
											$genotype->{$genotypeNumber}->{'genotype name'} = $_;
											$genotypeNumber++;
										}
										$metadataLineCol++;
									}
								}
								else
								{
									my $metadataSNPCol = 0;
									my $dartSNP;
									my $dartSNPId = 0;
									my $dartCloneId = '';
									$genotypeNumber = 0;
									foreach (@dartLine)
									{
										if($metadataSNPCol < $metadataCol)
										{
											$dartCloneId = $_ if (!$dartCloneId);
											$dartSNP->{$metadata[$metadataSNPCol]} = $_;
										}
										else
										{
											if(!$dartSNPId)
											{
												my $dartSNPEncoded = $json->encode($dartSNP);
												my $insertDartSNP=$dbh->prepare("INSERT INTO matrix VALUES ('', 'dartSNP', ?, 0, 0, 0, ?, 0, ?, ?, NOW())");
												$insertDartSNP->execute($dartCloneId,$dartId,$dartSNPEncoded,$userName);
												$dartSNPId = $dbh->{mysql_insertid};
											}
											$genotype->{$genotypeNumber}->{$dartSNPId} = $_;
											$genotypeNumber++;
										}
										$metadataSNPCol++;
									}
									$snpNumber++;
								}
							}
						}
						close INFILE;
						for (sort {$a <=> $b} keys %$genotype)
						{
							my $dartGenotypeEncoded = $json->encode($genotype->{$_});
							my $insertDartGenotype=$dbh->prepare("INSERT INTO matrix VALUES ('', 'dartGenotype', ?, 0, 0, 0, ?, 0, ?, ?, NOW())");
							$insertDartGenotype->execute($genotype->{$_}->{'genotype name'},$dartId,$dartGenotypeEncoded,$userName);
						}
						unlink ($dartInfile);
						my $updateDartToLoaded = $dbh->do("UPDATE matrix SET o = $snpNumber, x = $genotypeNumber, barcode = 1, creationDate = NOW() WHERE id = $dartId");			
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
	parent.errorPop("The given dart name is existing!");
</script>	
END
		}
	}
	else
	{
		if($dartFile || $dartFilePath)
		{
			my $checkDartName = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'dart' AND name LIKE ?");
			$checkDartName->execute($dartName);
			if($checkDartName->rows < 1)
			{
				my $insertDart=$dbh->prepare("INSERT INTO matrix VALUES ('', 'dart', ?, 0, 0, 0, ?, 0, ?, ?, NOW())");
				$insertDart->execute($dartName,$genebankId,$dartDescription,$userName);
				$dartId = $dbh->{mysql_insertid};
				my $pid = fork();
				if ($pid)
				{
					print <<END;
<script>
parent.closeDialog();
parent.refresh("general");
</script>	
END
				}
				elsif($pid == 0)
				{
					close (STDOUT);
					my $updateDartToRunning = $dbh->do("UPDATE matrix SET barcode = -1 WHERE id = $dartId");			
					if($dartFilePath)
					{
						`cp $dartFilePath $dartInfile`;
					}
					else
					{
						open (FILE, ">$dartInfile");
						while (read ($dartFile, my $Buffer, 1024)) {
							print FILE $Buffer;
						}
						close FILE;
					}
					`perl -p -i -e 's/\r/\n/g' $dartInfile`; #convert CR line terminators (MAC OS) into LF line terminators (Unix)
					my @headerRows = ("Order number","DArT plate barcode","client plate barcode","well row position","well column position","sample comments","genotype name");
					my $snpNumber = 0;
					my $genotype;
					my $genotypeNumber = 0;
					my $headerRow = 0;
					my $metadataCol = 0;
					my @metadata;
					open (INFILE, "$dartInfile");
					while (<INFILE>) {
						chop;
						/^\s*$/ and next; #skip blank line;
						my @dartLine =  split("\t",$_);
						if($dartLine[0] eq '*')
						{
							$metadataCol = 0;
							$genotypeNumber = 0;
							foreach (@dartLine)
							{
								if($_ eq '*')
								{
									$metadataCol++;
								}
								else
								{
									$genotype->{$genotypeNumber}->{$headerRows[$headerRow]} = $_;
									$genotypeNumber++;
								}
								
							}
							$headerRow++;
						}
						else
						{
							if($dartLine[0] =~ /ID/) # any words with 'ID', e.g. CloneID, SNP-ID, etc.
							{
								my $metadataLineCol = 0;
								$genotypeNumber = 0;
								foreach (@dartLine)
								{
									if($metadataLineCol < $metadataCol)
									{
										push @metadata, $_;
									}
									else
									{
										$genotype->{$genotypeNumber}->{'genotype name'} = $_;
										$genotypeNumber++;
									}
									$metadataLineCol++;
								}
							}
							else
							{
								my $metadataSNPCol = 0;
								my $dartSNP;
								my $dartSNPId = 0;
								my $dartCloneId = '';
								$genotypeNumber = 0;
								foreach (@dartLine)
								{
									if($metadataSNPCol < $metadataCol)
									{
										$dartCloneId = $_ if (!$dartCloneId);
										$dartSNP->{$metadata[$metadataSNPCol]} = $_;
									}
									else
									{
										if(!$dartSNPId)
										{
											my $dartSNPEncoded = $json->encode($dartSNP);
											my $insertDartSNP=$dbh->prepare("INSERT INTO matrix VALUES ('', 'dartSNP', ?, 0, 0, 0, ?, 0, ?, ?, NOW())");
											$insertDartSNP->execute($dartCloneId,$dartId,$dartSNPEncoded,$userName);
											$dartSNPId = $dbh->{mysql_insertid};
										}
										$genotype->{$genotypeNumber}->{$dartSNPId} = $_;
										$genotypeNumber++;
									}
									$metadataSNPCol++;
								}
								$snpNumber++;
							}
						}
					}
					close INFILE;
					for (sort {$a <=> $b} keys %$genotype)
					{
						my $dartGenotypeEncoded = $json->encode($genotype->{$_});
						my $insertDartGenotype=$dbh->prepare("INSERT INTO matrix VALUES ('', 'dartGenotype', ?, 0, 0, 0, ?, 0, ?, ?, NOW())");
						$insertDartGenotype->execute($genotype->{$_}->{'genotype name'},$dartId,$dartGenotypeEncoded,$userName);
					}
					unlink ($dartInfile);
					my $updateDartToLoaded = $dbh->do("UPDATE matrix SET o = $snpNumber, x = $genotypeNumber, barcode = 1, creationDate = NOW() WHERE id = $dartId");			
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
	parent.errorPop("The given dart name is existing!");
</script>	
END
				exit;
			}
		}
		else
		{
			print <<END;
<script>
	parent.errorPop("No Dart file found!");
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
	parent.errorPop("Please give a dart name!");
</script>	
END
}

