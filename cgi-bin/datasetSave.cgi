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

my $datasetId = param ('datasetId') || '';
my $datasetName = param ('name') || '';
my $datasetType = param ('datasetType') || '0';
my $parentId = param ('parentId') || '0';
my $idColumn = param ('idColumn') || '1';
$idColumn = $idColumn - 1;

my $datasetDescription = param('description') || '';
my $datasetFile = upload ('datasetFile');
my $datasetFilePath = param ('datasetFilePath') || '';

my $safeFilenameCharacters = "a-zA-Z0-9_.-";
my $datasetInfile = "$commoncfg->{TMPDIR}/$$.dataset";
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

if($datasetName)
{
	if($datasetId)
	{
		my $checkDatasetName = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'dataset' AND name LIKE ? AND id != ?");
		$checkDatasetName->execute($datasetName, $datasetId);
		if($checkDatasetName->rows < 1)
		{
			my $checkCreator = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
			$checkCreator->execute($datasetId);
			my @checkCreator=$checkCreator->fetchrow_array();
			if($checkCreator[9] eq $userName || exists $userPermission->{$userId}->{'dataset'})
			{
				my $updateDataset=$dbh->prepare("UPDATE matrix SET name = ?, o = ?, z = ?, note = ? WHERE id = ?");
				$updateDataset->execute($datasetName,$datasetType,$parentId,$datasetDescription,$datasetId);
				if($datasetFile || $datasetFilePath)
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
						my $updateDatasetToRunning = $dbh->do("UPDATE matrix SET barcode = -1 WHERE id = $datasetId");			
						if($datasetFilePath)
						{
							`cp $datasetFilePath $datasetInfile`;
						}
						else
						{
							open (FILE, ">$datasetInfile");
							while (read ($datasetFile, my $Buffer, 1024)) {
								print FILE $Buffer;
							}
							close FILE;
						}		
						`perl -p -i -e 's/\r/\n/g' $datasetInfile`; #convert CR line terminators (MAC OS) into LF line terminators (Unix)
						#delete old data set
						my $deleteTable = $dbh->do("DELETE FROM matrix WHERE container LIKE 'record' AND z = $datasetId");
						my $bodyRow;
						my $line = 0;
						open (INFILE, "$datasetInfile");
						while (<INFILE>) {
							chop;
							$_ .= "\tEND"; #add a col
							my @tableLine =  split("\t",$_);
							pop @tableLine; #remove added col
							next unless ($tableLine[$idColumn]);
							my $column = 0;
							if($line < 1)
							{
								for (@tableLine)
								{
									$bodyRow->{$column}->{'field'} = $_;
									$column++;
								}
								$line++;
								next;
							}
							else
							{
								foreach (@tableLine)
								{
									s/^"+//g;
									s/"+$//g;
									$bodyRow->{$column}->{'value'} = $_;
									$column++;
								}
								$tableLine[$idColumn] =~ s/^"+//g;
								$tableLine[$idColumn] =~ s/"+$//g;
								my $tableRowEncoded = $json->encode($bodyRow);
								my $insertTableRow=$dbh->prepare("INSERT INTO matrix VALUES ('', 'record', ?, 0, 0, 0, ?, 0, ?, ?, NOW())");
								$insertTableRow->execute($tableLine[$idColumn],$datasetId,$tableRowEncoded,$userName);
								$line++;
							}
						}
						close INFILE;
						$line--;
						unlink ($datasetInfile);
						my $updateDatasetToLoaded = $dbh->do("UPDATE matrix SET x = $line, barcode = 1, creationDate = NOW() WHERE id = $datasetId");			
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
	parent.errorPop("The given dataset name is existing!");
</script>	
END
		}
	}
	else
	{
		if($datasetFile || $datasetFilePath)
		{
			my $checkDatasetName = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'dataset' AND name LIKE ?");
			$checkDatasetName->execute($datasetName);
			if($checkDatasetName->rows < 1)
			{
				my $insertDataset=$dbh->prepare("INSERT INTO matrix VALUES ('', 'dataset', ?, ?, 0, 0, ?, 0, ?, ?, NOW())");
				$insertDataset->execute($datasetName,$datasetType,$parentId,$datasetDescription,$userName);
				$datasetId = $dbh->{mysql_insertid};
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
					my $updateDatasetToRunning = $dbh->do("UPDATE matrix SET barcode = -1 WHERE id = $datasetId");			
					if($datasetFilePath)
					{
						`cp $datasetFilePath $datasetInfile`;
					}
					else
					{
						open (FILE, ">$datasetInfile");
						while (read ($datasetFile, my $Buffer, 1024)) {
							print FILE $Buffer;
						}
						close FILE;
					}
					`perl -p -i -e 's/\r/\n/g' $datasetInfile`; #convert CR line terminators (MAC OS) into LF line terminators (Unix)
					my $bodyRow;
					my $line = 0;
					open (INFILE, "$datasetInfile");
					while (<INFILE>) {
						chop;
						$_ .= "\tEND"; #add a col
						my @tableLine =  split("\t",$_);
						pop @tableLine; #remove added col
						next unless ($tableLine[$idColumn]);
						my $column = 0;
						if($line < 1)
						{
							for (@tableLine)
							{
								$bodyRow->{$column}->{'field'} = $_;
								$column++;
							}
							$line++;
							next;
						}
						else
						{
							foreach (@tableLine)
							{
								s/^"+//g;
								s/"+$//g;
								$bodyRow->{$column}->{'value'} = $_;
								$column++;
							}
							$tableLine[$idColumn] =~ s/^"+//g;
							$tableLine[$idColumn] =~ s/"+$//g;
							my $tableRowEncoded = $json->encode($bodyRow);
							my $insertTableRow=$dbh->prepare("INSERT INTO matrix VALUES ('', 'record', ?, 0, 0, 0, ?, 0, ?, ?, NOW())");
							$insertTableRow->execute($tableLine[$idColumn],$datasetId,$tableRowEncoded,$userName);
							$line++;
						}
					}
					close INFILE;
					unlink ($datasetInfile);
					my $countTable = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'record' AND z = $datasetId");
					$countTable->execute();
					$line = $countTable->rows;
					my $updateDatasetToLoaded = $dbh->do("UPDATE matrix SET x = $line, barcode = 1, creationDate = NOW() WHERE id = $datasetId");			
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
	parent.errorPop("The given dataset name is existing!");
</script>	
END
				exit;
			}
		}
		else
		{
			print <<END;
<script>
	parent.errorPop("No Dataset file found!");
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
	parent.errorPop("Please give a dataset name!");
</script>	
END
}

