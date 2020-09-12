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

my $type = param ('type') || '';
my $parentId = param ('parentId') || '0';

my $tableFile = upload ('tableFile');
my $tableFilePath = param ('tableFilePath') || '';
my $replace = param ('replace') || '0';
my $idColumn = param ('idColumn') || '1';
$idColumn = $idColumn - 1;
my $refresh = param ('refresh') || 'menu';

my $safeFilenameCharacters = "a-zA-Z0-9_.-";
my $tableInfile = "$commoncfg->{TMPDIR}/$$.table";
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

if($type)
{
	if($parentId)
	{
		if($tableFile || $tableFilePath)
		{
			my $pid = fork();
			if ($pid)
			{
				print <<END;
<script>
parent.closeDialog();
parent.refresh("$refresh");
</script>	
END
			}
			elsif($pid == 0)
			{
				close (STDOUT);
				if($tableFilePath)
				{
					#$tableInfile = $tableFilePath;
					`cp $tableFilePath $tableInfile`;
				}
				else
				{
					open (FILE, ">$tableInfile");
					while (read ($tableFile, my $Buffer, 1024)) {
						print FILE $Buffer;
					}
					close FILE;
				}
				`perl -p -i -e 's/\r/\n/g' $tableInfile`; #convert CR line terminators (MAC OS) into LF line terminators (Unix)

				my $deleteTable = $dbh->do("DELETE FROM matrix WHERE container LIKE '$type' AND z = $parentId") if ($replace);
				my $line = 0;
				my $bodyRow;
				open (INFILE, "$tableInfile");
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
						my $insertTableRow=$dbh->prepare("INSERT INTO matrix VALUES ('', ?, ?, 0, 0, 0, ?, 0, ?, ?, NOW())");
						$insertTableRow->execute($type,$tableLine[$idColumn],$parentId,$tableRowEncoded,$userName);
						$line++;
					}
				}
				close INFILE;
				unlink ($tableInfile);
				my $countTable = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE '$type' AND z = $parentId");
				$countTable->execute();
				$line = $countTable->rows;
				my $updateParent = $dbh->do("UPDATE matrix SET o = $line, creationDate = NOW() WHERE id = $parentId");
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
	parent.errorPop("No Table file found!");
</script>
END
			exit;
		}
	}
	else
	{
	print <<END;
<script>
	parent.errorPop("Invalid operation!");
</script>	
END
	}
}
else
{
	print <<END;
<script>
	parent.errorPop("Please provide a table type!");
</script>	
END
}

