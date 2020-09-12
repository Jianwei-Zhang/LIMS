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
use userCookie;

my $userCookie = new userCookie;
my $userId = (cookie('cid')) ? $userCookie->checkCookie(cookie('cid')) : 0;
exit if (!$userId);

my $user = new user;
my $userDetail = $user->getAllFieldsWithUserId($userId);
my $userName = $userDetail->{"userName"};

my $commoncfg = readConfig("main.conf");
my $libraryId = param ('libraryId') || '';
my $reformat = param ('reformat') || '0';

print header;
if($libraryId)
{
	my $pid = fork();
	if ($pid)
	{
		print <<END;
	<script>
		parent.closeDialog();
		parent.informationPop("It's running!");
	</script>	
END
	}
	elsif($pid == 0)
	{
		close (STDOUT);
		#connect to the mysql server
		my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});
		my $delete = $dbh->do("DELETE FROM clones WHERE libraryId = $libraryId");
		my %wellx = ( 96 => ['01' .. '12'], 384 =>['01' .. '24']);
		my %welly = ( 96 => ['A' .. 'H'], 384 =>['A' .. 'P']);
		my $library=$dbh->prepare("SELECT * FROM matrix WHERE id = $libraryId");
		$library->execute();
		my @library = $library->fetchrow_array();
		my $libraryDetails = decode_json $library[8];
		my $plateId;
		my $plateDescription;
		my $plateInLibrary=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'plate' AND z = $libraryId ORDER BY o");
		$plateInLibrary->execute();
		while (my @plateInLibrary = $plateInLibrary->fetchrow_array())
		{
			$plateId->{$plateInLibrary[2]}->{$plateInLibrary[4]} = $plateInLibrary[0];
			$plateDescription->{$plateInLibrary[2]} = $plateInLibrary[8];
		}

		open(TMP,">$commoncfg->{TMPDIR}/$$.clone");
		if($library[5]) #rearrayed lib
		{
			my $sourceLibrary=$dbh->prepare("SELECT * FROM matrix WHERE id = $library[5]");
			$sourceLibrary->execute();
			my @sourceLibrary = $sourceLibrary->fetchrow_array();
			my $sourceLibraryDetails = decode_json $sourceLibrary[8];

			for (sort keys %{$plateId})
			{
				my $plate = $_;
				my @sourceClone = split /\s+/, $plateDescription->{$plate};
				my $wellNumber = 0;
				for(sort @{$welly{$libraryDetails->{'format'}}})
				{
					my $y = $_;
					for(sort @{$wellx{$libraryDetails->{'format'}}})
					{
						$sourceClone[$wellNumber] =~ s/[^a-zA-Z0-9]//g;
						if($sourceClone[$wellNumber] =~ /(\d+)(\D+)(\d+)$/)
						{
							my $plateName =  sprintf "%0*d", 4, $1;
							my $row =  uc $2;
							my $col =  sprintf "%0*d", 2, $3;
							my $originCloneName = ($reformat) ? "$sourceLibrary[2]$plateName$row$col" : $sourceClone[$wellNumber];
							my $cloneName = "$library[2]$plate$y$_";
							my $sequence=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND (name LIKE ? OR name LIKE ?) AND o < 50 ORDER BY o,x DESC LIMIT 1");
							$sequence->execute($cloneName,$originCloneName);
							my @sequence = $sequence->fetchrow_array();
							my $sequenced = $sequence[3] || 0;
							print  TMP "$cloneName\t$libraryId\t$plate\t$y$_\t$originCloneName\t$sequenced\n";
						}
						else
						{
							my $originCloneName = "EMPTY";
							my $cloneName = "$library[2]$plate$y$_";
							my $sequence=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND (name LIKE ? OR name LIKE ?) AND o < 50 ORDER BY o,x DESC LIMIT 1");
							$sequence->execute($cloneName,$originCloneName);
							my @sequence = $sequence->fetchrow_array();
							my $sequenced = $sequence[3] || 0;
							print  TMP "$cloneName\t$libraryId\t$plate\t$y$_\t$originCloneName\t$sequenced\n";
						}
						$wellNumber++;
					}
				}
			}		
		}
		else
		{
			for (sort keys %{$plateId})
			{
				my $plate = $_;
				for(sort @{$welly{$libraryDetails->{'format'}}})
				{
					my $y = $_;
					for(sort @{$wellx{$libraryDetails->{'format'}}})
					{
						my $cloneName = "$library[2]$plate$y$_";
						my $sequence=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND name LIKE ? AND o < 50 ORDER BY o,x DESC LIMIT 1");
						$sequence->execute($cloneName);
						my @sequence = $sequence->fetchrow_array();
						my $sequenced = $sequence[3] || 0;
						print  TMP "$cloneName\t$libraryId\t$plate\t$y$_\t\t$sequenced\n";
					}
				}
			}
		}
		close(TMP);
		my $loadClones = $dbh->do("LOAD DATA LOCAL INFILE '$commoncfg->{TMPDIR}/$$.clone' INTO TABLE clones (name, libraryId, plate, well, origin, sequenced)");
		my $toDelete = 1;
		do {
			$toDelete = 0;
			unlink("$commoncfg->{TMPDIR}/$$.clone");
			$toDelete = 1 if (-e "$commoncfg->{TMPDIR}/$$.clone");
			sleep 3;
		} while($toDelete);
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
	parent.errorPop("Not a valid operation!");
</script>	
END
}
