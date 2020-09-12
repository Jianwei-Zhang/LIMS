#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
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
my $replace = param ('replace') || '0';
my $fpcFile = upload ('fpcFile');
my $infile = "$commoncfg->{TMPDIR}/$libraryId-$$.fpc";

print header;

if($fpcFile)
{
	my $pid = fork();
	if ($pid)
	{
		print <<END;
	<script>
	parent.closeDialog();
	</script>
END
	}
	elsif($pid == 0)
	{
		open (FILE, ">$infile");
		while (read ($fpcFile, my $Buffer, 1024)) {
			print FILE $Buffer;
		}
		close FILE;
		#connect to the mysql server
		my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});
		my $library=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$library->execute($libraryId);
		my @library = $library->fetchrow_array();

		my $fpc = '';
		my $version = 1;
		my $info;
		my $fpcClone;
		my $fpcMarker;
		my $Markerdata;
		my $fpcCtg;
		my $Contigdata;
		my $fpcSeq;
		my $Seqdata;
		my $bac;
		my $marker;
		my $ctg;
		my $seq;
		my $bacFlag = 1;
		my $markerFlag = 0;
		my $ctgFlag = 0;
		my $seqFlag = 0;	

		open (INFILE, "$infile") or die "can't open INFILE: $!";
		while (<INFILE>)
		{
			chop;
			/^\s*$/ and next; #skip blank line
			if(/^\/\//) #read remark part 
			{
				if(/fpc project (.*)/)
				{
					$fpc = $1;
					print <<END;
	<script>
	parent.informationPop("Loading process is running! This might take several minutes.");
	parent.refresh("menu");
	</script>
END
					close (STDOUT);
				}
				$info .= "$_ \n";
			}
			else #process data part
			{
				unless ($fpc)
				{
					print <<END;
	<script>
	parent.errorPop("Sorry, your input file is not in FPC format!");
	</script>
END
					unlink ($infile);
					exit;
				}
			
				if(/^Markerdata (.*)/)
				{
					$Markerdata = $1;
					$bacFlag = 0;
					$markerFlag = 1;
					$ctgFlag = 0;		
					$seqFlag = 0;	
					next;
				}
				elsif(/^Contigdata (.*)/)
				{
					$Contigdata = $1;
					$bacFlag = 0;
					$markerFlag = 0;
					$ctgFlag = 1;
					$seqFlag = 0;	
					next;
				}
				elsif(/^Seqdata (.*)/)
				{
					$Seqdata = $1;
					$bacFlag = 0;
					$markerFlag = 0;
					$ctgFlag = 0;	
					$seqFlag = 1;
					next;
				}

				if($bacFlag)
				{
					if(/^BAC : \"(.*)\"/)
					{
						$bac = $1;
						$bac =~ s/[^a-zA-Z0-9]//g;
						$bac =~ /(\d+)(\D+)(\d+)$/;
						my $plateName =  sprintf "%0*d", 4, $1;
						my $row =  uc $2;
						my $col =  sprintf "%0*d", 2, $3;
						$bac = "$library[2]$plateName$row$col";
						$fpcClone->{$bac} = "BAC : \"$bac\"\n";
						next;
					}
					$fpcClone->{$bac} .= "$_\n";
				}
	
				if($markerFlag)
				{
					#do nothing, to be added;
				}

				if($ctgFlag)
				{
					if(/^(Ctg\d+) (.*)/)
					{
						$ctg = $1;
					}
					$fpcCtg->{$ctg} .= "$_\n";
				}

				if($seqFlag)
				{
					#do nothing, to be added;
				}
			}
		}
		close (INFILE);

		my $fpcId;
		my @fpcId = ();
		my $existsFpc = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'fpc' AND name LIKE ? AND barcode = ? ORDER BY o");
		$existsFpc->execute($fpc,$libraryId);
		while(my @row = $existsFpc->fetchrow_array())
		{
			push(@fpcId, $row[0]);
			$version = $row[3] + 1 if($row[3] >= $version);
		}

		if(@fpcId > 0 && $replace > 0)
		{
				$fpcId = pop @fpcId;		
				my $updateFpc = $dbh->prepare("UPDATE matrix SET note = ? WHERE id = ?");
				$updateFpc->execute($info,$fpcId);
				my $deleteFpcClone=$dbh->do("DELETE FROM matrix WHERE o = $fpcId AND container LIKE 'fpcClone'");	
				my $deleteFpcCtg=$dbh->do("DELETE FROM matrix WHERE o = $fpcId AND container LIKE 'fpcCtg'");	
		}
		else
		{
			my $insertFpc=$dbh->prepare("INSERT INTO matrix VALUES ('', 'fpc', ?, ?, 0, 0, 0, ?, ?, ?, NOW())");
			$insertFpc->execute($fpc,$version,$libraryId,$info,$userName);
			$fpcId = $dbh->{mysql_insertid};
		}
		my $updateFpcMarkerdata = $dbh->do("UPDATE matrix SET x = $Markerdata WHERE id = $fpcId") if ($Markerdata);
		my $updateFpcContigdata = $dbh->do("UPDATE matrix SET y = $Contigdata WHERE id = $fpcId") if ($Contigdata);
		my $updateFpcSeqdata = $dbh->do("UPDATE matrix SET z = $Seqdata WHERE id = $fpcId") if ($Seqdata);

		foreach my $clone (sort keys %$fpcClone)
		{
			my $cloneInLibrary=$dbh->prepare("SELECT * FROM clones WHERE name LIKE ? OR origin LIKE ? LIMIT 1");
			$cloneInLibrary->execute($clone,$clone);
			my @cloneInLibrary = $cloneInLibrary->fetchrow_array();
			my $sequenced = $cloneInLibrary[6] || 0;
			my $insertFpcClone = $dbh->prepare("INSERT INTO matrix VALUES ('','fpcClone', ?, ?, ?, 0, 0, 0, ?, ?, NOW())");
			$insertFpcClone->execute($clone,$fpcId,$sequenced,$fpcClone->{$clone},$userName);
		}

		foreach my $contig (sort keys %$fpcCtg)
		{
			my $fpcCloneNumber =  $dbh->prepare("SELECT COUNT(*) FROM matrix WHERE container LIKE 'fpcClone' AND o = ? AND MATCH (note) AGAINST (?)");
			$fpcCloneNumber->execute($fpcId,$contig);
			my @fpcCloneNumber = $fpcCloneNumber->fetchrow_array();
			my $insertFpcCtg = $dbh->prepare("INSERT INTO matrix VALUES ('','fpcCtg', ?, ?, ?, 0, 0, 0, ?, ?, NOW())");
			$insertFpcCtg->execute($contig,$fpcId,$fpcCloneNumber[0],$fpcCtg->{$contig},$userName);
		}
		unlink ($infile);
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
	parent.errorPop("No FPC file found!");
</script>
END
}


