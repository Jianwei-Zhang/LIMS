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
my $replace = param ('replace') || 0;
my $tagFile = upload ('tagFile');
my $infile = "$commoncfg->{TMPDIR}/$libraryId-$$.tag";

print header;
if($tagFile)
{
	my $pid = fork();
	if ($pid)
	{
		print <<END;
	<script>
	parent.informationPop("Loading process is running! This might take several minutes.");
	parent.closeDialog();
	parent.refresh("menu");
	</script>
END
	}
	elsif($pid == 0)
	{
		close (STDOUT);
		open (FILE, ">$infile");
		while (read ($tagFile, my $Buffer, 1024)) {
			print FILE $Buffer;
		}
		close FILE;
		#connect to the mysql server
		my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});
		my $library=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$library->execute($libraryId);
		my @library = $library->fetchrow_array();
		if($replace)
		{
			my $deleteTag=$dbh->do("DELETE FROM matrix WHERE x = $libraryId AND container LIKE 'tag'");	
		}

		open (INFILE, "$infile") or die "can't open INFILE: $!";
		while (<INFILE>)
		{
			chop;
			my ($tagSeq,$tagId,$bacIds) = split(" ",$_);
			my @bacIds = split(",",$bacIds);
			for (@bacIds)
			{
				my $bac = $_;
				$bac =~ s/[^a-zA-Z0-9]//g;
				$bac =~ /(\d+)(\D+)(\d+)$/;
				my $plateName =  sprintf "%0*d", 4, $1;
				my $row =  uc $2;
				my $col =  sprintf "%0*d", 2, $3;
				$bac = "$library[2]$plateName$row$col";
				my $insertTag=$dbh->prepare("INSERT INTO matrix VALUES ('', 'tag', ?, ?, ?, 0, 0, 0, ?, ?, NOW())");
				$insertTag->execute($bac,$tagId,$libraryId,$tagSeq,$userName);
			}
		}
		close (INFILE);
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
	parent.errorPop("No WGP Tag file found!");
</script>
END
}

