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
my $replace = param ('replace') || 0;
my $jobFile = upload ('jobFile');
my $infile = "$commoncfg->{TMPDIR}/$$.job";

print header;
if($jobFile)
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
		while (read ($jobFile, my $Buffer, 1024)) {
			print FILE $Buffer;
		}
		close FILE;
		#connect to the mysql server
		my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});
		my $job;
		my $jobList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'job'");
		$jobList->execute();
		while (my @jobList = $jobList->fetchrow_array())
		{
			if($jobList[2] =~ /^-/)
			{
				$jobList[2] =~ s/-//g;
			}
			@{$job->{$jobList[2]}} = @jobList;
		}
		open (INFILE, "$infile") or die "can't open INFILE: $!";
		my @allJobsInFile = <INFILE>;
		my $allJobs  = decode_json (join "\n", @allJobsInFile);
		for my $eachJob (@{$allJobs->{'rows'}})
		{
			next if ($eachJob->{'jobStatus'}  ne 'Completed');
			my $jobDetails = encode_json $eachJob;
			my $jobId = $eachJob->{'jobId'};
			if (exists $job->{$jobId})
			{
				my $updateJob=$dbh->do("UPDATE matrix SET note = '$jobDetails' WHERE container LIKE 'job' AND (name LIKE '$jobId' OR name LIKE '-$jobId')") if($replace);
			}
			else
			{
				my $insertJob=$dbh->prepare("INSERT INTO matrix VALUES ('', 'job', ?, 0, 0, 0, 0, -3, ?, ?, NOW())");
				$insertJob->execute($jobId,$jobDetails,$userName);
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
	parent.errorPop("No Job List Found!");
</script>
END
}

