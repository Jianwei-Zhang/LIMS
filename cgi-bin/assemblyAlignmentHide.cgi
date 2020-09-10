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
my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});

my $alignmentId = param ('alignmentId') || '';
my $assemblyId = param ('assemblyId') || '';
my $assemblyCtgId = param ('assemblyCtgId') || '';
my $openAssemblyId = param ('openAssemblyId') || '';
my $chr = param ('chr') || '0';
my $scrollLeft = param ('scrollLeft') || '0';

print header;

if($alignmentId)
{
	my @alignmentId = split "-", $alignmentId;
	my $queryDir;
	for (my $position = 0; $position < length($alignmentId[0]); $position += 2)
	{
		$queryDir .= "/q". substr($alignmentId[0],$position,2);
	}
	my $subjectDir;
	for (my $position = 0; $position < length($alignmentId[1]); $position += 2)
	{
		$subjectDir .= "/s". substr($alignmentId[1],$position,2);
	}
	my @getAlignment;
	my $alignmentCount = 0;
	open (NEW,">$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$alignmentId[0]-$alignmentId[1].new") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$alignmentId[0]-$alignmentId[1].new";
	open (TBL,"$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$alignmentId[0]-$alignmentId[1].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$alignmentId[0]-$alignmentId[1].tbl";
	while(<TBL>)
	{
		if (/^#/)
		{
			print NEW $_;
			next;
		}
		$alignmentCount++;
		@getAlignment = split("\t",$_);
		$getAlignment[12] =~ s/\W//g;
		$getAlignment[12] = 1 if ($alignmentCount == $alignmentId[2]);
		print NEW join "\t", @getAlignment;
		print NEW "\n";
	}
	close(TBL);
	close(NEW);
	unlink("$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$alignmentId[0]-$alignmentId[1].tbl");
	rename("$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$alignmentId[0]-$alignmentId[1].new", "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$alignmentId[0]-$alignmentId[1].tbl");

	if($openAssemblyId)
	{
			print <<END;
<script>
closeViewer();
openViewer("assemblyChrView.cgi?assemblyId=$openAssemblyId&chr=$chr&scrollLeft=$scrollLeft");
</script>	
END
	}
	else
	{
		if($chr)
		{
			print <<END;
<script>
closeViewer();
openViewer("assemblyChrView.cgi?assemblyId=$assemblyId&chr=$chr&scrollLeft=$scrollLeft");
</script>	
END
		}
		else
		{
			if($assemblyCtgId)
			{
			print <<END;
<script>
closeViewer();
openViewer("assemblyCtgView.cgi?assemblyCtgId=$assemblyCtgId&scrollLeft=$scrollLeft");
</script>	
END
			}
			else
			{
				print <<END;
<script>
	informationPop('Alignment is hidden now!');
</script>	
END
			}
		}
	}
}
else
{
	print <<END;
<script>
	errorPop('Not a valid operation!');
</script>	
END
}
