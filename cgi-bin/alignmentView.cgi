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

unless (-e $commoncfg->{TMPDIR})
{
	mkdir $commoncfg->{TMPDIR};
}
my %seqType = (
	0=>'Assembled',
	1=>'BAC-Insert',
	2=>'BAC-Circularized',
	3=>'BAC-NonVector',
	4=>'BAC-Gapped',
	5=>'Partial',
	6=>'Vector/Mixer',
	7=>'Mixer',
	8=>'SHORT',
	97=>'Piece',
	98=>'BES',
	99=>'Genome'
	);
my %bacAssignType = (
	0=>'',
	1=>'TagValid',
	2=>'BesValid',
	3=>'TagValid+BesValid',
	4=>'TagForced'
	);
my %seqDir = (
	0 => "NA",
	1 => "f",
	2 => "r"
	);

print header;

if($alignmentId)
{
	unlink ("$commoncfg->{TMPDIR}/$alignmentId.aln.html") if(-e "$commoncfg->{TMPDIR}/$alignmentId.aln.html");
	if(!-e "$commoncfg->{TMPDIR}/$alignmentId.aln.html")
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
		open (TBL,"$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$alignmentId[0]-$alignmentId[1].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$alignmentId[0]-$alignmentId[1].tbl";
		while(<TBL>)
		{
			chop;
			/^#/ and next;
			$alignmentCount++;
			@getAlignment = split("\t",$_);
			last if ($alignmentCount == $alignmentId[2]);
		}
		close(TBL);
		
		my $queryDirection = ($getAlignment[6] < $getAlignment[7]) ? "+" : "-";
		my $subjectDirection = ($getAlignment[8] < $getAlignment[9]) ? "+" : "-";
		my $getSequenceA = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$getSequenceA->execute($getAlignment[0]);
		my @getSequenceA =  $getSequenceA->fetchrow_array();
		$getSequenceA[5] = commify($getSequenceA[5]);
		my $getSequenceB = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$getSequenceB->execute($getAlignment[1]);
		my @getSequenceB = $getSequenceB->fetchrow_array();
		$getSequenceB[5] = commify($getSequenceB[5]);
		$getAlignment[3] = commify($getAlignment[3]);				
		$getAlignment[6] = commify($getAlignment[6]);
		$getAlignment[7] = commify($getAlignment[7]);
		$getAlignment[8] = commify($getAlignment[8]);
		$getAlignment[9] = commify($getAlignment[9]);
		$getAlignment[6] = "<a class='ui-state-error-text' title='Seq End'>$getAlignment[6]</a>" if($getAlignment[6] eq 1 || $getAlignment[6] eq $getSequenceA[5]);
		$getAlignment[7] = "<a class='ui-state-error-text' title='Seq End'>$getAlignment[7]</a>" if($getAlignment[7] eq 1 || $getAlignment[7] eq $getSequenceA[5]);
		$getAlignment[8] = "<a class='ui-state-error-text' title='Seq End'>$getAlignment[8]</a>" if($getAlignment[8] eq 1 || $getAlignment[8] eq $getSequenceB[5]);
		$getAlignment[9] = "<a class='ui-state-error-text' title='Seq End'>$getAlignment[9]</a>" if($getAlignment[9] eq 1 || $getAlignment[9] eq $getSequenceB[5]);
		
		open (ALN,">$commoncfg->{TMPDIR}/$alignmentId.aln.html") or die "can't open file: $commoncfg->{TMPDIR}/$alignmentId.aln.html";
		print ALN <<END;
	<table id='alignment$alignmentId$$' class='display'>
	<thead><tr><th></th><th>Seqs</th><th>Length</th><th>Identity</th><th>Alignment Length</th><th>Start</th><th></th><th>End</th></tr></thead>
	<tbody>
		<tr>
			<td><b>Query</b></td>
			<td>$getSequenceA[2] <sup onclick='closeDialog();openDialog("seqView.cgi?seqId=$getSequenceA[0]")' title='View this sequence'>$seqType{$getSequenceA[3]}</sup> $bacAssignType{$getSequenceA[7]}</td>
			<td>$getSequenceA[5]</td>
			<td rowspan='2'><a title='E-value:$getAlignment[10] \nBit-score:$getAlignment[11]'>$getAlignment[2]%</a></td>
			<td rowspan='2'>$getAlignment[3]</td>
			<td>$getAlignment[6]</td>
			<td>$queryDirection</td>
			<td>$getAlignment[7]</td>
		</tr>
		<tr>
			<td><b>Subject</b></td>
			<td>$getSequenceB[2] <sup onclick='closeDialog();openDialog("seqView.cgi?seqId=$getSequenceB[0]")' title='View this sequence'>$seqType{$getSequenceB[3]}</sup> $bacAssignType{$getSequenceB[7]}</td>
			<td>$getSequenceB[5]</td>
			<td>$getAlignment[8]</td>
			<td>$subjectDirection</td>
			<td>$getAlignment[9]</td>
		</tr>
	</tbody>
	</table>
<script>
\$('#dialog').dialog("option", "title", "Selected Alignment Between $getSequenceA[2]($getSequenceA[0]) & $getSequenceB[2]($getSequenceB[0])");
\$( "#dialog" ).dialog( "option", "buttons", [ { text: "OK", click: function() {closeDialog(); } } ] );
</script>
END
		close (ALN);
	}
	print <<END;
	<script>
	parent.loadingHide();
	parent.closeDialog();
	parent.openDialog('$commoncfg->{TMPURL}/$alignmentId.aln.html');
	</script>
END
	
}
else
{
	print <<END;
<script>
	parent.loadingHide();
	parent.errorPop('Please give an alignment Id!');
</script>	
END
}