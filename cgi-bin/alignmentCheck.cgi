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

my $assemblyId = param ('assemblyId') || '';
my $seqOne = param ('seqOne') || '';
my $seqTwo = param ('seqTwo') || '';
my $redo = param ('redo') || '0';
my $hidden = param ('hidden') || '1';
my $filter = param ('filter') || '0';

my $alignmentCheckFormUrl = "alignmentCheckForm.cgi";
if($assemblyId)
{
	$alignmentCheckFormUrl .= "?assemblyId=$assemblyId";
}
my $blastTwoseqFormUrl = "blastTwoseqForm.cgi";
if($assemblyId)
{
	$blastTwoseqFormUrl .= "?assemblyId=$assemblyId";
}

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

if($seqOne || $seqTwo)
{
	if($seqOne && $seqTwo)
	{
		unlink ("$commoncfg->{TMPDIR}/$seqOne-$seqTwo.aln.html") if $redo;
		if(!-e "$commoncfg->{TMPDIR}/$seqOne-$seqTwo.aln.html")
		{
			my $getSequenceA = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
			$getSequenceA->execute($seqOne);
			my @getSequenceA =  $getSequenceA->fetchrow_array();
			$getSequenceA[5] = commify($getSequenceA[5]);
			my $getSequenceB = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
			$getSequenceB->execute($seqTwo);
			my @getSequenceB = $getSequenceB->fetchrow_array();
			$getSequenceB[5] = commify($getSequenceB[5]);
			open (ALN,">$commoncfg->{TMPDIR}/$seqOne-$seqTwo.aln.html") or die "can't open file: $commoncfg->{TMPDIR}/$seqOne-$seqTwo.aln.html";
			print ALN <<END;
	<table id='alignment$$' class='display'><thead><tr><th>$getSequenceA[2]</th><th>Query Length</th><th>$getSequenceB[2]</th><th>Subject Length</th><th>Identity %</th><th>Alignment Length</th><th>Query Start</th><th>Query End</th><th></th><th>Subject Start</th><th>Subject End</th></tr></thead><tbody>
END
			my $queryDir;
			for (my $position = 0; $position < length($getSequenceA[0]); $position += 2)
			{
				$queryDir .= "/q". substr($getSequenceA[0],$position,2);
			}
			my $subjectDir;
			for (my $position = 0; $position < length($getSequenceB[0]); $position += 2)
			{
				$subjectDir .= "/s". substr($getSequenceB[0],$position,2);
			}
			open (TBL,"$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$getSequenceA[0]-$getSequenceB[0].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$getSequenceA[0]-$getSequenceB[0].tbl";
			while(<TBL>)
			{
				/^#/ and next;
				my @hit = split("\t",$_);
				next if ($hit[12] > 0);
				my $direction = ($hit[8] < $hit[9]) ? "+":"-";
				$hit[3] = commify($hit[3]);				
				$hit[6] = commify($hit[6]);
				$hit[7] = commify($hit[7]);
				$hit[8] = commify($hit[8]);
				$hit[9] = commify($hit[9]);
				$hit[6] = "<a class='ui-state-error-text' title='Seq End'>$hit[6]</a>" if($hit[6] eq 1 || $hit[6] eq $getSequenceA[5]);
				$hit[7] = "<a class='ui-state-error-text' title='Seq End'>$hit[7]</a>" if($hit[7] eq 1 || $hit[7] eq $getSequenceA[5]);
				$hit[8] = "<a class='ui-state-error-text' title='Seq End'>$hit[8]</a>" if($hit[8] eq 1 || $hit[8] eq $getSequenceB[5]);
				$hit[9] = "<a class='ui-state-error-text' title='Seq End'>$hit[9]</a>" if($hit[9] eq 1 || $hit[9] eq $getSequenceB[5]);
				print ALN "<tr>
					<td><a onclick='closeDialog();openDialog(\"seqView.cgi?seqId=$getSequenceA[0]\")' title='View this sequence'>$seqType{$getSequenceA[3]}</a> ($getSequenceA[4]) $bacAssignType{$getSequenceA[7]}</td>
					<td>$getSequenceA[5]</td>
					<td><a onclick='closeDialog();openDialog(\"seqView.cgi?seqId=$getSequenceB[0]\")' title='View this sequence'>$seqType{$getSequenceB[3]}</a> ($getSequenceB[4]) $bacAssignType{$getSequenceB[7]}</td>
					<td>$getSequenceB[5]</td>
					<td><a title='E-value:$hit[10] \nBit-score:$hit[11]'>$hit[2]</a></td>
					<td>$hit[3]</td>
					<td>$hit[6]</td>
					<td>$hit[7]</td>
					<td>$direction</td>
					<td>$hit[8]</td>
					<td>$hit[9]</td>
					</tr>";
			}
			close(TBL);

			$alignmentCheckFormUrl .= "&seqOne=$seqOne&seqTwo=$seqTwo";
			$blastTwoseqFormUrl .= "&seqOne=$seqOne&seqTwo=$seqTwo";
			print ALN <<END;
</tbody></table>
<script>
\$('#dialog').dialog("option", "title", "Alignment Between $getSequenceA[2]($seqOne) & $getSequenceB[2]($seqTwo)");
\$( "#dialog" ).dialog( "option", "buttons", [ { text: "New Check", click: function() {  closeDialog();openDialog('$alignmentCheckFormUrl'); } }, { text: "reRun BLAST", click: function() { closeDialog();openDialog('$blastTwoseqFormUrl'); } }, { text: "OK", click: function() {closeDialog(); } } ] );
\$('#dialog').dialog("option", "width", 1000);
\$( "#alignment$$" ).dataTable({
	"scrollY": "400px",
	"scrollCollapse": true,
	"paging": false,
	"searching": false
});
</script>
END
			close (ALN);
		}
		print <<END;
		<script>
		parent.loadingHide();
		parent.closeDialog();
		parent.openDialog('$commoncfg->{TMPURL}/$seqOne-$seqTwo.aln.html');
		</script>
END
	}
	else
	{
		my $querySeq = ($seqOne) ? $seqOne : $seqTwo;
		unlink ("$commoncfg->{TMPDIR}/$querySeq.aln.html") if $redo;
		if(!-e "$commoncfg->{TMPDIR}/$querySeq.aln.html")
		{
			my $getSequenceA = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
			$getSequenceA->execute($querySeq);
			my @getSequenceA =  $getSequenceA->fetchrow_array();
			$getSequenceA[5] = commify($getSequenceA[5]);
			open (ALN,">$commoncfg->{TMPDIR}/$querySeq.aln.html") or die "can't open file: $commoncfg->{TMPDIR}/$querySeq.aln.html";
			print ALN <<END;
<table id='alignment$$' class='display'><thead><tr><th>Query</th><th>Query Length</th><th>Subject</th><th>Subject Length</th><th>Identity %</th><th>Alignment Length</th><th>Query Start</th><th>Query End</th><th></th><th>Subject Start</th><th>Subject End</th></tr></thead><tbody>
END

			my $queryDir;
			for (my $position = 0; $position < length($getSequenceA[0]); $position += 2)
			{
				$queryDir .= "/q". substr($getSequenceA[0],$position,2);
			}

			my @dirs;
			push @dirs, "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir";
			my @dirFile;

			do{
				my $dir = shift @dirs;
				opendir(DIR, $dir) or die "can't opendir $dir: $!";
				my @files = readdir(DIR);
				closedir DIR;
				foreach my $file (sort @files)
				{
					next if ($file =~ /^\./);
					if(-f "$dir/$file")
					{
						push @dirFile, "$dir/$file";
					}
					elsif(-d "$dir/$file")
					{
						push @dirFile, "$dir/$file";
						unshift @dirs, "$dir/$file";
					}
				}
			}while (@dirs);

			foreach my $file (sort @dirFile)
			{
				my $getSequenceBline;
				open (TBL, $file) or die "can't open file: $file";
				while(<TBL>)
				{
					/^#/ and next;
					my @hit = split("\t",$_);
					next if ($hit[12] > 0);
					my $direction = ($hit[8] < $hit[9]) ? "+":"-";
					$hit[3] = commify($hit[3]);
					unless (exists $getSequenceBline->{$hit[1]})
					{
						my $getSequenceB = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
						$getSequenceB->execute($hit[1]);
						@{$getSequenceBline->{$hit[1]}} = $getSequenceB->fetchrow_array();
						$getSequenceBline->{$hit[1]}[5] = commify($getSequenceBline->{$hit[1]}[5]);
					}
					next if ($filter && ($getSequenceA[2] eq $getSequenceBline->{$hit[1]}[2]));
					$hit[6] = commify($hit[6]);
					$hit[7] = commify($hit[7]);
					$hit[8] = commify($hit[8]);
					$hit[9] = commify($hit[9]);
					$hit[6] = "<a class='ui-state-error-text' title='Seq End'>$hit[6]</a>" if($hit[6] eq 1 || $hit[6] eq $getSequenceA[5]);
					$hit[7] = "<a class='ui-state-error-text' title='Seq End'>$hit[7]</a>" if($hit[7] eq 1 || $hit[7] eq $getSequenceA[5]);
					$hit[8] = "<a class='ui-state-error-text' title='Seq End'>$hit[8]</a>" if($hit[8] eq 1 || $hit[8] eq $getSequenceBline->{$hit[3]}[5]);
					$hit[9] = "<a class='ui-state-error-text' title='Seq End'>$hit[9]</a>" if($hit[9] eq 1 || $hit[9] eq $getSequenceBline->{$hit[3]}[5]);
					print ALN "<tr>
						<td><a onclick='closeDialog();openDialog(\"seqView.cgi?seqId=$getSequenceA[0]\")' title='View this sequence'>$seqType{$getSequenceA[3]}</a> ($getSequenceA[4]) $bacAssignType{$getSequenceA[7]}</td>
						<td>$getSequenceA[5]</td>
						<td><a onclick='closeDialog();openDialog(\"seqView.cgi?seqId=$getSequenceBline->{$hit[1]}[0]\")' title='View this sequence'>$getSequenceBline->{$hit[1]}[2]</a> $seqType{$getSequenceBline->{$hit[1]}[3]} ($getSequenceBline->{$hit[1]}[4]) $bacAssignType{$getSequenceBline->{$hit[1]}[7]}</td>
						<td>$getSequenceBline->{$hit[1]}[5]</td>
						<td><a title='E-value:$hit[10] \nBit-score:$hit[11]'>$hit[2]</a></td>
						<td>$hit[3]</td>
						<td>$hit[6]</td>
						<td>$hit[7]</td>
						<td>$direction</td>
						<td>$hit[8]</td>
						<td>$hit[9]</td>
						</tr>";
				}
				close(TBL);
			}

			$alignmentCheckFormUrl .= "&seqOne=$querySeq";
			$blastTwoseqFormUrl .= "&seqOne=$querySeq";
			print ALN <<END;
</tbody></table>
<script>
\$('#dialog').dialog("option", "title", "Alignment of $getSequenceA[2] ($querySeq)");
\$( "#dialog" ).dialog( "option", "buttons", [ { text: "New Check", click: function() {  closeDialog();openDialog('$alignmentCheckFormUrl'); } }, { text: "reRun BLAST", click: function() { closeDialog();openDialog('$blastTwoseqFormUrl'); } }, { text: "OK", click: function() {closeDialog(); } } ] );
\$('#dialog').dialog("option", "width", 1000);
\$( "#alignment$$" ).dataTable({
	"scrollY": "400px",
	"scrollCollapse": true,
	"paging": false
});
</script>
END
			close (ALN);
		}

		print <<END;
		<script>
		parent.loadingHide();
		parent.closeDialog();
		parent.openDialog('$commoncfg->{TMPURL}/$querySeq.aln.html');
		</script>
END
	}
}
else
{
	print <<END;
<script>
	parent.loadingHide();
	parent.errorPop('Please give at least one seq name!');
</script>	
END
}