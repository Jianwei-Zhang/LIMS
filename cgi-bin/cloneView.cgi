#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use DBI;
use lib "lib/";
use lib "lib/pangu";
use pangu;
use userCookie;

my $userCookie = new userCookie;
my $userId = (cookie('cid')) ? $userCookie->checkCookie(cookie('cid')) : 0;
exit if (!$userId);

my $commoncfg = readConfig("main.conf");
my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});

my $cloneName = param ('cloneName') || '';

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
	0=>'NA',
	1=>'f',
	2=>'r'
	);

print header;
if ($cloneName)
{
	my $getClone = $dbh->prepare("SELECT * FROM clones WHERE name LIKE ?");
	$getClone->execute($cloneName);
	if($getClone->rows < 1)
	{
		print 'No valid clone Found!';
		exit;
	}
	else
	{
		my $relatedClone = '';
		my $wgpTags;
		my $fpcClone;
		my $fpcView = 'None.';
		my @getClone = $getClone->fetchrow_array();
		$getClone[5] = 'NA' unless ($getClone[5]);
		my $originalClone = ($getClone[5] ne 'NA' ) ? "<tr><td style='text-align:right;white-space: nowrap;'><b>Original Clone</b></td><td><a onclick='closeDialog();openDialog(\"cloneView.cgi?cloneName=$getClone[5]\")'>$getClone[5]</a></td></tr>" : "";
		my $getRelatedClone = $dbh->prepare("SELECT * FROM clones WHERE origin LIKE ?");
		$getRelatedClone->execute($cloneName);
		if ($getRelatedClone->rows > 0)
		{
			$relatedClone = "<tr><td style='text-align:right;white-space: nowrap;'><b>Related Clone</b></td><td>";
			while( my @getRelatedClone =  $getRelatedClone->fetchrow_array())
			{
				$relatedClone .= "<a onclick='closeDialog();openDialog(\"cloneView.cgi?cloneName=$getRelatedClone[1]\")'>$getRelatedClone[1]</a>\n";
			}
			$relatedClone .= "</td></tr>";
		}
# 		my $getTags = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'tag' AND (name LIKE ? OR name LIKE ?)");
# 		$getTags->execute($cloneName,$getClone[5]);
# 		$wgpTags = $getTags->rows;
		my $fpcViewTabs;
		my $getFpcClone = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'fpcClone' AND (name LIKE ? OR name LIKE ?)");
		$getFpcClone->execute($cloneName,$getClone[5]);
		while (my @getFpcClone = $getFpcClone->fetchrow_array())
		{
			$fpcView = "<div id='fpcCloneTabs$$'><ul>" unless ($fpcView ne 'None.');
			my $contig = 'Ctg0';
			if($getFpcClone[8] =~ /Map "(.*)" Ends Left/)
			{
				$contig = ucfirst ($1);
			}
			my $fpcList=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
			$fpcList->execute($getFpcClone[3]);
			my @fpcList = $fpcList->fetchrow_array();
			$fpcView .= "<li><a href='#fpcTabs-$fpcList[0]-$contig'>$fpcList[2] v.$fpcList[3] $contig</a></li>";
			my $fpcCtg=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'fpcCtg' AND o = ? AND Name LIKE ?");
			$fpcCtg->execute($fpcList[0],$contig);
			my @fpcCtg = $fpcCtg->fetchrow_array();
			$fpcViewTabs .= ($contig ne 'Ctg0') ? "<div id='fpcTabs-$fpcList[0]-$contig'>
				<a style='float: right;' onclick='closeDialog();closeViewer();openViewer(\"fpcCtgView.cgi?fpcCtgId=$fpcCtg[0]&highlight=$getFpcClone[2]\")'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-image' title='View'></span>View $contig</a><br>
				<pre>$getFpcClone[8]</pre></div>" : "<div id='fpcTabs-$fpcList[0]-$contig'><pre>$getFpcClone[8]</pre></div>";
		}
		$fpcView .= "</ul>$fpcViewTabs<div>" if($fpcView ne 'None.');

		my $sequenceList = "<table id='sequenceList$$' name='sequenceList$$' class='display'>
				<thead>
					<tr>
						<th>Sequence Type</th>
						<th>Length (bp)</th>
						<th>Related</th>
					</tr>
				</thead><tbody>";
		my $getSequences = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND (name LIKE ? OR name LIKE ?) ORDER BY o, x DESC");
		$getSequences->execute($cloneName,$getClone[5]);
		while(my @getSequences =  $getSequences->fetchrow_array())
		{
			$getSequences[5] = commify($getSequences[5]);
			$sequenceList .= ($getSequences[3] == 98) ? "<tr>
								<td><a onclick='closeDialog();openDialog(\"seqView.cgi?seqId=$getSequences[0]\")' title='View this sequence'>$seqType{$getSequences[3]}.$seqDir{$getSequences[6]}-end</a></td>
								<td><a href='download.cgi?seqId=$getSequences[0]' title='Download this sequence' target='hiddenFrame'>$getSequences[5]<span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-disk'></span></a></td>
								<td><a onclick='closeDialog();openDialog(\"libraryView.cgi?libraryId=$getSequences[4]\")'>LibraryId: $getSequences[4]</a>
								</td>
							</tr>" : "<tr>
								<td><a onclick='closeDialog();openDialog(\"seqView.cgi?seqId=$getSequences[0]\")' title='View this sequence'>$seqType{$getSequences[3]}</a> ($bacAssignType{$getSequences[7]}: <a title='mapped tags'>$getSequences[6]</a>)</td>
								<td><a href='download.cgi?seqId=$getSequences[0]' title='Download this sequence' target='hiddenFrame'>$getSequences[5]<span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-disk'></span></a></td>
								<td><a onclick='closeDialog();openDialog(\"jobView.cgi?jobId=$getSequences[4]\")'>JobId: $getSequences[4]</a>
								</td>
							</tr>";
		}
		$sequenceList .= "</tbody></table>";
		$getClone[6] = ($getSequences->rows > 0) ? $sequenceList : 'None.';
		
		undef $/;# enable slurp mode
		my $html = <DATA>;
		$html =~ s/\$cloneName/$cloneName/g;
		$html =~ s/\$originalClone/$originalClone/g;
		$html =~ s/\$relatedClone/$relatedClone/g;
		$html =~ s/\$fpcView/$fpcView/g;
		$html =~ s/\$sequenced/$getClone[6]/g;
		$html =~ s/\$\$/$$/g;

		print $html;
	}
}
else
{
	print 'No valid clone Found!';
	exit;
}

__DATA__
<div id="viewClone$$" name="viewClone$$">
	<table width='100%'>
	<tr><td style='text-align:right;white-space: nowrap;'><b>Clone Name</b></td><td>$cloneName</td></tr>
	$originalClone
	$relatedClone
	<tr><td style='text-align:right'><b>FPC</b></td><td>$fpcView</td></tr>
	<tr><td style='text-align:right'><b>Sequence</b></td><td>$sequenced</td></tr>
	</table>
</div>
<script>
buttonInit();
$( "#fpcCloneTabs$$" ).tabs();
$( "#sequenceList$$" ).dataTable({
	"scrollY": "200px",
	"scrollCollapse": true,
	"paging": false,
	"searching": false,
	"info": false
});
$('#dialog').dialog("option", "title", "View Clone $cloneName");
$( "#dialog" ).dialog( "option", "buttons", [ { text: "Print", click: function() {printDiv('viewClone$$'); } }, { text: "OK", click: function() {closeDialog(); } } ] );
</script>
