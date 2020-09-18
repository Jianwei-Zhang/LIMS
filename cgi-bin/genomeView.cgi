#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use JSON::XS;
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

my $genomeId = param ('genomeId') || '';
my $genome = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$genome->execute($genomeId);
my @genome=$genome->fetchrow_array();
my %yesOrNo = ( 0=>'No', 1=>'Yes' );
my $genomeStatus;
$genomeStatus->{0} = "Sequences not";
$genomeStatus->{-1} = "Sequences are being";
$genomeStatus->{1} = "$genome[3] sequences";

my $objectComponent;
$objectComponent->{0} = "Unknown";
$objectComponent->{1} = "Chr-Seq";
$objectComponent->{2} = "Ctg-Seq";

my $relatedLibraries = 'None.';
if ($genome[6])
{
	my $relatedLibrary=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$relatedLibrary->execute($genome[6]);
	my @relatedLibrary = $relatedLibrary->fetchrow_array();
	$relatedLibraries = "<a onclick='closeDialog();openDialog(\"libraryView.cgi?libraryId=$genome[6]\")'>$relatedLibrary[2]</a> ";
}

my $sequences = "
	<table id='genome$$' class='display' style='width: 100%;'>
		<thead>
			<tr>
				<th><b>Seq Name</b></th>
				<th><b title='numbers in this column may be inaccurate and are editable.'>Chr</b></th>
				<th><b>Length (bp)</b></th>
				<th><b>Gap#</b></th>
				<th><b>Action</b></th>
			</tr>
		</thead>
		<tbody>";
my @sequenceLength = ();
my $totalLength = 0;
my $totalGaps = 0;
my $gapCount;
my $getSequences = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 99 AND x = ?");
$getSequences->execute($genomeId);
while (my @getSequences =  $getSequences->fetchrow_array())
{
	push @sequenceLength, $getSequences[5];
	$totalLength += $getSequences[5];
	$getSequences[5] = commify($getSequences[5]);
	my $sequenceDetails = decode_json $getSequences[8];
	$sequenceDetails->{'id'} = '' unless (exists $sequenceDetails->{'id'});
	$sequenceDetails->{'description'} = '' unless (exists $sequenceDetails->{'description'});
	my $descriptionLength = length ($sequenceDetails->{'description'});
	if ($descriptionLength > 35)
	{
		$sequenceDetails->{'description'} = "<a title='$sequenceDetails->{'description'}'>". substr($sequenceDetails->{'description'},0,25). "...". substr($sequenceDetails->{'description'},-5). "</a>";
	}
	$sequenceDetails->{'gapList'} = '' unless (exists $sequenceDetails->{'gapList'});
	my $gaps = 0;
	my @gaps = split (",",$sequenceDetails->{'gapList'});
	if(@gaps)
	{
		$gaps = @gaps;
	}
	$totalGaps += $gaps;
	if(exists $gapCount->{$gaps})
	{
		$gapCount->{$gaps}++
	}
	else
	{
		$gapCount->{$gaps} = 1;
	}
	$sequences .= "<tr>
	<td><a onclick='closeDialog();openDialog(\"seqView.cgi?seqId=$getSequences[0]\")' title='View $getSequences[2]'>$getSequences[2]</a><br><sup>$sequenceDetails->{'description'}</sup></td>
	<td><div id='seqChr$getSequences[0]' style='position: relative;'><a id='seqChr$getSequences[0]$$' onmouseover='editIconShow(\"seqChr$getSequences[0]$$\")' onmouseout='editIconHide(\"seqChr$getSequences[0]$$\")' onclick='loaddiv(\"seqChr$getSequences[0]\",\"seqChrEdit.cgi?seqId=$getSequences[0]\")' title='Edit this chromosome number'>$getSequences[6]</a></div></td>
	<td>$getSequences[5]</td>
	<td>$gaps</td>
	<td><a onclick='closeDialog();closeViewer();openViewer(\"seqBrowse.cgi?seqId=$getSequences[0]\")'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-image' title='Browse'></span></a></td>
	</tr>";
}

@sequenceLength = sort {$b <=> $a} @sequenceLength;
my $medianLength = 0;
if($#sequenceLength % 2 == 1)
{
	my $median = int ($#sequenceLength/2);
	$medianLength = ($sequenceLength[$median]+$sequenceLength[$median+1])/2;
}
else
{
	my $median = $#sequenceLength/2;
	$medianLength = $sequenceLength[$median];
}

my $n50Length = 0;
my $subtotal = 0;
my $assemblyCtgLengthCount;
foreach (@sequenceLength)
{
	$subtotal += $_;
	if($subtotal > $totalLength/2 && $n50Length == 0)
	{
		$n50Length = $_;
	}
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$genomeId/$genomeId/g;
$html =~ s/\$genomeName/$genome[2]/g;
if($genome[3] > 1)
{
	$html =~ s/\$sequenceNumber/$genome[3] Sequences/g;
}
else
{
	$html =~ s/\$sequenceNumber/$genome[3] Sequence/g;
}

$totalLength = commify($totalLength);
$html =~ s/\$totalLength/$totalLength/g;
$n50Length = commify($n50Length);
$html =~ s/\$n50Length/$n50Length/g;
$medianLength = commify($medianLength);
$html =~ s/\$medianLength/$medianLength/g;
$totalGaps = commify($totalGaps);

my $gapStats = "<table id='gap$$' class='display' style='width: 100%;'>
		<thead>
			<tr>
				<th><b>Gap#</b></th>
				<th><b>Count</b></th>
			</tr>
		</thead>
		<tbody>";
foreach (sort {$a <=> $b} keys %$gapCount)
{
	$gapStats .= "<tr><td>Gap-$_</td><td>$gapCount->{$_}</td></tr>";
}
$gapStats .= "</tbody>
	<tfoot>
		<tr>
			<th><b>Total</b></th>
			<th>$totalGaps</th>
		</tr>
	</tfoot>
</table>\n";

$html =~ s/\$gapStats/$gapStats/g;

$sequences .= "</tbody>
	<tfoot>
		<tr>
			<th><b>Total</b></th>
			<th></th>
			<th>$totalLength</th>
			<th>$totalGaps</th>
			<th></th>
		</tr>
	</tfoot>
</table>\n";

$html =~ s/\$genomeForGPM/$yesOrNo{$genome[4]}/g;
$html =~ s/\$genomeAsReference/$yesOrNo{$genome[5]}/g;
my $agpAvailable = "";
my $agpList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'agp' AND x = ? ORDER BY o");
$agpList->execute($genome[0]);
while (my @agpList = $agpList->fetchrow_array())
{
	$agpAvailable .= ($agpAvailable) ? "<li><a href='download.cgi?agpId=$agpList[0]' target='hiddenFrame' title='$objectComponent->{$agpList[5]} AGP v.$agpList[3]'>$agpList[2]</a></li>" : "Available AGP:<ul><li><a href='download.cgi?agpId=$agpList[0]' target='hiddenFrame' title='$objectComponent->{$agpList[5]} AGP v.$agpList[3]'>$agpList[2]</a></li>";
}
$agpAvailable .= ($agpAvailable) ? '</ul>' : '';

$html =~ s/\$agpAvailable/$agpAvailable/g;
$html =~ s/\$relatedLibraries/$relatedLibraries/g;
$html =~ s/\$genomeStatus/$genomeStatus->{$genome[7]}/g;
$genome[8] = escapeHTML($genome[8]);
$genome[8] =~ s/\n/<br>/g;
$html =~ s/\$genomeDescription/$genome[8]/g;
$html =~ s/\$genomeCreator/$genome[9]/g;
$html =~ s/\$genomeCreationDate/$genome[10]/g;
$html =~ s/\$sequences/$sequences/g;
$html =~ s/\$\$/$$/g;
print header;
print $html;

__DATA__
<table>
	<tr><td style='text-align:right;white-space: nowrap;'><b>Genome</b></td><td>$genomeName<br>$genomeStatus <sup class='ui-state-disabled'>loaded by $genomeCreator on $genomeCreationDate</sup></td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'>For Reassembly?</td><td>$genomeForGPM. $agpAvailable</td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'>As Reference?</td><td>$genomeAsReference.</td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'>Link to library:</td><td>$relatedLibraries</td></tr>
</table>
<div id='genomeTabs$$'>
	<ul>
		<li><a href='#genomeTabs-sequences$$'>Sequences</a></li>
		<li><a href='#genomeTabs-stats$$'>Stats</a></li>
		<li><a href='#genomeTabs-description$$'>Description</a></li>
	</ul>
	<div id='genomeTabs-sequences$$'>$sequences</div>
	<div id='genomeTabs-stats$$'><b>$sequenceNumber</b><br>Total: $totalLength bp<br>Median: $medianLength bp<br>N50: $n50Length bp<br><hr>$gapStats</div>
	<div id='genomeTabs-description$$'>$genomeDescription</div>
</div>
<script>
$('#dialog').dialog("option", "title", "View Genome");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Alignment as Query", click: function() { closeDialog();openDialog("alignmentForm.cgi?queryId=$genomeId"); } }, { text: "Alignment as Reference", click: function() { closeDialog();openDialog("alignmentForm.cgi?subjectId=$genomeId"); } }, { text: "Edit", click: function() { closeDialog();openDialog("genomeEdit.cgi?genomeId=$genomeId"); } }, { text: "OK", click: function() {closeDialog(); } } ] );
$( "#genome$$" ).dataTable({
	"scrollY": "400px",
	"scrollCollapse": false,
	"paging": false,
	"columnDefs": [
    { "orderable": false, "targets": 4 }
  ]
});
$( "#gap$$" ).dataTable({
	"dom": 'lfrtip',
	"scrollY": "400px",
	"scrollCollapse": true,
	"paging": false,
	"searching": false,
	"info": false
});
$( "#genomeTabs$$" ).tabs();
</script>