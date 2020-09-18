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

my $sequences = "
	<table id='genome$$' class='display'>
		<thead>
			<tr>
				<th style='white-space: nowrap;'><b>Seq Name</b></th>
				<th style='white-space: nowrap;'><b title='numbers in this column may be inaccurate and are editable.'>Chr</b></th>
				<th style='white-space: nowrap;'><b>Length (bp)</b></th>
				<th style='white-space: nowrap;'><b>Gap#</b></th>
				<th style='white-space: nowrap;'><b>Action</b></th>
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
	$sequences .= "<tr><td style='white-space: nowrap;'><a onclick='closeDialog();openDialog(\"seqView.cgi?seqId=$getSequences[0]\")' title='View $getSequences[2]'>$getSequences[2]</a><br><sup>$sequenceDetails->{'description'}</sup></td><td><div id='seqChr$getSequences[0]' style='position: relative;'><a id='seqChr$getSequences[0]$$' onmouseover='editIconShow(\"seqChr$getSequences[0]$$\")' onmouseout='editIconHide(\"seqChr$getSequences[0]$$\")' onclick='loaddiv(\"seqChr$getSequences[0]\",\"seqChrEdit.cgi?seqId=$getSequences[0]\")' title='Edit this chromosome number'>$getSequences[6]</a></div></td><td>$getSequences[5]</td><td>$gaps</td><td><a style='float: right;' onclick='closeDialog();closeViewer();openViewer(\"seqBrowse.cgi?seqId=$getSequences[0]\")'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-image' title='Browse'></span></a></td></tr>";
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

my $gapStats = "<table id='gap$$' class='display'>
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
$agpAvailable .= ($agpAvailable) ? '</ul>' : '<br>';

$html =~ s/\$agpAvailable/$agpAvailable/g;
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
	<tr><td style='text-align:left;white-space: nowrap;'><h3><a onclick='openDialog("genomeView.cgi?genomeId=$genomeId")' title='View'>$genomeName</a> <sub class='ui-state-disabled'>$sequenceNumber</sub></h3></td><td></td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'>$sequences</td><td>For Reassembly? $genomeForGPM. $agpAvailable 
		As Reference? $genomeAsReference. <br>
		Total: $totalLength bp;
		Median: $medianLength bp;
		N50: $n50Length bp;<br>
		$gapStats<br><b>Description:</b> $genomeDescription</td></tr>
</table>
<script>
$( "#gap$$" ).dataTable({
	"dom": 'lfrtip',
	"scrollY": "300px",
	"scrollCollapse": true,
	"paging": false,
	"searching": false,
	"info": false
});
$( "#genome$$" ).dataTable({
	"scrollY": "400px",
	"scrollCollapse": true,
	"paging": false
});
</script>