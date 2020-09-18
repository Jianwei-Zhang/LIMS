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

my $libraryId = param ('libraryId') || '';
my $library=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$library->execute($libraryId);
my @library = $library->fetchrow_array();
my $libraryDetails = decode_json $library[8];

my $plateId;
my $plateInLibrary=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'plate' AND z = ? ORDER BY o");
$plateInLibrary->execute($libraryId);
while (my @plateInLibrary = $plateInLibrary->fetchrow_array())
{
	$plateId->{$plateInLibrary[2]} = 1;
}
my $clones;
my $totalPlates = keys %{$plateId};
my $totalClones = 0;
my $totalPooled = 0;
my $totalSequenced = 0;
my $totalNotSequenced = 0;
for my $plate (sort keys %{$plateId})
{
	$clones .= "
		<table class='ui-widget-content ui-corner-all' style='page-break-inside: avoid;'><tr><th class='ui-state-highlight ui-corner-all'></th>";
	my @clones;
	my $originalName;
	my $hasSequence;
	my $poolClone;
	my $poolClones = $dbh->prepare("SELECT link.* FROM link,clones WHERE link.type LIKE 'poolClone' AND link.child = clones.name AND clones.libraryId = ? AND clones.plate LIKE ? ORDER BY link.child");
	$poolClones->execute($libraryId,$plate);
	while(my @poolClones = $poolClones->fetchrow_array())
	{
		$poolClone->{$poolClones[1]} = $poolClones[0];
	}
	my $pooledNumber = keys %$poolClone;
	$totalPooled += $pooledNumber;
	my $getClones = $dbh->prepare("SELECT * FROM clones WHERE libraryId = ? AND plate LIKE ? ORDER BY name");
	$getClones->execute($libraryId,$plate);
	my $sequencedNumber = 0;
	my $notSequencedNumber = 0;
	while(my @getClones = $getClones->fetchrow_array())
	{
		push @clones, $getClones[1];
		$originalName->{$getClones[1]} = $getClones[5];
		$hasSequence->{$getClones[1]} = $getClones[6];
		($getClones[6]) ? $sequencedNumber++ : $notSequencedNumber++;
	}
	$totalSequenced += $sequencedNumber;
	$totalNotSequenced += $notSequencedNumber;
	$totalClones = $totalSequenced + $totalNotSequenced;
	my %wellx = ( 96 => ['01' .. '12'], 384 =>['01' .. '24']);
	my %welly = ( 96 => ['A' .. 'H'], 384 =>['A' .. 'P']);
	my $rowspan = @{$welly{$libraryDetails->{'format'}}} + 1;
	for (sort @{$wellx{$libraryDetails->{'format'}}})
	{
		$clones .= "<th style='text-align:center;' class='ui-state-highlight ui-corner-all'>$_</th>";
	}
	$clones .= "<th style='text-align:left;padding: 10px;' rowspan='$rowspan'>
		Plate $plate<hr>
		<div style='position: relative;'>
		<span style='position: absolute;left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-triangle-1-nw' title='Pooled clone'></span><input type='checkbox' name='DEMO' disabled='disabled'>Pooled:$pooledNumber
		</div>
		<div style='position: relative;'>
		<span style='position: absolute;left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-check' title='Sequenced clone'></span><input type='checkbox' name='DEMO' disabled='disabled'>Sequenced:$sequencedNumber<br>
		Not sequenced:$notSequencedNumber
		</div>
		</th></tr>";
	
	for(sort @{$welly{$libraryDetails->{'format'}}})
	{
		my $y = $_;
		$clones .= "<tr><td style='text-align:center;' class='ui-state-highlight ui-corner-all'>$y</td>";
		for(sort @{$wellx{$libraryDetails->{'format'}}})
		{
			my $currentClone = shift @clones;
			my $pooledClone = (exists $poolClone->{$currentClone}) ? "<span style='position: absolute;left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-triangle-1-nw' title='Pooled clone'></span>" : "";
			my $sequencedClone = ($hasSequence->{$currentClone} > 0) ? "<span style='position: absolute;left: 2px;top: 2px;display:inline-block;' class='ui-icon ui-icon-check' title='Sequenced clone'></span>" : "";
			$clones .= "<td style='text-align:center;' class='ui-widget-content'><div style='position: relative;'>$pooledClone$sequencedClone<input type='checkbox' class='row$y col$_' id='cloneName$currentClone' name='cloneName' value='$currentClone' disabled='disabled'></div></td>";
		}
		$clones .= "</tr>";
	}
	$clones .= "</table><p></p>";
}

my $progressRate = ($totalClones > 0) ? int ($totalPooled*100/$totalClones) : '0';
my $successRate = ($totalPooled > 0) ? int ($totalSequenced*100/$totalPooled) : '0';
my $notSequencedRate = ($totalClones > 0) ? int ($totalNotSequenced*100/$totalClones) : '0';

my $total = "<h3>Library '$library[2]' Sequencing Status</h3>
		<p>
		<b>Total Plates:</b> $totalPlates<br>
		<b>Total Clones:</b> $totalClones<br>
		<b>Pooled Clones:</b> $totalPooled ($progressRate% of total)<br>
		<b>Sequenced Clones:</b> $totalSequenced ($successRate% of pooled)<br>
		<b>Not Sequenced:</b> $totalNotSequenced ($notSequencedRate% of total)<br>
		</p>";

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$total/$total/g;
$html =~ s/\$clones/$clones/g;
$html =~ s/\$\$/$$/g;
print header;
print $html;

__DATA__
<button style="float: right; margin-top: .3em; margin-right: .3em;" onclick="printDiv('seqLibraryReport$$')">Print</button>
<div id="seqLibraryReport$$"  name="seqLibraryReport$$">
$total
$clones
</div>
<script>
buttonInit();
$('#dialog').dialog("option", "title", "Sequencing Report");
$( "#dialog" ).dialog( "option", "buttons", [ { text: "Print", click: function() {printDiv('seqLibraryReport$$'); } }, { text: "OK", click: function() {closeDialog(); } } ] );
</script>