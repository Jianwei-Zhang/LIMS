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
	98=>'BES',
	99=>'Genome'
	);

my $libraryId = param ('libraryId') || '';
my $plate = param ('plate') || '';
my $button = "<button style='float: right; margin-top: .3em; margin-right: .3em;' onclick='openDialogForm(\"poolNew.cgi\", \"newPoolClone$$\")'>New pool</button>
	<button style='float: right; margin-top: .3em; margin-right: .3em;' onclick='printDiv(\"newPoolClone$$\")'>Print</button>";
my $clones = "
	<form id='newPoolClone$$' name='newPoolClone$$' style='margin-top: .3em; margin-right: .3em;'>
	<input name='libraryId' id='seqLibrary$$' type='hidden' value='$libraryId' />
	<table class='ui-widget-content ui-corner-all'><tr><th class='ui-state-error ui-corner-all'><a onClick='checkAll(\"cloneName\")' title='Check All'><span class='ui-icon ui-icon-circle-check'></span></a></th>";
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

my $getClones = $dbh->prepare("SELECT * FROM clones WHERE libraryId = ? AND plate LIKE ? ORDER BY name");
$getClones->execute($libraryId,$plate);
if($getClones->rows < 1)
{
	$clones = "No Clones Found.<a onclick='openDialog(\"cloneNew.cgi?libraryId=$libraryId\")'>Generate Clone List</a>";
}
else
{
	my $sequencedNumber = 0;
	my $notSequencedNumber = 0;
	while(my @getClones = $getClones->fetchrow_array())
	{
		push @clones, $getClones[1];
		$originalName->{$getClones[1]} = $getClones[5];
		$hasSequence->{$getClones[1]} = $getClones[6];
		($getClones[6]) ? $sequencedNumber++ : $notSequencedNumber++;
	}
	my %wellx = ( 96 => ['01' .. '12'], 384 =>['01' .. '24']);
	my %welly = ( 96 => ['A' .. 'H'], 384 =>['A' .. 'P']);
	my $sourceLibrary=$dbh->prepare("SELECT * FROM matrix WHERE id = $libraryId");
	$sourceLibrary->execute();
	my @sourceLibrary = $sourceLibrary->fetchrow_array();
	my $sourceLibraryDetails = decode_json $sourceLibrary[8];
    my $rowspan = @{$welly{$sourceLibraryDetails->{'format'}}} + 2;
	my $tableBottom;
	for (sort @{$wellx{$sourceLibraryDetails->{'format'}}})
	{
		$clones .= "<th style='text-align:center;' class='ui-state-error ui-corner-all'><a onclick='checkClass(\"col$_\")' title='Check Column $_'>$_</a></th>";
		$tableBottom .= "<th style='text-align:center;' class='ui-widget-content ui-state-disabled ui-corner-all'><a onclick='uncheckClass(\"col$_\")' title='UnCheck Column $_'>$_</a></th>";
	}
	$clones .= "<th></th><th style='text-align:left;padding: 10px;' rowspan='$rowspan'>
		Plate $plate<hr>
		<div style='position: relative;'>
		<span style='position: absolute;left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-triangle-1-nw' title='Pooled clone'></span><input type='checkbox' name='DEMO' disabled='disabled'>Pooled:$pooledNumber
		</div>
		<div style='position: relative;'>
		<span style='position: absolute;left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-check' title='Sequenced clone'></span><input type='checkbox' name='DEMO' disabled='disabled'>Sequenced:$sequencedNumber<br>
		Not sequenced:$notSequencedNumber
		</div>
		
		</th></tr>";
	
	for(sort @{$welly{$sourceLibraryDetails->{'format'}}})
	{
		my $y = $_;
		$clones .= "<tr><td style='text-align:center;' class='ui-state-error ui-corner-all'><a onclick='checkClass(\"row$y\")' title='Check Row $y'>$y</a></td>";
		for(sort @{$wellx{$sourceLibraryDetails->{'format'}}})
		{
			my $currentClone = shift @clones;
			my $disabledWell = ($originalName->{$currentClone} eq "EMPTY") ? " title='Empty' disabled='disabled'" : " title='$originalName->{$currentClone}'" ;
			my $pooledClone = (exists $poolClone->{$currentClone}) ? "<label for='cloneName$currentClone'><span style='position: absolute;left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-triangle-1-nw' title='Pooled clone'></span></label>" : "";
			my $sequencedClone = ($hasSequence->{$currentClone} > 0) ? "<label for='cloneName$currentClone'><span style='position: absolute;left: 2px;top: 2px;display:inline-block;' class='ui-icon ui-icon-check' title='$originalName->{$currentClone}:$seqType{$hasSequence->{$currentClone}}'></span></label>" : "";
			$clones .= ($hasSequence->{$currentClone} > 0) ? "<td style='text-align:center;' class='ui-state-highlight'><div style='position: relative;'>$pooledClone$sequencedClone<input type='checkbox' class='row$y col$_' id='cloneName$currentClone' name='cloneName' value='$currentClone'$disabledWell></div></td>" :
					"<td style='text-align:center;' class='ui-widget-content'><div style='position: relative;'>$pooledClone$sequencedClone<input type='checkbox' class='row$y col$_' id='cloneName$currentClone' name='cloneName' value='$currentClone'$disabledWell></div></td>";
		}
		$clones .= "<td style='text-align:center;' class='ui-widget-content ui-state-disabled ui-corner-all'><a onclick='uncheckClass(\"row$y\")' title='UnCheck Row $y'>$y</a></td></tr>";
	}
	$clones .= "<tr><th></th>$tableBottom<th class='ui-widget-content ui-state-disabled ui-corner-all'><a onClick='uncheckAll(\"cloneName\")' title='Uncheck All'><span class='ui-icon ui-icon-circle-close'></span></a></th></tr>";
	$clones .= "</table></form>";
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$button/$button/g;
$html =~ s/\$clones/$clones/g;

print header(-cookie=>cookie(-name=>'plate',-value=>$libraryId.$plate));;
print $html;

__DATA__
$button
$clones
<script>
buttonInit();
loadingHide();
</script>