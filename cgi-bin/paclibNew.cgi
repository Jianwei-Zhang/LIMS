#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use JSON::XS;
use DBI;
use POSIX qw(strftime);
use lib "lib/";
use lib "lib/pangu";
use pangu;
use userCookie;

my $userCookie = new userCookie;
my $userId = (cookie('cid')) ? $userCookie->checkCookie(cookie('cid')) : 0;
exit if (!$userId);

my $commoncfg = readConfig("main.conf");
my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});

my %status = (
	0=>'Na',
	1=>'Status One',
	2=>'Status Two'
	);
my $paclibStatus = '';
foreach (sort {$a <=> $b} keys %status)
{
	$paclibStatus .= "<option value='$_'>$status{$_}</option>";
}

my $sampleId = param ('sampleId') || '';
my $noun = 'word';

my $sample=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$sample->execute($sampleId);
my @sample = $sample->fetchrow_array();
my $sampleToService=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$sampleToService->execute($sample[6]);
my @sampleToService = $sampleToService->fetchrow_array();
my $serviceToProject=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$serviceToProject->execute($sampleToService[6]);
my @serviceToProject = $serviceToProject->fetchrow_array();

my $datetimeString = strftime "%y %Y %m %d %H:%M:%S", localtime;
my ($tdYear,$year,$month,$dayofmonth,$time)=split(/\s+/,$datetimeString);
my $autoPaclibDate = "$year-$month-$dayofmonth";

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$sampleId/$sampleId/g;
$html =~ s/\$sampleName/$sample[2]/g;
$html =~ s/\$serviceName/$sampleToService[2]/g;
$html =~ s/\$projectName/$serviceToProject[2]/g;
$html =~ s/\$autoPaclibDate/$autoPaclibDate/g;
$html =~ s/\$paclibStatus/$paclibStatus/g;
$html =~ s/\$noun/$noun/g;

print header;
print $html;

__DATA__
<form id="newPaclib" name="newPaclib" action="paclibSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<table>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="newPaclibSampleId"><b>For<b></label></td><td><a title='Project'>$projectName</a> > <a title='Service'>$serviceName</a> > <a title='Sample'>$sampleName</a> <input name="sampleId" id="newPaclibSampleId" type="hidden" value="$sampleId" /></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="newPaclibName"><b>Paclib Name</b></label></td><td><input class='ui-widget-content ui-corner-all' name="name" id="newPaclibName" size="10" type="text" maxlength="32" value="" /></td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="newPaclibDate"><b>Library Date</b></label></td><td><input class='ui-widget-content ui-corner-all' name="libraryDate" id="newPaclibDate" size="10" type="text" maxlength="10" placeholder="YYYY-MM-DD" value="$autoPaclibDate" /></td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><b>Shearing</b></td><td>
		<label for="newPaclibShearingInput">DNA input</label><input class='ui-widget-content ui-corner-all' name="shearingInput" id="newPaclibShearingInput" size="10" type="text" maxlength="10" placeholder="" value="" />ug<br>
		<label for="newPaclibShearingRpm">RPM</label><input class='ui-widget-content ui-corner-all' name="shearingRpm" id="newPaclibShearingRpm" size="10" type="text" maxlength="10" placeholder="" value="" /><br>
		<label for="newPaclibShearingOutput">DNA output</label><input class='ui-widget-content ui-corner-all' name="shearingOutput" id="newPaclibShearingOutput" size="10" type="text" maxlength="10" placeholder="" value="" />ug<br>
		<label for="newPaclibShearingBeadsSteps">PB beads steps</label><input class='ui-widget-content ui-corner-all' name="shearingBeadsSteps" id="newPaclibShearingBeadsSteps" size="10" type="text" maxlength="10" placeholder="" value="" /><br>
		<label for="newPaclibShearingOperator">Operator</label><input class='ui-widget-content ui-corner-all' name="shearingOperator" id="newPaclibShearingOperator" size="10" type="text" maxlength="10" placeholder="Operator" value="" /><br>
		</td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><b>BluePippin</b></td><td>
		<label for="newPaclibBluepippinInput">DNA input</label><input class='ui-widget-content ui-corner-all' name="bluepippinInput" id="newPaclibBluepippinInput" size="10" type="text" maxlength="10" placeholder="" value="" />ug<br>
		<label for="newPaclibBluepippinSize">Size selection</label><input class='ui-widget-content ui-corner-all' name="bluepippinSize" id="newPaclibBluepippinSize" size="10" type="text" maxlength="10" placeholder="" value="" />kb<br>
		<label for="newPaclibBluepippinOutput">DNA output</label><input class='ui-widget-content ui-corner-all' name="bluepippinOutput" id="newPaclibBluepippinOutput" size="10" type="text" maxlength="10" placeholder="" value="" />ug<br>
		<label for="newPaclibBluepippinConcentration">DNA concentration</label><input class='ui-widget-content ui-corner-all' name="bluepippinConcentration" id="newPaclibBluepippinConcentration" size="10" type="text" maxlength="10" placeholder="" value="" />ng/uL<br>
		<label for="newPaclibBluepippinOperator">Operator</label><input class='ui-widget-content ui-corner-all' name="bluepippinOperator" id="newPaclibBluepippinOperator" size="10" type="text" maxlength="10" placeholder="Operator" value="" /><br>
		</td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="newPaclibStatus"><b>Status</b></label></td><td><select class='ui-widget-content ui-corner-all' name="status" id="newPaclibStatus">$paclibStatus</select></td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="newPaclibDescription"><b>Description</b></label><br><sub id='newPaclibDescription_count' style='display:none'></sub></td><td><textarea class='ui-widget-content ui-corner-all word_count' name="description" id="newPaclibDescription" cols="50" rows="5"></textarea></td></tr>
	</table>
</form>
<script>

$('#dialog').dialog("option", "title", "New Paclib");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('newPaclib'); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
$( "#newPaclibDate" ).datepicker({dateFormat:"yy-mm-dd"});
wordCount('$noun');
</script>