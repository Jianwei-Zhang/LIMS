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

my %status = (
	0=>'Na',
	1=>'Status One',
	2=>'Status Two'
	);
my $paclibId = param ('paclibId') || '';
my $paclib = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$paclib->execute($paclibId);
my @paclib=$paclib->fetchrow_array();
my $paclibDetails = decode_json $paclib[8];
$paclibDetails->{'libraryDate'} = '' unless (exists $paclibDetails->{'libraryDate'});
$paclibDetails->{'shearingInput'} = '' unless (exists $paclibDetails->{'shearingInput'});
$paclibDetails->{'shearingRpm'} = '' unless (exists $paclibDetails->{'shearingRpm'});
$paclibDetails->{'shearingOutput'} = '' unless (exists $paclibDetails->{'shearingOutput'});
$paclibDetails->{'shearingBeadsSteps'} = '' unless (exists $paclibDetails->{'shearingBeadsSteps'});
$paclibDetails->{'shearingOperator'} = '' unless (exists $paclibDetails->{'shearingOperator'});
$paclibDetails->{'bluepippinInput'} = '' unless (exists $paclibDetails->{'bluepippinInput'});
$paclibDetails->{'bluepippinSize'} = '' unless (exists $paclibDetails->{'bluepippinSize'});
$paclibDetails->{'bluepippinOutput'} = '' unless (exists $paclibDetails->{'bluepippinOutput'});
$paclibDetails->{'bluepippinConcentration'} = '' unless (exists $paclibDetails->{'bluepippinConcentration'});
$paclibDetails->{'bluepippinOperator'} = '' unless (exists $paclibDetails->{'bluepippinOperator'});
$paclibDetails->{'description'} = '' unless (exists $paclibDetails->{'description'});

my $paclibStatus = '';
foreach (sort {$a <=> $b} keys %status)
{
	$paclibStatus .= ($paclib[4] == $_) ? "<option value='$_' selected>$status{$_}</option>" : "<option value='$_'>$status{$_}</option>";
}
my $noun = 'word';
my $paclibToSample=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$paclibToSample->execute($paclib[6]);
my @paclibToSample = $paclibToSample->fetchrow_array();
my $sampleToService=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$sampleToService->execute($paclibToSample[6]);
my @sampleToService = $sampleToService->fetchrow_array();
my $serviceToProject=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$serviceToProject->execute($sampleToService[6]);
my @serviceToProject = $serviceToProject->fetchrow_array();

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$paclibId/$paclibId/g;
$html =~ s/\$sampleId/$paclib[6]/g;
$html =~ s/\$project/$serviceToProject[2]/g;
$html =~ s/\$serviceName/$sampleToService[2]/g;
$html =~ s/\$sampleName/$paclibToSample[2]/g;
$html =~ s/\$paclibName/$paclib[2]/g;
$html =~ s/\$paclibStatus/$paclibStatus/g;
$html =~ s/\$paclibBarcode/$paclib[7]/g;
$html =~ s/\$paclibDate/$paclibDetails->{'libraryDate'}/g;
$html =~ s/\$shearingInput/$paclibDetails->{'shearingInput'}/g;
$html =~ s/\$shearingRpm/$paclibDetails->{'shearingRpm'}/g;
$html =~ s/\$shearingOutput/$paclibDetails->{'shearingOutput'}/g;
$html =~ s/\$shearingBeadsSteps/$paclibDetails->{'shearingBeadsSteps'}/g;
$html =~ s/\$shearingOperator/$paclibDetails->{'shearingOperator'}/g;
$html =~ s/\$bluepippinInput/$paclibDetails->{'bluepippinInput'}/g;
$html =~ s/\$bluepippinSize/$paclibDetails->{'bluepippinSize'}/g;
$html =~ s/\$bluepippinOutput/$paclibDetails->{'bluepippinOutput'}/g;
$html =~ s/\$bluepippinConcentration/$paclibDetails->{'bluepippinConcentration'}/g;
$html =~ s/\$bluepippinOperator/$paclibDetails->{'bluepippinOperator'}/g;
$html =~ s/\$paclibDescription/$paclibDetails->{'description'}/g;
$html =~ s/\$paclibCreator/$paclib[9]/g;
$html =~ s/\$paclibCreationDate/$paclib[10]/g;
$html =~ s/\$noun/$noun/g;

print header;
print $html;

__DATA__
<form id="editPaclib" name="editPaclib" action="paclibSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="paclibId" id="editPaclibId" type="hidden" value="$paclibId" />
	<table>
	<tr><td style='text-align:right;white-space: nowrap;'><b>For</b></td><td><a title='Project'>$project</a> > <a title='Service'>$serviceName</a> > <a title='Sample'>$sampleName</a><input name="sampleId" id="editPaclibSampleId" type="hidden" value="$sampleId" /></td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="editPaclibName"><b>Paclib Name</b></label></td><td><input class='ui-widget-content ui-corner-all' name="name" id="editPaclibName" size="10" type="text" maxlength="32" value="$paclibName" /><img alt='$paclibBarcode' src='barcode.cgi?code=$paclibBarcode'/><sup class='ui-state-disabled'>by $paclibCreator on $paclibCreationDate</sup></td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="editPaclibDate"><b>Library Date</b></label></td><td><input class='ui-widget-content ui-corner-all' name="libraryDate" id="editPaclibDate" size="10" type="text" maxlength="10" placeholder="YYYY-MM-DD" value="$paclibDate" /></td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><b>Shearing</b></td><td>
		<label for="editPaclibShearingInput">DNA input</label><input class='ui-widget-content ui-corner-all' name="shearingInput" id="editPaclibShearingInput" size="10" type="text" maxlength="10" placeholder="" value="$shearingInput" />ug<br>
		<label for="editPaclibShearingRpm">RPM</label><input class='ui-widget-content ui-corner-all' name="shearingRpm" id="editPaclibShearingRpm" size="10" type="text" maxlength="10" placeholder="" value="$shearingRpm" /><br>
		<label for="editPaclibShearingOutput">DNA output</label><input class='ui-widget-content ui-corner-all' name="shearingOutput" id="editPaclibShearingOutput" size="10" type="text" maxlength="10" placeholder="" value="$shearingOutput" />ug<br>
		<label for="editPaclibShearingBeadsSteps">PB beads steps</label><input class='ui-widget-content ui-corner-all' name="shearingBeadsSteps" id="editPaclibShearingBeadsSteps" size="10" type="text" maxlength="10" placeholder="" value="$shearingBeadsSteps" /><br>
		<label for="editPaclibShearingOperator">Operator</label><input class='ui-widget-content ui-corner-all' name="shearingOperator" id="editPaclibShearingOperator" size="10" type="text" maxlength="10" placeholder="Operator" value="$shearingOperator" /><br>
		</td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><b>BluePippin</b></td><td>
		<label for="editPaclibBluepippinInput">DNA input</label><input class='ui-widget-content ui-corner-all' name="bluepippinInput" id="editPaclibBluepippinInput" size="10" type="text" maxlength="10" placeholder="" value="$bluepippinInput" />ug<br>
		<label for="editPaclibBluepippinSize">Size selection</label><input class='ui-widget-content ui-corner-all' name="bluepippinSize" id="editPaclibBluepippinSize" size="10" type="text" maxlength="10" placeholder="" value="$bluepippinSize" />kb<br>
		<label for="editPaclibBluepippinOutput">DNA output</label><input class='ui-widget-content ui-corner-all' name="bluepippinOutput" id="editPaclibBluepippinOutput" size="10" type="text" maxlength="10" placeholder="" value="$bluepippinOutput" />ug<br>
		<label for="editPaclibBluepippinConcentration">DNA concentration</label><input class='ui-widget-content ui-corner-all' name="bluepippinConcentration" id="editPaclibBluepippinConcentration" size="10" type="text" maxlength="10" placeholder="" value="$bluepippinConcentration" />ng/uL<br>
		<label for="editPaclibBluepippinOperator">Operator</label><input class='ui-widget-content ui-corner-all' name="bluepippinOperator" id="editPaclibBluepippinOperator" size="10" type="text" maxlength="10" placeholder="Operator" value="$bluepippinOperator" /><br>
		</td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="editPaclibStatus"><b>Status</b></label></td><td><select class='ui-widget-content ui-corner-all' name="status" id="editPaclibStatus">$paclibStatus</select></td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="editPaclibDescription"><b>Description</b></label><br><sub id='editPaclibDescription_count' style='display:none'></sub></td><td><textarea class='ui-widget-content ui-corner-all word_count' name="description" id="editPaclibDescription" cols="50" rows="5">$paclibDescription</textarea></td></tr>
	</table>
</form>
<script>
$('#dialog').dialog("option", "title", "Edit Paclib");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('editPaclib'); } }, { text: "Delete", click: function() { deleteItem($paclibId); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
$( "#editPaclibDate" ).datepicker({dateFormat:"yy-mm-dd"});
wordCount('$noun');
</script>