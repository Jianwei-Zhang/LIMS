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

my %purpose = (
	0=>'General Item',
	1=>'Sample for PacBio Library'
	);

my %type = (
	0=>'Other:',
	1=>'Genomic DNA',
	2=>'Total RNA',
	3=>'BAC clones, as bacteria',
	4=>'BAC plasmids',
	5=>'PCR amplicons',
	6=>'cDNA',
	7=>'Tissue(s):'
	);

my %status = (
	0=>'Na',
	1=>'Status One',
	2=>'Status Two'
	);
my $noun = 'word';
my $sampleId = param ('sampleId') || '';
my $sample = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$sample->execute($sampleId);
my @sample=$sample->fetchrow_array();

my $sampleDetails = decode_json $sample[8];
$sampleDetails->{'sampleTypeOther'} = '' unless (exists $sampleDetails->{'sampleTypeOther'});
$sampleDetails->{'description'} = '' unless (exists $sampleDetails->{'description'});

my $samplePurpose = '';
foreach (sort {$a <=> $b} keys %purpose)
{
	$samplePurpose .= ($sample[3] eq $_) ? "<option value='$_' selected>$purpose{$_}</option>" : "<option value='$_'>$purpose{$_}</option>";
}

my $sampleTypeList = '';
foreach (sort {$a <=> $b} keys %type)
{
	next unless ($_);
	$sampleTypeList .= ($sample[4] == $_) ? "<option value='$_' selected>$_. $type{$_}</option>" : "<option value='$_'>$_. $type{$_}</option>";
}
$sampleTypeList .= ($sample[4] == 0) ? "<option value='0' selected>$type{0}</option>" :  "<option value='0'>$type{0}</option>";

my $sampleStatus = '';
foreach (sort {$a <=> $b} keys %status)
{
	$sampleStatus .= ($sample[5] eq $_) ? "<option value='$_' selected>$status{$_}</option>" : "<option value='$_'>$status{$_}</option>";
}
my $sampleToService=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$sampleToService->execute($sample[6]);
my @sampleToService = $sampleToService->fetchrow_array();
my $serviceToProject=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$serviceToProject->execute($sampleToService[6]);
my @serviceToProject = $serviceToProject->fetchrow_array();

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$sampleTypeOther/$sampleDetails->{'sampleTypeOther'}/g;
$html =~ s/\$sampleType/$sampleTypeList/g;
$html =~ s/\$sampleId/$sampleId/g;
$html =~ s/\$projectName/$serviceToProject[2]/g;
$html =~ s/\$serviceId/$sample[6]/g;
$html =~ s/\$serviceName/$sampleToService[2]/g;
$html =~ s/\$sampleName/$sample[2]/g;
$html =~ s/\$samplePurpose/$samplePurpose/g;
$html =~ s/\$sampleStatus/$sampleStatus/g;
$html =~ s/\$sampleBarcode/$sample[7]/g;
$html =~ s/\$sampleDescription/$sampleDetails->{'description'}/g;
$html =~ s/\$sampleCreator/$sample[9]/g;
$html =~ s/\$sampleCreationDate/$sample[10]/g;
$html =~ s/\$noun/$noun/g;

print header;
print $html;

__DATA__
<form id="editSample" name="editSample" action="sampleSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="sampleId" id="editSampleId" type="hidden" value="$sampleId" />
	<table>
	<tr><td style='text-align:right;white-space: nowrap;'>For</td><td><a title='Project'>$projectName</a> > <a title='Service'>$serviceName</a><input name="serviceId" id="editSampleServiceId" type="hidden" value="$serviceId" /></td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="editSampleName"><b>Sample Name</b></label></td><td><input class='ui-widget-content ui-corner-all' name="name" id="editSampleName" size="10" type="text" maxlength="32" value="$sampleName" /><img alt='$sampleBarcode' src='barcode.cgi?code=$sampleBarcode'/><sup class='ui-state-disabled'>by $sampleCreator on $sampleCreationDate</sup></td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="editSamplePurpose"><b>Purpose</b></label></td><td><select class='ui-widget-content ui-corner-all' name="purpose" id="editSamplePurpose">$samplePurpose</select></td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="editSampleType"><b>Type</b></label></td>
		<td>
			<select class='ui-widget-content ui-corner-all' name="type" id="editSampleType" onchange="hideShowTypeOther();">$sampleType</select>
			<br><textarea class='ui-widget-content ui-corner-all' name="sampleTypeOther" id="editSampleTypeOther" cols="40" rows="2" placeholder="">$sampleTypeOther</textarea>
		</td>
	</tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="editSampleStatus"><b>Status</b></label></td><td><select class='ui-widget-content ui-corner-all' name="status" id="editSampleStatus">$sampleStatus</select></td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="editSampleDescription"><b>Description</b></label><br><sub id='editSampleDescription_count' style='display:none'></sub></td><td><textarea class='ui-widget-content ui-corner-all word_count' name="description" id="editSampleDescription" cols="50" rows="10">$sampleDescription</textarea></td></tr>
	</table>
</form>
<script>
hideShowTypeOther();
$('#dialog').dialog("option", "title", "Edit Sample");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('editSample'); } }, { text: "Delete", click: function() { deleteItem($sampleId); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
wordCount('$noun');
function hideShowTypeOther()
{
	if($('#editSampleType').val() == 0 || $('#editSampleType').val() == 7)
	{
		$('#editSampleTypeOther').show();
	}
	else
	{
		$('#editSampleTypeOther').hide();
	}
}
</script>