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

my %purpose = (
	0=>'General Item',
	1=>'Sample for PacBio Library'
	);
my $samplePurpose = '';
foreach (sort {$a <=> $b} keys %purpose)
{
	$samplePurpose .= "<option value='$_'>$purpose{$_}</option>";
}
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
my $sampleTypeList = '';
foreach (sort {$a <=> $b} keys %type)
{
	next unless ($_);
	$sampleTypeList .= "<option value='$_'>$_. $type{$_}</option>";
}
$sampleTypeList .= "<option value='0'>$type{0}</option>";
my %status = (
	0=>'Na',
	1=>'Status One',
	2=>'Status Two'
	);
my $sampleStatus = '';
foreach (sort {$a <=> $b} keys %status)
{
	$sampleStatus .= "<option value='$_'>$status{$_}</option>";
}
my $noun = 'word';

my $serviceId = param ('serviceId') || '';
my $service=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$service->execute($serviceId);
my @service = $service->fetchrow_array();
my $serviceToProject=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$serviceToProject->execute($service[6]);
my @serviceToProject = $serviceToProject->fetchrow_array();

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$serviceId/$serviceId/g;
$html =~ s/\$serviceName/$service[2]/g;
$html =~ s/\$projectName/$serviceToProject[2]/g;
$html =~ s/\$samplePurpose/$samplePurpose/g;
$html =~ s/\$sampleType/$sampleTypeList/g;
$html =~ s/\$sampleStatus/$sampleStatus/g;
$html =~ s/\$noun/$noun/g;

print header;
print $html;

__DATA__
<form id="newSample" name="newSample" action="sampleSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	
	<table>
	<tr><td style='text-align:right'><label for="newSampleServiceId">For</label></td><td><a title='Project'>$projectName</a> > <a title='Service'>$serviceName</a><input name="serviceId" id="newSampleServiceId" type="hidden" value="$serviceId" /></tr>
	<tr><td style='text-align:right'><label for="newSampleName"><b>Sample Name</b></label></td><td><input class='ui-widget-content ui-corner-all' name="name" id="newSampleName" size="10" type="text" maxlength="32" value="" /></td></tr>
	<tr><td style='text-align:right'><label for="newSamplePurpose"><b>Purpose</b></label></td><td><select class='ui-widget-content ui-corner-all' name="purpose" id="newSamplePurpose">$samplePurpose</select></td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="newSampleType"><b>Type</b></label></td>
		<td>
			<select class='ui-widget-content ui-corner-all' name="type" id="newSampleType" onchange="hideShowTypeOther();">$sampleType</select>
			<br><textarea class='ui-widget-content ui-corner-all' name="sampleTypeOther" id="newSampleTypeOther" cols="40" rows="2" placeholder=""></textarea>
		</td>
	</tr>
	<tr><td style='text-align:right'><label for="newSampleStatus"><b>Status</b></label></td><td><select class='ui-widget-content ui-corner-all' name="status" id="newSampleStatus">$sampleStatus</select></td></tr>
	<tr><td style='text-align:right'><label for="newSampleDescription"><b>Description</b></label><br><sub id='newSampleDescription_count' style='display:none'></sub></td><td><textarea class='ui-widget-content ui-corner-all word_count' name="description" id="newSampleDescription" cols="50" rows="10"></textarea></td></tr>
	</table>
</form>
<script>
hideShowTypeOther();
$('#dialog').dialog("option", "title", "New Sample");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('newSample'); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
wordCount('$noun');
function hideShowTypeOther()
{
	if($('#newSampleType').val() == 0 || $('#newSampleType').val() == 7)
	{
		$('#newSampleTypeOther').show();
	}
	else
	{
		$('#newSampleTypeOther').hide();
	}
}

</script>