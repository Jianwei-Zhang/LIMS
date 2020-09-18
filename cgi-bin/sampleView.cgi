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
my $sampleId = param ('sampleId') || '';
my $sample = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$sample->execute($sampleId);
my @sample=$sample->fetchrow_array();
my $sampleDetails = decode_json $sample[8];
$sampleDetails->{'sampleTypeOther'} = '' unless (exists $sampleDetails->{'sampleTypeOther'});
$sampleDetails->{'description'} = '' unless (exists $sampleDetails->{'description'});
$sampleDetails->{'description'} = escapeHTML($sampleDetails->{'description'});
$sampleDetails->{'description'} =~ s/\n/<br>/g;

my $sampleToService=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$sampleToService->execute($sample[6]);
my @sampleToService = $sampleToService->fetchrow_array();
my $serviceToProject=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$serviceToProject->execute($sampleToService[6]);
my @serviceToProject = $serviceToProject->fetchrow_array();

my $sampleLocation = 'Unknown.';
my $sampleToBox=$dbh->prepare("SELECT matrix.* FROM link,matrix WHERE link.type LIKE 'box' AND link.parent = matrix.id AND link.child = ?");
$sampleToBox->execute($sample[0]);
my @sampleToBox = $sampleToBox->fetchrow_array();
if($sampleToBox[0])
{
	my $boxToFreezer=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$boxToFreezer->execute($sampleToBox[3]);
	my @boxToFreezer = $boxToFreezer->fetchrow_array();
	my $freezerToRoom=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$freezerToRoom->execute($boxToFreezer[3]);
	my @freezerToRoom = $freezerToRoom->fetchrow_array();
	$sampleLocation = "$freezerToRoom[2] > Freezer $boxToFreezer[2] > Box $sampleToBox[2]";
}
my $relatedLibs = '';
if($sample[3] == 1) #for paclib
{
	my %paclibStatus = (
		0=>'Na',
		1=>'Status One',
		2=>'Status Two'
		);
	my $paclibForSample=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'paclib' AND z = ? ORDER BY o");
	$paclibForSample->execute($sample[0]);
	while (my @paclibForSample = $paclibForSample->fetchrow_array())
	{
		$relatedLibs .= "<li><a onclick='closeDialog();openDialog(\"paclibView.cgi?paclibId=$paclibForSample[0]\")' title='View'>$paclibForSample[2]</a> (Status: $paclibStatus{$paclibForSample[4]})</li>";	
	}
	$relatedLibs = "<tr><td style='text-align:right'><label for='viewSamplePaclibs'><b>Related Paclibs</b></label></td><td><ul style='margin: 0;'>$relatedLibs</ul></td></tr>";
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$sampleId/$sampleId/g;
$html =~ s/\$project/$serviceToProject[2]/g;
$html =~ s/\$serviceId/$sample[6]/g;
$html =~ s/\$serviceName/$sampleToService[2]/g;
$html =~ s/\$sampleName/$sample[2]/g;
$html =~ s/\$samplePurpose/$purpose{$sample[3]}/g;
$html =~ s/\$sampleTypeOther/$sampleDetails->{'sampleTypeOther'}/g;
$html =~ s/\$sampleType/$type{$sample[4]}/g;
$html =~ s/\$sampleStatus/$status{$sample[5]}/g;
$html =~ s/\$sampleLocation/$sampleLocation/g;
$html =~ s/\$relatedLibs/$relatedLibs/g;
$html =~ s/\$sampleBarcode/$sample[7]/g;
$html =~ s/\$sampleDescription/$sampleDetails->{'description'}/g;
$html =~ s/\$sampleCreator/$sample[9]/g;
$html =~ s/\$sampleCreationDate/$sample[10]/g;

print header;
print $html;

__DATA__
<table>
	<tr><td style='text-align:right;white-space: nowrap;'><b>For</b></td><td><a onclick='closeDialog();openDialog("serviceView.cgi?serviceId=$serviceId")' title='View Service'>$project > $serviceName</a></td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="viewSampleName"><b>Sample Name</b></label></td><td>$sampleName<br><img alt='$sampleBarcode' src='barcode.cgi?code=$sampleBarcode'/><sup class='ui-state-disabled'>by $sampleCreator on $sampleCreationDate</sup></td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="viewSamplePurpose"><b>Purpose</b></label></td><td>$samplePurpose</td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="viewSampleType"><b>Type</b></label></td><td>$sampleType$sampleTypeOther</td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="viewSampleStatus"><b>Status</b></label></td><td>$sampleStatus</td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="viewSampleLocation"><b>Container</b></label></td><td>$sampleLocation</td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="viewSampleDescription"><b>Description</b></label><br></td><td>$sampleDescription</td></tr>
	$relatedLibs
</table>
<script>
$('#dialog').dialog("option", "title", "View Sample");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Edit", click: function() { closeDialog();openDialog("sampleEdit.cgi?sampleId=$sampleId"); } }, { text: "OK", click: function() {closeDialog(); } } ] );
</script>