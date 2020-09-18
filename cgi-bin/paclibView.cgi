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
$paclibDetails->{'shearingOperator'} = 'Unknown.' unless (exists $paclibDetails->{'shearingOperator'});
$paclibDetails->{'bluepippinInput'} = '' unless (exists $paclibDetails->{'bluepippinInput'});
$paclibDetails->{'bluepippinSize'} = '' unless (exists $paclibDetails->{'bluepippinSize'});
$paclibDetails->{'bluepippinOutput'} = '' unless (exists $paclibDetails->{'bluepippinOutput'});
$paclibDetails->{'bluepippinConcentration'} = '' unless (exists $paclibDetails->{'bluepippinConcentration'});
$paclibDetails->{'bluepippinOperator'} = 'Unknown.' unless (exists $paclibDetails->{'bluepippinOperator'});
$paclibDetails->{'description'} = '' unless (exists $paclibDetails->{'description'});

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

my $paclibLocation = 'Unknown.';
my $paclibToBox=$dbh->prepare("SELECT matrix.* FROM link,matrix WHERE link.type LIKE 'box' AND link.parent = matrix.id AND link.child = ?");
$paclibToBox->execute($paclib[0]);
my @paclibToBox = $paclibToBox->fetchrow_array();
if($paclibToBox[0])
{
	my $boxToFreezer=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$boxToFreezer->execute($paclibToBox[3]);
	my @boxToFreezer = $boxToFreezer->fetchrow_array();
	my $freezerToRoom=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$freezerToRoom->execute($boxToFreezer[3]);
	my @freezerToRoom = $freezerToRoom->fetchrow_array();
	$paclibLocation = "$freezerToRoom[2] > Freezer $boxToFreezer[2] > Box $paclibToBox[2]";
}

my %smrtrunStatus = (
	0=>'Initialized',
	1=>'Started',
	2=>'Completed'
	);
my $relatedSmrtwells = '';
my $smrtwellForPaclib=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'smrtwell' AND o = ? ORDER BY name");
$smrtwellForPaclib->execute($paclib[0]);
while (my @smrtwellForPaclib = $smrtwellForPaclib->fetchrow_array())
{
	my $smrtrun=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$smrtrun->execute($smrtwellForPaclib[6]);
	my @smrtrun = $smrtrun->fetchrow_array();
	$relatedSmrtwells .= "<li><a onclick='closeDialog();openDialog(\"smrtwellView.cgi?smrtwellId=$smrtwellForPaclib[0]\")' title='View'>SMRT Run $smrtrun[2] > $smrtwellForPaclib[2]</a> (Status:$smrtrunStatus{$smrtrun[7]})</li>";	
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$paclibId/$paclibId/g;
$html =~ s/\$sampleId/$paclib[6]/g;
$html =~ s/\$project/$serviceToProject[2]/g;
$html =~ s/\$serviceName/$sampleToService[2]/g;
$html =~ s/\$sampleName/$paclibToSample[2]/g;
$html =~ s/\$paclibName/$paclib[2]/g;
$html =~ s/\$paclibStatus/$status{$paclib[4]}/g;
$html =~ s/\$paclibLocation/$paclibLocation/g;
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
$html =~ s/\$relatedSmrtwells/$relatedSmrtwells/g;
$html =~ s/\$paclibCreator/$paclib[9]/g;
$html =~ s/\$paclibCreationDate/$paclib[10]/g;
$html =~ s/\$noun/$noun/g;

print header;
print $html;

__DATA__
<table>
	<tr><td style='text-align:right;white-space: nowrap;'><b>For</b></td><td><a onclick='closeDialog();openDialog("sampleView.cgi?sampleId=$sampleId")' title='View Sample'>$project > $serviceName > $sampleName</a></td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="viewPaclibName"><b>Paclib Name</b></label></td><td>$paclibName<br><img alt='$paclibBarcode' src='barcode.cgi?code=$paclibBarcode'/><sup class='ui-state-disabled'>by $paclibCreator on $paclibCreationDate</sup></td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="viewPaclibDate"><b>Library Date</b></label></td><td>$paclibDate</td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><b>Shearing</b></td><td>
		DNA input: $shearingInput ug<br>
		RPM: $shearingRpm<br>
		DNA output: $shearingOutput ug<br>
		PB beads steps: $shearingBeadsSteps<br>
		Operator: $shearingOperator
		</td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><b>BluePippin</b></td><td>
		DNA input: $bluepippinInput ug<br>
		Size selection: $bluepippinSize kb<br>
		DNA output: $bluepippinOutput ug<br>
		DNA concentration: $bluepippinConcentration ng/uL<br>
		Operator: $bluepippinOperator
		</td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="viewPaclibStatus"><b>Status</b></label></td><td>$paclibStatus</td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="viewPaclibLocation"><b>Container</b></label></td><td>$paclibLocation</td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="viewPaclibDescription"><b>Description</b></label><br></td><td>$paclibDescription</td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="viewPaclibSmrtrun"><b>Related Wells</b></label></td><td><ul style='margin: 0;'>$relatedSmrtwells</ul></td></tr>
</table>
<script>
$('#dialog').dialog("option", "title", "View Paclib");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Edit", click: function() { closeDialog();openDialog("paclibEdit.cgi?paclibId=$paclibId"); } }, { text: "OK", click: function() {closeDialog(); } } ] );
wordCount('$noun');
</script>