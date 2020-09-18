#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use JSON::XS;
use DBI;
use lib "lib/";
use lib "lib/pangu";
use pangu;
use config;
use userCookie;

my $userCookie = new userCookie;
my $userId = (cookie('cid')) ? $userCookie->checkCookie(cookie('cid')) : 0;
exit if (!$userId);

my $commoncfg = readConfig("main.conf");
my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});
my $config = new config;
my $RSDashboardURL = $config->getFieldValueWithFieldName("RSDashboardURL");

my $smrtrunId = param ('smrtrunId') || '';
my %status = (
	0=>'Initialized',
	1=>'Started',
	2=>'Completed'
	);
my %movieLength = (
	30=>'30 mins',
	45=>'45 mins',
	60=>'60 mins',
	90=>'90 mins',
	120=>'120 mins',
	150=>'150 mins',
	180=>'180 mins',
	240=>'240 mins'
	);

print header;
if ($smrtrunId)
{
	my $smrtrun=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$smrtrun->execute($smrtrunId);
	my @smrtrun = $smrtrun->fetchrow_array();
	my $smrtrunDetails = decode_json $smrtrun[8];
	$smrtrunDetails->{'comments'} = '' unless (exists $smrtrunDetails->{'comments'});
	$smrtrunDetails->{'comments'} = escapeHTML($smrtrunDetails->{'comments'});
	$smrtrunDetails->{'comments'} =~ s/\n/<br>/g;
	my $relatedSmrtwells = '';
	my $smrtwellInRun=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'smrtwell' AND z = ? ORDER BY name");
	$smrtwellInRun->execute($smrtrunId);
	while (my @smrtwellInRun = $smrtwellInRun->fetchrow_array())
	{
		my $smrtwellDetails = decode_json $smrtwellInRun[8];
		$smrtwellDetails->{'loadingName'} = '' unless (exists $smrtwellDetails->{'loadingName'});
		$smrtwellDetails->{'customizedMovieLength'} = '' unless (exists $smrtwellDetails->{'customizedMovieLength'});
		$relatedSmrtwells .= ($smrtwellInRun[5] == 0) ? 
		"<li><a onclick='closeDialog();openDialog(\"smrtwellView.cgi?smrtwellId=$smrtwellInRun[0]\")' title='View'>$smrtwellInRun[2] ($smrtwellDetails->{'loadingName'} x $smrtwellInRun[4], movie length: $smrtwellDetails->{'customizedMovieLength'})</a></li>"
		: "<li><a onclick='closeDialog();openDialog(\"smrtwellView.cgi?smrtwellId=$smrtwellInRun[0]\")' title='View'>$smrtwellInRun[2] ($smrtwellDetails->{'loadingName'} x $smrtwellInRun[4], movie length: $movieLength{$smrtwellInRun[5]})</a></li>";	
	}
	my $RSDashboard = ($smrtrun[7])? "<a style='float: right;' href='$RSDashboardURL/Metrics/RSRunReport?instrument=&run=$smrtrun[2]&from=20000101-000000&to=20000101-000000' title='Go to RSDashboard' target='_blank'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-extlink'></span>Go to RSDashboard</a>" : "";

	undef $/;# enable slurp mode
	my $html = <DATA>;
	$html =~ s/\$smrtrunId/$smrtrunId/g;
	$html =~ s/\$smrtrunName/$smrtrun[2]/g;
	$html =~ s/\$smrtrunSmrtcell/$smrtrun[3]/g;
	$html =~ s/\$relatedSmrtwells/$relatedSmrtwells/g;
	$html =~ s/\$comments/$smrtrunDetails->{'comments'}/g;
	$html =~ s/\$status/$status{$smrtrun[7]}/g;
	$html =~ s/\$RSDashboard/$RSDashboard/g;
	$html =~ s/\$smrtrunCreator/$smrtrun[9]/g;
	$html =~ s/\$smrtrunEnteredDate/$smrtrun[10]/g;
	$html =~ s/\$\$/$$/g;

	print $html;
}
else
{
	print '402 Invalid operation';
	exit;
}

__DATA__
<table>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="viewSmrtrunName"><b>SMRT Run</b></label></td>
		<td>$smrtrunName $RSDashboard</td>
	</tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="viewSmrtrunSmrtcell"><b>SMRTCells</b></label></td>
		<td>$smrtrunSmrtcell</td>
	</tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="viewSmrtrunWells"><b>Sample Plate Wells</b></label></td>
		<td><ul style='margin: 0;'>$relatedSmrtwells</ul></td>
	</tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="viewSmrtrunStatus"><b>Status</b></label></td>
		<td>$status</td>
	</tr>
	<tr><td style='text-align:right;white-space: nowrap;'><label for="viewSmrtrunComments"><b>Comments</b></label></td>
		<td>$comments</td>
	</tr>
	<tr><td style='text-align:right;white-space: nowrap;'></td>
		<td><sup class='ui-state-disabled'>Last changed by $smrtrunCreator on $smrtrunEnteredDate</sup></td>
	</tr>
</table>
<script>
$('#dialog').dialog("option", "title", "View SMRT Run $smrtrunName");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Edit", click: function() { closeDialog();openDialog("smrtrunEdit.cgi?smrtrunId=$smrtrunId"); } }, { text: "OK", click: function() {closeDialog(); } } ] );
</script>
