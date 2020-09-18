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
if(!$userId)
{
	print header;
	print <<eof;
	<h2>Error</h2>
	Please <a onclick="openDialog('login.cgi');">login</a> first!
	<script type="text/javascript">
		openDialog('login.cgi');
		loadingHide();
	</script>
eof
    exit;
}

my $commoncfg = readConfig("main.conf");
my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});
my $config = new config;
my $RSDashboardURL = $config->getFieldValueWithFieldName("RSDashboardURL");

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
my %chemistry = (
	P4=>'P4',
	P5=>'P5',
	P6=>'P6'
	);
my %condition = (
	0=>'Customized:',
	1=>'100%',
	2=>'75%'
	);

my $smrtruns = '';
my $smrtrun=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'smrtrun' ORDER BY name");
$smrtrun->execute();
while (my @smrtrun = $smrtrun->fetchrow_array())
{
	my $smrtrunDetails = decode_json $smrtrun[8];
	$smrtrunDetails->{'comments'} = '' unless (exists $smrtrunDetails->{'comments'});
	$smrtrunDetails->{'comments'} = escapeHTML($smrtrunDetails->{'comments'});
	$smrtrunDetails->{'comments'} =~ s/\n/<br>/g;
	$smrtruns = "<form id='smrtrunList$$' name='smrtrunList$$'>
		<table id='smrtruns$$' class='display' style='width: 100%;'>
			<thead>
				<tr>
					<th style='text-align:left'><b>SMRT Run</b></th>
					<th style='text-align:left'><b>Sample Plate Wells</b></th>
					<th style='text-align:left'><b>Status</b></th>
					<th style='text-align:left'><b>Creator</b></th>
				</tr>
			</thead>
			<tbody>" unless($smrtruns);
	my $relatedSmrtwells = '';
	my $smrtwellInRun=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'smrtwell' AND z = ? ORDER BY name");
	$smrtwellInRun->execute($smrtrun[0]);
	while (my @smrtwellInRun = $smrtwellInRun->fetchrow_array())
	{
		my $smrtwellDetails = decode_json $smrtwellInRun[8];
		$smrtwellDetails->{'loadingName'} = '' unless (exists $smrtwellDetails->{'loadingName'});
		$smrtwellDetails->{'customizedMovieLength'} = '' unless (exists $smrtwellDetails->{'customizedMovieLength'});
		$movieLength{$smrtwellInRun[5]} = $smrtwellDetails->{'customizedMovieLength'} unless ($smrtwellInRun[5]);
		$smrtwellDetails->{'concentration'} = '' unless (exists $smrtwellDetails->{'concentration'});
		$smrtwellDetails->{'polRatio'} = '' unless (exists $smrtwellDetails->{'polRatio'});
		$smrtwellDetails->{'chemistry'} = '' unless (exists $smrtwellDetails->{'chemistry'});
		$smrtwellDetails->{'customizedCondition'} = '' unless (exists $smrtwellDetails->{'customizedCondition'});
		$condition{$smrtwellInRun[7]} = $smrtwellDetails->{'customizedCondition'} unless ($smrtwellInRun[7]);
		my $paclib = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$paclib->execute($smrtwellInRun[3]);
		my @paclib=$paclib->fetchrow_array();
		my $sample;
		my $service;
		my $project;	
		my $paclibToSample=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$paclibToSample->execute($paclib[6]);
		my @paclibToSample = $paclibToSample->fetchrow_array();
		$sample = $paclibToSample[2] if ($paclibToSample[2]);
		if($paclibToSample[0])
		{
			my $sampleToService=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
			$sampleToService->execute($paclibToSample[6]);
			my @sampleToService = $sampleToService->fetchrow_array();
			$service = $sampleToService[2] if ($sampleToService[2]);
			if($sampleToService[0])
			{
				my $serviceToProject=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
				$serviceToProject->execute($sampleToService[6]);
				my @serviceToProject = $serviceToProject->fetchrow_array();
				$project = $serviceToProject[2] if ($serviceToProject[2]);
			}
		}
		$relatedSmrtwells = "<tr><th>Well</th><th>Sample name</th><th>Cell number</th><th>Movie length</th><th>DNA concentration</th><th>Polymerase ratio</th><th>Chemistry</th><th>Condition</th></tr>" unless ($relatedSmrtwells);
		$relatedSmrtwells .= "<tr>
			<td><a onclick='openDialog(\"smrtwellView.cgi?smrtwellId=$smrtwellInRun[0]\")' title='View'>$smrtwellInRun[2]</a></td>
			<td title='SMRT sample name'>$smrtwellDetails->{'loadingName'}</td>
			<td title='Cell number'>$smrtwellInRun[4]</td>
			<td title='Movie legnth'>$movieLength{$smrtwellInRun[5]}</td>
			<td title='DNA concentration'>$smrtwellDetails->{'concentration'} nM</td>
			<td title='Polymerase ratio'>$smrtwellDetails->{'polRatio'}</td>
			<td title='Chemistry'>$chemistry{$smrtwellDetails->{'chemistry'}}</td>
			<td title='Condition'>$condition{$smrtwellInRun[7]}</td>
			</tr>";	
	}
	my $RSDashboard = ($smrtrun[7])? "<a style='float: right;' href='$RSDashboardURL/Metrics/RSRunReport?instrument=&run=$smrtrun[2]&from=20000101-000000&to=20000101-000000' title='Go to RSDashboard' target='_blank'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-extlink'></span></a>" : "";
	$smrtruns .= "<tr>
		<td title='Run Name'><div style='position: relative;'><a id='smrtrunId$smrtrun[0]$$' onclick='openDialog(\"smrtrunView.cgi?smrtrunId=$smrtrun[0]\")' title='View'>$smrtrun[2]</a></div></td>
		<td title='Wells'>
			<table style='margin: 0;'>
				$relatedSmrtwells
			</table>
		</td>
		<td>$status{$smrtrun[7]} $RSDashboard</td>
		<td title='Creation Date: $smrtrun[10]'>$smrtrun[9]</td>
		</tr>";
}
$smrtruns .= "</tbody></table></form>\n" if($smrtruns);
my $button = "<div class='ui-state-highlight ui-corner-all' style='padding: 0 .7em;'><button style='float: right; margin-top: .3em; margin-right: .3em;' onclick='refresh(\"general\")'>Refresh</button>
			<h2>SMRT Runs</h2>";

unless($smrtruns)
{
	$button .= "<p class='ui-state-error ui-corner-all' style='padding: .7em;'><span class='ui-icon ui-icon-alert' style='float: left; margin-right: .3em;'></span>
			<strong>No runs, please create one!</strong></p>";
}
$button .= "</div>\n";

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$button/$button/g;
$html =~ s/\$smrtruns/$smrtruns/g;
$html =~ s/\$\$/$$/g;

print header(-cookie=>cookie(-name=>'general',-value=>2));
print $html;

__DATA__
$button
$smrtruns
<script>
buttonInit();
$( "#smrtruns$$" ).dataTable({
	"scrollY": "600px",
	"scrollCollapse": true,
	"paging": false
});
loadingHide();
</script>