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

my %yesOrNo = ( 0=>'<span style="float: right;" class="ui-icon ui-icon-cancel" title="Not yet"></span>', 1=>'<span style="float: right;" class="ui-icon ui-icon-check" title="Paid"></span>' );
my %status = ( 1=>'Active', 0=>'Inactive' );
my %serviceType = (
	0=>'Other',
	1=>'DE: DNA extraction',
	2=>'BL: BAC library construction',
	3=>'BS: BAC (pool) sequencing',
	4=>'BE: BAC end sequencing',
	5=>'WS: WGS sequencing',
	6=>'RE: RNA extraction',
	7=>'IS: Iso-Seq',
	8=>'MP: Illumina mate pair library',
	9=>'PE: Illumina paired end library',
	10=>'NS: Illumina sequencing',
	11=>'WP: Whole Genome Profiling',
	12=>'PA: PCR amplicon sequencing'
	);
my %paymentType = (
	0=>'NA',
	1=>'Check',
	2=>'Credit card, 3% fee added',
	3=>'Purchase order',
	4=>'Wire transfer, fee of $50 added'
	);

my $services = '';
my $service=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'service'");# ORDER BY name
$service->execute();
if ($service->rows > 0)
{
	my $serviceResult;
	while (my @service = $service->fetchrow_array())
	{
		@{$serviceResult->{$service[2]}} = @service;
	}
	foreach (sort {uc ($a) cmp uc($b)} keys %$serviceResult)
	{
		my @service = @{$serviceResult->{$_}};
		my $serviceDetails = decode_json $service[8];
		$serviceDetails->{'serviceTypeOther'} = '' unless (exists $serviceDetails->{'serviceTypeOther'});
		$serviceDetails->{'amount'} = '' unless (exists $serviceDetails->{'amount'});
		$serviceDetails->{'assembly'} = '' unless (exists $serviceDetails->{'assembly'});
		$serviceDetails->{'analysis'} = '' unless (exists $serviceDetails->{'analysis'});
		$serviceDetails->{'serviceNote'} = '' unless (exists $serviceDetails->{'serviceNote'});
		$serviceDetails->{'genus'} = '' unless (exists $serviceDetails->{'genus'});
		$serviceDetails->{'species'} = '' unless (exists $serviceDetails->{'species'});
		$serviceDetails->{'subspecies'} = '' unless (exists $serviceDetails->{'subspecies'});
		$serviceDetails->{'commonName'} = '' unless (exists $serviceDetails->{'commonName'});
		$serviceDetails->{'genomeSize'} = '' unless (exists $serviceDetails->{'genomeSize'});
		$serviceDetails->{'sampleTypeOther'} = '' unless (exists $serviceDetails->{'sampleTypeOther'});
		$serviceDetails->{'sampleName'} = '' unless (exists $serviceDetails->{'sampleName'});
		$serviceDetails->{'sampleNote'} = '' unless (exists $serviceDetails->{'sampleNote'});
		$serviceDetails->{'paymentType'} = '0' unless (exists $serviceDetails->{'paymentType'});
		$serviceDetails->{'poNumber'} = '' unless (exists $serviceDetails->{'poNumber'});
		$serviceDetails->{'paymentDate'} = '' unless (exists $serviceDetails->{'paymentDate'});
		$serviceDetails->{'piName'} = '' unless (exists $serviceDetails->{'piName'});
		$serviceDetails->{'piEmail'} = '' unless (exists $serviceDetails->{'piEmail'});
		$serviceDetails->{'institution'} = '' unless (exists $serviceDetails->{'institution'});
		$serviceDetails->{'contactName'} = '' unless (exists $serviceDetails->{'contactName'});
		$serviceDetails->{'contactPhone'} = '' unless (exists $serviceDetails->{'contactPhone'});
		$serviceDetails->{'contactEmail'} = '' unless (exists $serviceDetails->{'contactEmail'});
		$serviceDetails->{'contactAddress'} = '' unless (exists $serviceDetails->{'contactAddress'});
		$serviceDetails->{'submitDate'} = '' unless (exists $serviceDetails->{'submitDate'});
		$serviceDetails->{'submitPerson'} = '' unless (exists $serviceDetails->{'submitPerson'});
		$serviceDetails->{'startDate'} = '' unless (exists $serviceDetails->{'startDate'});
		$serviceDetails->{'endDate'} = '' unless (exists $serviceDetails->{'endDate'});
		$serviceDetails->{'comments'} = '' unless (exists $serviceDetails->{'comments'});
		$serviceDetails->{'serviceNote'} = escapeHTML($serviceDetails->{'serviceNote'});
		$serviceDetails->{'serviceNote'} =~ s/\n/<br>/g;
		$serviceDetails->{'sampleNote'} = escapeHTML($serviceDetails->{'sampleNote'});
		$serviceDetails->{'sampleNote'} =~ s/\n/<br>/g;
		$serviceDetails->{'contactAddress'} = escapeHTML($serviceDetails->{'contactAddress'});
		$serviceDetails->{'contactAddress'} =~ s/\n/<br>/g;
		$serviceDetails->{'comments'} = escapeHTML($serviceDetails->{'comments'});
		$serviceDetails->{'comments'} =~ s/\n/<br>/g;
		$serviceDetails->{'serviceTypeOther'} = '' if ($service[3] > 0);
		my $project=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$project->execute($service[6]);
		my @project = $project->fetchrow_array();

		my $newSample = "<a style='float: right;' class='ui-state-highlight ui-corner-all' onclick='openDialog(\"sampleNew.cgi?serviceId=$service[0]\")' title='Add New Sample'>
				<span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-plus' title='Add New Sample'></span></a>";
		my $relatedSamples = '';
		my $sampleInService=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sample' AND z = ? ORDER BY o");
		$sampleInService->execute($service[0]);
		while (my @sampleInService = $sampleInService->fetchrow_array())
		{
			my %sampleType = (
				0=>'Other:',
				1=>'Genomic DNA',
				2=>'Total RNA',
				3=>'BAC clones, as bacteria',
				4=>'BAC plasmids',
				5=>'PCR amplicons',
				6=>'cDNA',
				7=>'Tissue(s):'
				);
			my $sampleDetails = decode_json $sampleInService[8];
			$sampleDetails->{'sampleTypeOther'} = '' unless (exists $sampleDetails->{'sampleTypeOther'});
			$sampleDetails->{'description'} = '' unless (exists $sampleDetails->{'description'});
			$sampleDetails->{'description'} = escapeHTML($sampleDetails->{'description'});
			$sampleDetails->{'description'} =~ s/\n/<br>/g;
			$sampleInService[4] = ($sampleInService[4] == 0 || $sampleInService[4] == 7) ? "$sampleType{$sampleInService[4]}$sampleDetails->{'sampleTypeOther'}" : "$sampleType{$sampleInService[4]}";
			$relatedSamples .= "<li><a onclick='openDialog(\"sampleView.cgi?sampleId=$sampleInService[0]\")' title='View'>$sampleInService[2]-$sampleInService[4]</a></li>";	
		}
		$services = "
			<table id='services$$' class='display' style='width: 100%;'>
				<thead>
					<tr>
						<th style='text-align:left'><b>Project</b></th>
						<th style='text-align:left'><b>Service</b></th>
						<th style='text-align:left'><b>Type</b></th>
						<th style='text-align:left'><b>Related samples</b></th>
						<th style='text-align:left'><b>P.I.</b></th>
						<th style='text-align:left'><b>Status</b></th>
						<th style='text-align:left'><b>Paid</b></th>
						<th style='text-align:left'><b>Creator</b></th>
					</tr>
				</thead>
				<tbody>" unless($services);
		my $invoiceList = "";
		for (split /\s*,\s*/ ,$serviceDetails->{'invoice'})
		{
			$invoiceList .= ($invoiceList) ? ", <a title='Invoice#$_'>$_</a>" : "<a title='Invoice#$_'>$_</a>";
		}
		$services .= "<tr>
			<td title='Project id: $project[0]'>$project[2]</td>
			<td title='Service'><div style='position: relative;'><a id='serviceId$service[0]$$' onclick='openDialog(\"serviceView.cgi?serviceId=$service[0]\")' title='View'>$service[2]</a></div></td>
			<td title='$serviceDetails->{'serviceTypeOther'}'>$serviceType{$service[3]}</td>
			<td title='Related samples'>$newSample
				<ul style='margin: 0;'>$relatedSamples</ul>
			</td>
			<td title='Send an Email'><a href='mailto:$serviceDetails->{'piEmail'}'>$serviceDetails->{'piName'}</a></td>
			<td title='Duration: $serviceDetails->{'startDate'} - $serviceDetails->{'endDate'}'>$status{$service[7]}</td>
			<td title='via $paymentType{$serviceDetails->{'paymentType'}}'>$invoiceList$yesOrNo{$service[5]}</td>
			<td title='Creation Date: $service[10]'>$service[9]</td>
			</tr>";
	}
}
$services .= "</tbody></table>\n" if($services);
my $button = "<div class='ui-state-highlight ui-corner-all' style='padding: 0 .7em;'>
	<button style='float: right; margin-top: .3em; margin-right: .3em;' onclick='refresh(\"general\")'>Refresh</button>
			<h2>Services</h2>";

unless($services)
{
	$button .= "<p class='ui-state-error ui-corner-all' style='padding: .7em;'><span class='ui-icon ui-icon-alert' style='float: left; margin-right: .3em;'></span>
			<strong>No service, please create one!</strong></p>";
}
$button .= "</div>\n";

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$button/$button/g;
$html =~ s/\$services/$services/g;
$html =~ s/\$\$/$$/g;


print header(-cookie=>cookie(-name=>'general',-value=>1));
print $html;

__DATA__
$button
$services
<script>
buttonInit();
$( "#services$$" ).dataTable({
	"scrollY": "600px",
	"scrollCollapse": true,
	"paging": false
});
loadingHide();
</script>