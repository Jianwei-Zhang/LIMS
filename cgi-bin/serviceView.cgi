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

my $serviceId = param ('serviceId') || '';
my %yesOrNo = ( 0=>'No', 1=>'Yes' );
my %status = ( 1=>'Active', 0=>'Inactive');
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
	2=>'Crview card, 3% fee added',
	3=>'Purchase order',
	4=>'Wire transfer, fee of $50 added'
	);
print header;
if ($serviceId)
{
	my $service=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$service->execute($serviceId);
	my @service = $service->fetchrow_array();
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
	$serviceDetails->{'paymentType'} = '0' unless (exists $serviceDetails->{'paymentType'});
	$serviceDetails->{'poNumber'} = '' unless (exists $serviceDetails->{'poNumber'});
	$serviceDetails->{'paymentDate'} = '' unless (exists $serviceDetails->{'paymentDate'});
	$serviceDetails->{'invoice'} = '' unless (exists $serviceDetails->{'invoice'});
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
	$serviceDetails->{'contactAddress'} = escapeHTML($serviceDetails->{'contactAddress'});
	$serviceDetails->{'contactAddress'} =~ s/\n/<br>/g;
	$serviceDetails->{'comments'} = escapeHTML($serviceDetails->{'comments'});
	$serviceDetails->{'comments'} =~ s/\n/<br>/g;

	my $project=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$project->execute($service[6]);
	my @project = $project->fetchrow_array();

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
		$relatedSamples .= "<li><a onclick='closeDialog();openDialog(\"sampleView.cgi?sampleId=$sampleInService[0]\")' title='View'>$sampleInService[2]-$sampleInService[4]</a></li>";	
	}

	undef $/;# enable slurp mode
	my $html = <DATA>;
	$html =~ s/\$serviceId/$serviceId/g;
	$html =~ s/\$serviceName/$service[2]/g;
	if($service[3] == 0)
	{
		$html =~ s/\$serviceTypeOther/$serviceDetails->{'serviceTypeOther'}/g;
	}
	else
	{
		$html =~ s/\$serviceTypeOther//g;
	}
	$html =~ s/\$serviceType/$serviceType{$service[3]}/g;
	$html =~ s/\$amount/$serviceDetails->{'amount'}/g;
	$html =~ s/\$assembly/$serviceDetails->{'assembly'}/g;
	$html =~ s/\$analysis/$serviceDetails->{'analysis'}/g;
	$html =~ s/\$serviceNote/$serviceDetails->{'serviceNote'}/g;
	$html =~ s/\$genus/$serviceDetails->{'genus'}/g;
	$html =~ s/\$species/$serviceDetails->{'species'}/g;
	$html =~ s/\$subspecies/$serviceDetails->{'subspecies'}/g;
	$html =~ s/\$commonName/$serviceDetails->{'commonName'}/g;
	$html =~ s/\$genomeSize/$serviceDetails->{'genomeSize'}/g;
	$html =~ s/\$relatedSamples/$relatedSamples/g;
	$html =~ s/\$paymentType/$paymentType{$serviceDetails->{'paymentType'}}/g;
	if($serviceDetails->{'paymentType'} == 3)
	{
		$html =~ s/\$poNumber/:$serviceDetails->{'poNumber'}/g;
	}
	else
	{
		$html =~ s/\$poNumber//g;
	}
	if($service[5])
	{
		$html =~ s/\$paymentStatus/Fully paid on /g;
		$html =~ s/\$paymentDate/$serviceDetails->{'paymentDate'}/g;
	}
	else
	{
		$html =~ s/\$paymentStatus/Not paid./g;
		$html =~ s/\$paymentDate//g;
	}
	my $invoiceList = "";
	for (split /\s*,\s*/ ,$serviceDetails->{'invoice'})
	{
		$invoiceList .= ($invoiceList) ? ", <a title='Invoice#$_'>$_</a>" : "<a title='Invoice#$_'>$_</a>";
	}
	$html =~ s/\$invoice/$invoiceList/g;
	$html =~ s/\$piName/$serviceDetails->{'piName'}/g;
	$html =~ s/\$piEmail/$serviceDetails->{'piEmail'}/g;
	$html =~ s/\$institution/$serviceDetails->{'institution'}/g;
	$html =~ s/\$contactName/$serviceDetails->{'contactName'}/g;
	$html =~ s/\$contactPhone/$serviceDetails->{'contactPhone'}/g;
	$html =~ s/\$contactEmail/$serviceDetails->{'contactEmail'}/g;
	$html =~ s/\$contactAddress/$serviceDetails->{'contactAddress'}/g;
	$html =~ s/\$submitDate/$serviceDetails->{'submitDate'}/g;
	$html =~ s/\$submitPerson/$serviceDetails->{'submitPerson'}/g;
	$html =~ s/\$startDate/$serviceDetails->{'startDate'}/g;
	$html =~ s/\$endDate/$serviceDetails->{'endDate'}/g;
	$html =~ s/\$comments/$serviceDetails->{'comments'}/g;
	$html =~ s/\$status/$status{$service[7]}/g;
	$html =~ s/\$projectId/$service[6]/g;
	$html =~ s/\$projectName/$project[2]/g;
	$html =~ s/\$serviceCreator/$service[9]/g;
	$html =~ s/\$serviceEnteredDate/$service[10]/g;
	$html =~ s/\$\$/$$/g;

	print $html;
}
else
{
	print '402 Invalid operation';
	exit;
}

__DATA__
	<b>For Project</b> $projectName <sup class='ui-state-disabled'>(Project id: $projectId)</sup><br>
	<label for="viewServiceName"><b>Service Code</b></label> $serviceName
	<sup class='ui-state-disabled'>Last changed by $serviceCreator on $serviceEnteredDate</sup>
	<div id="viewServiceTabs$$">
		<ul>
		<li><a href="#viewServiceTabs-1">Submitter Information</a></li>
		<li><a href="#viewServiceTabs-2">Type of Service</a></li>
		<li><a href="#viewServiceTabs-3">Sample Information</a></li>
		<li><a href="#viewServiceTabs-4">Status</a></li>
		</ul>
		<div id="viewServiceTabs-1">
			<table>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="viewServicePiName"><b>P.I. Name</b></label></td>
				<td>$piName</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="viewServicePiEmail"><b>P.I. Email</b></label></td>
				<td>$piEmail</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="viewServiceInstitution"><b>Institution</b></label></td>
				<td>$institution</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="viewServiceContactName"><b>Contact Person's Name</b></label></td>
				<td>$contactName</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="viewServiceContactPhone"><b>Contact Person's Phone</b></label></td>
				<td>$contactPhone</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="viewServiceContactEmail"><b>Contact Person's Email</b></label></td>
				<td>$contactEmail</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="viewServiceContactAddress"><b>Contact Person's Address</b></label></td>
				<td>$contactAddress</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="viewServiceSubmitDate"><b>Date Submitted</b></label></td>
				<td>$submitDate</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="viewServiceSubmitPerson"><b>Submitted Person's Name</b></label></td>
				<td>$submitPerson</td>
			</tr>
			</table>
		</div>
		<div id="viewServiceTabs-2">
			<table>
			<tr><td style='text-align:right;white-space: nowrap;'><b>Service Requested</b></td>
				<td>
					$serviceType
					<br>$serviceTypeOther
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="viewServiceGenus"><b>Genus Species Subspecies</b></label><br>(Common Name)</td>
				<td>$genus $species $subspecies<br>$commonName
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="viewServiceGenomeSize"><b>Genome size</b></label></td>
				<td>$genomeSize</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><b>Notes about service desired</b></td>
				<td>
					<label for="viewServiceAmount">Amount or coverage:</label> $amount<br>
					<label for="viewServiceAssembly">Assembly:</label>$assembly<br>
					<label for="viewServiceAnalysis">Analysis:</label>$analysis<br>
					<label for="viewServiceNote">Other:</label><br>$serviceNote
				</td>
			</tr>
			</table>
		</div>
		<div id="viewServiceTabs-3">
			<table>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="viewServiceSamples"><b>Related Samples</b></label></td>
				<td><ul style='margin: 0;'>$relatedSamples</ul></td>
			</tr>
			</table>
		</div>
		<div id="viewServiceTabs-4">
			<table>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="viewServiceStatus"><b>Service Status</b></label></td>
				<td>$status
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="viewServiceStartDate"><b>Start Date</b></label></td>
				<td>$startDate</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="viewServiceEndDate"><b>End Date</b></label></td>
				<td>$endDate</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><b>Payment type</b></td>
				<td>
					$paymentType
					$poNumber
					($paymentStatus$paymentDate)
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="viewServiceInvoice"><b>Invoice number</b></label></td>
				<td>$invoice</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="viewServiceComments"><b>Comments</b></label></td>
				<td>$comments</td>
			</tr>
			</table>
		</div>
	</div>
<script>
$( "#viewServiceTabs$$" ).tabs();
$('#dialog').dialog("option", "title", "View Service $serviceName");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Edit", click: function() { closeDialog();openDialog("serviceEdit.cgi?serviceId=$serviceId"); } }, { text: "OK", click: function() {closeDialog(); } } ] );
</script>
