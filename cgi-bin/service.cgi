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

my $button;
my $samples;
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

	$button = "<ul id='serviceInfoMenu$serviceId$$' style='margin-top: .3em;width: 250px;'>
				<li><a><span class='ui-icon ui-icon-triangle-1-e'></span><b>Service '$service[2]'</b></a>
					<ul style='z-index: 1000;white-space: nowrap;'>
						<li><a onclick='openDialog(\"serviceEdit.cgi?serviceId=$serviceId\")' title='Edit/Delete $service[2]'><span class='ui-icon ui-icon-pencil'></span>Edit/Delete</a></li>
						<li><a onclick='openDialog(\"sampleNew.cgi?serviceId=$serviceId\")' title='New Sample'>New Sample</a></li>
						<li><a onclick='openDialogForm(\"itemDepotForm.cgi\",\"viewService$serviceId\")' title='Depot Item'>Depot Item</a></li>
						<li><a onclick='openDialogForm(\"itemDeleteForm.cgi\",\"viewService$serviceId\")' title='Delete Item'><span class='ui-icon ui-icon-trash'></span>Delete Item</a></li>
						<li><a title='Print'><span class='ui-icon ui-icon-print'></span>Print</a>
							<ul style='z-index: 1000;white-space: nowrap;'>
								<li><a onclick='openDialogForm(\"barcodePrintForm.cgi\",\"viewService$serviceId\")' title='Print Labels'><span class='ui-icon ui-icon-bullet'></span>Labels</a></li>
								<li><a onclick='printDiv(\"viewService$serviceId\")' title='Print Service'><span class='ui-icon ui-icon-bullet'></span>Service</a></li>
							</ul>
						</li>
					</ul>
				</li>
			</ul>";
	
	undef $/;# enable slurp mode
	my $html = <DATA>;
	$html =~ s/\$button/$button/g;

	my $sampleInService=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sample' AND z = ? ORDER BY o");
	$sampleInService->execute($serviceId);
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
		my %sampleStatus = (
			0=>'Na',
			1=>'Status One',
			2=>'Status Two'
			);
		my $sampleDetails = decode_json $sampleInService[8];
		$sampleDetails->{'sampleTypeOther'} = '' unless (exists $sampleDetails->{'sampleTypeOther'});
		$sampleDetails->{'description'} = '' unless (exists $sampleDetails->{'description'});
		$sampleDetails->{'description'} = escapeHTML($sampleDetails->{'description'});
		$sampleDetails->{'description'} =~ s/\n/<br>/g;
		$sampleInService[4] = ($sampleInService[4] == 0 || $sampleInService[4] == 7) ? "$sampleType{$sampleInService[4]}$sampleDetails->{'sampleTypeOther'}" : "$sampleType{$sampleInService[4]}";
		$samples = "<table id='sampleInService$serviceId$$' name='sampleInService$serviceId$$' class='display'>
			<thead><tr><th><input type='checkbox' id='checkAllSample$$' name='checkAllSample$$' value='Check all' checked='checked' onClick='checkClass(\"sampleAll\");return false;' title='Check all samples'>
			<input type='checkbox' id='uncheckAllSample$$' name='uncheckAllSample$$' value='Uncheck all' onClick='uncheckClass(\"sampleAll\");return false;' title='Uncheck all samples'></th>
			<th style='text-align:left'><b>Samples</b></th>
			<th><b>Details</b></th>
			</tr></thead><tbody>" unless ($samples);
		my $sampleDescription= '';
		if($sampleInService[3] == 1) #for paclib
		{
			my $newPaclib = "<a style='float: center;' class='ui-state-highlight ui-corner-all' onclick='openDialog(\"paclibNew.cgi?sampleId=$sampleInService[0]\")' title='New Paclib'>
				<span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-plus' title='New Paclib'></span></a>";
			my $relatedPaclibs = '';
			my $paclibForSample=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'paclib' AND z = ? ORDER BY o");
			$paclibForSample->execute($sampleInService[0]);
			while (my @paclibForSample = $paclibForSample->fetchrow_array())
			{
				my $paclibDetails = decode_json $paclibForSample[8];
				$paclibDetails->{'libraryDate'} = '' unless (exists $paclibDetails->{'libraryDate'});
				$paclibDetails->{'shearingInput'} = '' unless (exists $paclibDetails->{'shearingInput'});
				$paclibDetails->{'shearingRpm'} = '' unless (exists $paclibDetails->{'shearingRpm'});
				$paclibDetails->{'shearingOutput'} = '' unless (exists $paclibDetails->{'shearingOutput'});
				$paclibDetails->{'shearingBeadsSteps'} = '' unless (exists $paclibDetails->{'shearingBeadsSteps'});
				$paclibDetails->{'bluepippinInput'} = '' unless (exists $paclibDetails->{'bluepippinInput'});
				$paclibDetails->{'bluepippinSize'} = '' unless (exists $paclibDetails->{'bluepippinSize'});
				$paclibDetails->{'bluepippinOutput'} = '' unless (exists $paclibDetails->{'bluepippinOutput'});
				$paclibDetails->{'bluepippinConcentration'} = '' unless (exists $paclibDetails->{'bluepippinConcentration'});
				$paclibDetails->{'description'} = '' unless (exists $paclibDetails->{'description'});
				my %smrtrunStatus = (
				0=>'Initialized',
				1=>'Started',
				2=>'Completed'
				);

				my $newSmrtwell = "<a style='float: right;' class='ui-state-highlight ui-corner-all' onclick='openDialog(\"smrtwellNew.cgi?paclibId=$paclibForSample[0]\")' title='New SMRT Well'>
				<span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-circle-plus' title='New SMRT Well'></span></a>";
				my $relatedSmrtwells = '';
				my $smrtwellForPaclib=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'smrtwell' AND o = ? ORDER BY name");
				$smrtwellForPaclib->execute($paclibForSample[0]);
				while (my @smrtwellForPaclib = $smrtwellForPaclib->fetchrow_array())
				{
					my $smrtrun=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
					$smrtrun->execute($smrtwellForPaclib[6]);
					my @smrtrun = $smrtrun->fetchrow_array();
					$relatedSmrtwells .= ($relatedSmrtwells) ? "<br><a onclick='closeDialog();openDialog(\"smrtwellView.cgi?smrtwellId=$smrtwellForPaclib[0]\")' title='SMRT Run ($smrtrunStatus{$smrtrun[7]})'>$smrtrun[2]>$smrtwellForPaclib[2]</a>"
					: "<a onclick='closeDialog();openDialog(\"smrtwellView.cgi?smrtwellId=$smrtwellForPaclib[0]\")' title='SMRT Run ($smrtrunStatus{$smrtrun[7]})'>$smrtrun[2]>$smrtwellForPaclib[2]</a>";	
				}
				$relatedPaclibs .= "<tr>
					<td style='text-align:center'><input type='checkbox' class='paclibAll$sampleInService[0]' id='paclibList$sampleInService[0]$serviceId$$' name='itemId' value='$paclibForSample[0]'></td>
					<td style='text-align:center'><a onclick='openDialog(\"paclibView.cgi?paclibId=$paclibForSample[0]\")' title='View'>$paclibForSample[2]</a></td>
					<td style='text-align:center' title='Shearing DNA input'>$paclibDetails->{'shearingInput'} ug</td>
					<td style='text-align:center' title='Shearing RPM'>$paclibDetails->{'shearingRpm'}</td>
					<td style='text-align:center' title='Shearing DNA output'>$paclibDetails->{'shearingOutput'} ug</td>
					<td style='text-align:center' title='BluePippin DNA input'>$paclibDetails->{'bluepippinInput'} ug</td>
					<td style='text-align:center' title='BluePippin Size selection'>$paclibDetails->{'bluepippinSize'} kb</td>
					<td style='text-align:center' title='BluePippin DNA output'>$paclibDetails->{'bluepippinOutput'} ug</td>
					<td style='text-align:center' title='BluePippin DNA concentration'>$paclibDetails->{'bluepippinConcentration'} ng/uL</td>
					<td style='text-align:left' title='Related SMRT Wells'>$newSmrtwell$relatedSmrtwells</td>
					</tr>";	
			}
			$sampleDescription = "<table style='width:100%'>
					<tr>
						<th rowspan='2'><input type='checkbox' id='checkAllPaclib$sampleInService[0]$$' name='checkAllPaclib$sampleInService[0]$$' value='Check all' checked='checked' onClick='checkClass(\"paclibAll$sampleInService[0]\");return false;' title='Check all paclibs'>
							<input type='checkbox' id='uncheckAllPaclib$sampleInService[0]$$' name='uncheckAllPaclib$sampleInService[0]$$' value='Uncheck all' onClick='uncheckClass(\"paclibAll$sampleInService[0]\");return false;' title='Uncheck all paclibs'></th>
						<th>Name</th>
						<th colspan='3'>Shearing</th>
						<th colspan='4'>BluePippin</th>
						<th rowspan='2'>SMRT Sample Plate > Well</th>
					</tr>
					<tr>
						<th>$newPaclib</th>
						<th><sup>DNA input</sup></th>
						<th><sup>RPM</sup></th>
						<th><sup>DNA output</sup></th>
						<th><sup>DNA input</sup></th>
						<th><sup>Size selection</sup></th>
						<th><sup>DNA output</sup></th>
						<th><sup>DNA concentration</sup></th>
					</tr>
					$relatedPaclibs
				</table>";
		}
		else
		{
			$sampleDescription = $sampleDetails->{'description'};
		}

		$samples .= "<tr>
			<td style='text-align:center'><input type='checkbox' class='sampleAll' id='sampleList$serviceId$$' name='itemId' value='$sampleInService[0]'></td>
			<td style='text-align:left;'>
				<a onclick='openDialog(\"sampleView.cgi?sampleId=$sampleInService[0]\")' title='$sampleDetails->{'description'}'>$sampleInService[2]<br>
				<img alt='$sampleInService[7]' src='barcode.cgi?code=$sampleInService[7]' title='Barcode'/></a><br>
				<span title='Type'a>$sampleInService[4]</span><br>
				Status: $sampleStatus{$sampleInService[5]}
				</td>
			<td>$sampleDescription
			</td>
			</tr>";	
	}

	if ($samples)
	{
		$samples .= "</tbody></table>" 
	}
	else
	{
		$samples = "<p class='ui-state-error ui-corner-all' style='padding: .7em;'><span class='ui-icon ui-icon-alert' style='float: left; margin-right: .3em;'></span><strong>No sample found.</strong></p>";
	}
	$html =~ s/\$samples/$samples/g;
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
		$html =~ s/\$paymentStatus/Fully paid on/g;
		$html =~ s/\$paymentDate/$serviceDetails->{'paymentDate'}/g;
	}
	else
	{
		$html =~ s/\$paymentStatus/Not paid./g;
		$html =~ s/\$paymentDate//g;
	}
	$html =~ s/\$paymentDate/$serviceDetails->{'paymentDate'}/g;
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
	$html =~ s/\$serviceCreator/$service[9]/g;
	$html =~ s/\$serviceDateEntered/$service[10]/g;
	$html =~ s/\$\$/$$/g;

	print header(-cookie=>[cookie(-name=>'service',-value=>$serviceId),cookie(-name=>'library',-value=>0)]);
	print $html;
}
else
{
	print header(-type=>'text/html',-status=>'402 Invalid operation');
	exit;
}

__DATA__
$button
<form id='viewService$serviceId' name='viewService$serviceId'>
	<div id="viewServiceTabs$serviceId">
		<ul>
		<li><a href="#viewServiceTabs-1">Sample Details</a></li>
		<li><a href="#viewServiceTabs-2">About $serviceName</a></li>
		</ul>
		<div id="viewServiceTabs-1">
		$samples
		</div>
		<div id="viewServiceTabs-2">
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
			<tr><td style='text-align:right;white-space: nowrap;'><b>Invoice number</b></td>
				<td>
					$invoice
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="viewServiceComments"><b>Comments</b></label></td>
				<td>$comments</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'></td>
				<td><sup class='ui-state-disabled'>Last changed by $serviceCreator on $serviceDateEntered</sup></td>
			</tr>
			</table>
		</div>
	</div>
</form>
<script>
buttonInit();
$( "#viewServiceTabs$serviceId" ).tabs();
$( "#serviceInfoMenu$serviceId$$" ).menu();
$( "#sampleInService$serviceId$$" ).dataTable({
	"scrollCollapse": true,
	"paging": false
});
loadingHide();
</script>