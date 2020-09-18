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

my $serviceId = param ('serviceId') || '';
my $service = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$service->execute($serviceId);
my @service=$service->fetchrow_array();
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

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$serviceId/$serviceId/g;
$html =~ s/\$serviceName/$service[2]/g;
my $serviceTypeList = '';
foreach (sort {$a <=> $b} keys %serviceType)
{
	next unless ($_);
	$serviceTypeList .= ($service[3] eq $_) ? "<option value='$_' selected>$_. $serviceType{$_}</option>" :  "<option value='$_'>$_. $serviceType{$_}</option>";
}
$serviceTypeList .= ($service[3] == 0) ? "<option value='0' selected>$serviceType{0}</option>" : "<option value='0'>$serviceType{0}</option>";
$html =~ s/\$serviceTypeOther/$serviceDetails->{'serviceTypeOther'}/g;
$html =~ s/\$serviceType/$serviceTypeList/g;
$html =~ s/\$amount/$serviceDetails->{'amount'}/g;
$html =~ s/\$assembly/$serviceDetails->{'assembly'}/g;
$html =~ s/\$analysis/$serviceDetails->{'analysis'}/g;
$html =~ s/\$serviceNote/$serviceDetails->{'serviceNote'}/g;
$html =~ s/\$genus/$serviceDetails->{'genus'}/g;
$html =~ s/\$species/$serviceDetails->{'species'}/g;
$html =~ s/\$subspecies/$serviceDetails->{'subspecies'}/g;
$html =~ s/\$commonName/$serviceDetails->{'commonName'}/g;
$html =~ s/\$genomeSize/$serviceDetails->{'genomeSize'}/g;
my $paymentTypeList = '';
foreach (sort {$a <=> $b} keys %paymentType)
{
	$paymentTypeList .= ($serviceDetails->{'paymentType'} == $_) ? "<option value='$_' selected>$paymentType{$_}</option>" : "<option value='$_'>$paymentType{$_}</option>";
}
$html =~ s/\$paymentType/$paymentTypeList/g;
$html =~ s/\$poNumber/$serviceDetails->{'poNumber'}/g;
if($service[5])
{
	$html =~ s/\$editServicePaymentStatus/checked="checked"/g;
}
else
{
	$html =~ s/\$editServicePaymentStatus//g;
}
$html =~ s/\$paymentDate/$serviceDetails->{'paymentDate'}/g;
$html =~ s/\$invoice/$serviceDetails->{'invoice'}/g;
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
if($service[7])
{
	$html =~ s/\$editServiceStatusRadio1/checked="checked"/g;
	$html =~ s/\$editServiceStatusRadio2//g;
}
else
{
	$html =~ s/\$editServiceStatusRadio1//g;
	$html =~ s/\$editServiceStatusRadio2/checked="checked"/g;
}
$html =~ s/\$projectId/$service[6]/g;
my $project=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$project->execute($service[6]);
my @project = $project->fetchrow_array();
$html =~ s/\$projectName/$project[2]/g;

$html =~ s/\$serviceCreator/$service[9]/g;
$html =~ s/\$serviceEnteredDate/$service[10]/g;

print header;
print $html;

__DATA__
<form id="editService" name="editService" action="serviceSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="serviceId" id="editServiceId" type="hidden" value="$serviceId" />
	<label for="editServiceProjectId"><b>For Project $projectName</b></label><input name="projectId" id="editServiceProjectId" type="hidden" value="$projectId" /><br>
	<label for="editServiceName"><b>Service Code</b></label><input class='ui-widget-content ui-corner-all' name="name" id="editServiceName" placeholder="Service Name" size="15" type="text" maxlength="9" value="$serviceName" readonly />
	<sup class='ui-state-disabled'>Last changed by $serviceCreator on $serviceEnteredDate</sup>
	<div id="editServiceTabs">
		<ul>
		<li><a href="#editServiceTabs-1">Submitter Information</a></li>
		<li><a href="#editServiceTabs-2">Type of Service</a></li>
		<li><a href="#editServiceTabs-3">Status</a></li>
		</ul>
		<div id="editServiceTabs-1">
			<table>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editServicePiName"><b>P.I. Name</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="piName" id="editServicePiName" type="text" value="$piName" placeholder="P.I. Name" size="25" type="text" maxlength="255" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editServicePiEmail"><b>P.I. Email</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="piEmail" id="editServicePiEmail" type="text" value="$piEmail" placeholder="P.I. Email" size="45" type="text" maxlength="255" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editServiceInstitution"><b>Institution</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="institution" id="editServiceInstitution" type="text" value="$institution" placeholder="Institution" size="45" type="text" maxlength="255" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editServiceContactName"><b>Contact Person's Name</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="contactName" id="editServiceContactName" type="text" value="$contactName" placeholder="Contact Person" size="25" type="text" maxlength="255" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editServiceContactPhone"><b>Contact Person's Phone</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="contactPhone" id="editServiceContactPhone" type="text" value="$contactPhone" placeholder="Contact Phone" size="15" type="text" maxlength="25" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editServiceContactEmail"><b>Contact Person's Email</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="contactEmail" id="editServiceContactEmail" type="text" value="$contactEmail" placeholder="Contact Email" size="45" type="text" maxlength="255" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editServiceContactAddress"><b>Contact Person's Address</b></label></td>
				<td><textarea class='ui-widget-content ui-corner-all' name="contactAddress" id="editServiceContactAddress" cols="40" rows="3" placeholder="Address">$contactAddress</textarea></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editServiceSubmitDate"><b>Date Submitted</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="submitDate" id="editServiceSubmitDate" value="$submitDate" placeholder="YYYY-MM-DD" size="20" type="text" maxlength="19" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editServiceSubmitPerson"><b>Submitted Person's Name</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="submitPerson" id="editServiceSubmitPerson" value="$submitPerson" placeholder="Submitted Person" size="25" type="text" maxlength="255" /></td>
			</tr>
			</table>
		</div>
		<div id="editServiceTabs-2">
			<table>
			<tr><td style='text-align:right;white-space: nowrap;'><b>Service Requested</b></td>
				<td>
					<select class='ui-widget-content ui-corner-all' name="serviceType" id="editServiceType" onchange="hideShowOther();">$serviceType</select>
					<br><textarea class='ui-widget-content ui-corner-all' name="serviceTypeOther" id="editServiceTypeOther" cols="40" rows="2" placeholder="">$serviceTypeOther</textarea>
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editServiceGenus"><b>Genus Species Subspecies</b></label><br>(Common Name)</td>
				<td><INPUT class="ui-widget-content ui-corner-all" TYPE="text" name="genus" id="editServiceGenus" value="$genus" placeholder="Genus" maxlength="45" size="12">
					<INPUT class="ui-widget-content ui-corner-all" TYPE="text" name="species" id="editServiceSpecies" value="$species" placeholder="Species" maxlength="45" size="12">
					<INPUT class="ui-widget-content ui-corner-all" TYPE="text" name="subspecies" id="editServiceSubspecies" value="$subspecies" placeholder="Subspecies" maxlength="45" size="12">
					<INPUT class="ui-widget-content ui-corner-all" TYPE="text" name="commonName" id="editServiceCommonName" value="$commonName" placeholder="Common Name" maxlength="45" size="40">
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editServiceGenomeSize"><b>Genome size</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="genomeSize" id="editServiceGenomeSize" value="$genomeSize" placeholder="Genome Size" size="20" type="text" maxlength="45" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><b>Notes about service desired</b></td>
				<td>
					<label for="editServiceAmount">Amount or coverage:</label><input class='ui-widget-content ui-corner-all' name="amount" id="editServiceAmount" value="$amount" placeholder="" size="25" type="text" maxlength="255" /><br>
					<label for="editServiceAssembly">Assembly:</label><input class='ui-widget-content ui-corner-all' name="assembly" id="editServiceAssembly" value="$assembly" placeholder="" size="25" type="text" maxlength="255" /><br>
					<label for="editServiceAnalysis">Analysis:</label><input class='ui-widget-content ui-corner-all' name="analysis" id="editServiceAnalysis" value="$analysis" placeholder="" size="25" type="text" maxlength="255" /><br>
					<label for="editServiceNote">Other:</label><br><textarea class='ui-widget-content ui-corner-all' name="serviceNote" id="editServiceNote" cols="40" rows="3" placeholder="">$serviceNote</textarea>
				</td>
			</tr>
			</table>
		</div>
		<div id="editServiceTabs-3">
			<table>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editServiceStatus"><b>Service Status</b></label></td>
				<td><div id="editServiceStatus">
					<input type="radio" id="editServiceStatusRadio1" name="status" value="1" $editServiceStatusRadio1><label for="editServiceStatusRadio1">Active</label>
					<input type="radio" id="editServiceStatusRadio2" name="status" value="0" $editServiceStatusRadio2><label for="editServiceStatusRadio2">Inactive</label>
					</div>
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editServiceStartDate"><b>Start Date</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="startDate" id="editServiceStartDate" value="$startDate" placeholder="YYYY-MM-DD" size="20" type="text" maxlength="19" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editServiceEndDate"><b>End Date</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="endDate" id="editServiceEndDate" value="$endDate" placeholder="YYYY-MM-DD" size="20" type="text" maxlength="19" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><b>Payment type</b></td>
				<td>
					<select class='ui-widget-content ui-corner-all' name="paymentType" id="editServicePaymentType" onchange="hideShowOther();">$paymentType</select><input class='ui-widget-content ui-corner-all' name="poNumber" id="editServicePoNumber" type="text" value="$poNumber" placeholder="PO#" size="15" type="text" maxlength="25" />
					<br><input type="checkbox" id="editServicePaymentStatus" name="paymentStatus" value="1" onchange="hideShowOther();" $editServicePaymentStatus><label for="editServicePaymentStatus">Fully paid</label>
					<input class='ui-widget-content ui-corner-all' name="paymentDate" id="editServicePaymentDate" value="$paymentDate" placeholder="Paid on YYYY-MM-DD" size="20" type="text" maxlength="19" />
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editServiceInvoice"><b>Invoice number</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="invoice" id="editServiceInvoice" value="$invoice" placeholder="Invoice#" size="20" type="text" maxlength="256" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editServiceComments"><b>Comments</b></label></td>
				<td><textarea class='ui-widget-content ui-corner-all' name="comments" id="editServiceComments" cols="40" rows="5" placeholder="Please write anything relevant to the service that may not be included on this form">$comments</textarea></td>
			</tr>
			</table>
		</div>
	</div>
</form>
<script>
hideShowOther();
$( "#editServiceTabs" ).tabs();
$( "#editServiceStatus" ).buttonset();
$( "#editServiceSubmitDate" ).datepicker({dateFormat:"yy-mm-dd"});
$( "#editServiceStartDate" ).datepicker({
	dateFormat:"yy-mm-dd",
	onClose: function( selectedDate ) {
        $( "#editServiceEndDate" ).datepicker( "option", "minDate", selectedDate );
    }
});
$( "#editServiceEndDate" ).datepicker({
	dateFormat:"yy-mm-dd",
	onClose: function( selectedDate ) {
        $( "#editServiceStartDate" ).datepicker( "option", "maxDate", selectedDate );
    }
});
$( "#editServicePaymentDate" ).datepicker({dateFormat:"yy-mm-dd"});
$('#dialog').dialog("option", "title", "Edit Service $serviceName");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('editService'); } }, { text: "Delete", click: function() { deleteItem($serviceId); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
function hideShowOther()
{
	if($('#editServiceType').val() == 0)
	{
		$('#editServiceTypeOther').show();
	}
	else
	{
		$('#editServiceTypeOther').hide();
	}

	if($('#editServicePaymentType').val() == 3)
	{
		$('#editServicePoNumber').show();
	}
	else
	{
		$('#editServicePoNumber').hide();
	}

	if($('#editServicePaymentStatus').prop('checked'))
	{
		$('#editServicePaymentDate').show();
	}
	else
	{
		$('#editServicePaymentDate').hide();
	}
}
</script>


