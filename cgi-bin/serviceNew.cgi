#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use DBI;
use POSIX qw(strftime);
use lib "lib/";
use lib "lib/pangu";
use pangu;
use userCookie;

my $userCookie = new userCookie;
my $userId = (cookie('cid')) ? $userCookie->checkCookie(cookie('cid')) : 0;
exit if (!$userId);

my $commoncfg = readConfig("main.conf");
my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});

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
	2=>'Credit card, 3% fee added',
	3=>'Purchase order',
	4=>'Wire transfer, fee of $50 added'
	);

print header;
my $projectId = param ('projectId') || '';
if($projectId)
{
	my $project=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$project->execute($projectId);
	my @project = $project->fetchrow_array();
	my $serviceTypeList = '';
	foreach (sort {$a <=> $b} keys %serviceType)
	{
		next unless ($_);
		$serviceTypeList .= "<option value='$_'>$_. $serviceType{$_}</option>";
	}
	$serviceTypeList .= "<option value='0'>$serviceType{0}</option>";

	my $paymentTypeList = '';
	foreach (sort {$a <=> $b} keys %paymentType)
	{
		$paymentTypeList .= "<option value='$_'>$paymentType{$_}</option>";
	}
	my $datetimeString = strftime "%y %Y %m %d %H:%M:%S", localtime;
	my ($tdYear,$year,$month,$dayofmonth,$time)=split(/\s+/,$datetimeString);
	my $autoStartDate = "$year-$month-$dayofmonth";

	undef $/;# enable slurp mode
	my $html = <DATA>;
	$html =~ s/\$projectId/$projectId/g;
	$html =~ s/\$projectName/$project[2]/g;
	$html =~ s/\$serviceType/$serviceTypeList/g;
	$html =~ s/\$autoStartDate/$autoStartDate/g;
	$html =~ s/\$paymentType/$paymentTypeList/g;

	print $html;
}
else
{
	print <<END;
<script>
	errorPop("Not a valid operation: no projectId provided!");
</script>	
END
	exit;
}

__DATA__
<form id="newService" name="newService" action="serviceSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<label for="newServiceProjectId"><b>For Project $projectName</b></label><input name="projectId" id="newServiceProjectId" type="hidden" value="$projectId" /><br>
	<label for="newServiceName"><b>Service Code</b></label><input class='ui-widget-content ui-corner-all' name="name" id="newServiceName" placeholder="Service Coce" size="15" type="text" maxlength="9" />
	<sub class='ui-state-highlight'>This Code will NOT be editable after saving the form.</sub>
	<div id="newServiceTabs">
		<ul>
		<li><a href="#newServiceTabs-1">Submitter Information</a></li>
		<li><a href="#newServiceTabs-2">Type of Service</a></li>
		<li><a href="#newServiceTabs-3">Status</a></li>
		</ul>
		<div id="newServiceTabs-1">
			<table>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newServicePiName"><b>P.I. Name</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="piName" id="newServicePiName" type="text" value="" placeholder="P.I. Name" size="25" type="text" maxlength="255" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newServicePiEmail"><b>P.I. Email</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="piEmail" id="newServicePiEmail" type="text" value="" placeholder="P.I. Email" size="45" type="text" maxlength="255" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newServiceInstitution"><b>Institution</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="institution" id="newServiceInstitution" type="text" value="" placeholder="Institution" size="45" type="text" maxlength="255" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newServiceContactName"><b>Contact Person's Name</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="contactName" id="newServiceContactName" type="text" value="" placeholder="Contact Person" size="25" type="text" maxlength="255" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newServiceContactPhone"><b>Contact Person's Phone</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="contactPhone" id="newServiceContactPhone" type="text" value="" placeholder="Contact Phone" size="15" type="text" maxlength="25" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newServiceContactEmail"><b>Contact Person's Email</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="contactEmail" id="newServiceContactEmail" type="text" value="" placeholder="Contact Email" size="45" type="text" maxlength="255" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newServiceContactAddress"><b>Contact Person's Address</b></label></td>
				<td><textarea class='ui-widget-content ui-corner-all' name="contactAddress" id="newServiceContactAddress" cols="40" rows="3" placeholder="Address"></textarea></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newServiceSubmitDate"><b>Date Submitted</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="submitDate" id="newServiceSubmitDate" placeholder="YYYY-MM-DD" size="20" type="text" maxlength="19" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newServiceSubmitPerson"><b>Submitted Person's Name</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="submitPerson" id="newServiceSubmitPerson" placeholder="Submitted Person" size="25" type="text" maxlength="255" /></td>
			</tr>
			</table>
		</div>
		<div id="newServiceTabs-2">
			<table>
			<tr><td style='text-align:right;white-space: nowrap;'><b>Service Requested</b></td>
				<td>
					<select class='ui-widget-content ui-corner-all' name="serviceType" id="newServiceType" onchange="hideShowOther();">$serviceType</select>
					<br><textarea class='ui-widget-content ui-corner-all' name="serviceTypeOther" id="newServiceTypeOther" cols="40" rows="2" placeholder=""></textarea>
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newServiceGenus"><b>Genus Species Subspecies</b></label><br>(Common Name)</td>
				<td><INPUT class="ui-widget-content ui-corner-all" TYPE="text" name="genus" id="newServiceGenus" placeholder="Genus" maxlength="45" size="12">
					<INPUT class="ui-widget-content ui-corner-all" TYPE="text" name="species" id="newServiceSpecies" placeholder="Species" maxlength="45" size="12">
					<INPUT class="ui-widget-content ui-corner-all" TYPE="text" name="subspecies" id="newServiceSubspecies" placeholder="Subspecies" maxlength="45" size="12">
					<INPUT class="ui-widget-content ui-corner-all" TYPE="text" name="commonName" id="newServiceCommonName" placeholder="Common Name" maxlength="45" size="40">
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newServiceGenomeSize"><b>Genome size</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="genomeSize" id="newServiceGenomeSize" placeholder="Genome Size" size="20" type="text" maxlength="45" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><b>Notes about service desired</b></td>
				<td>
					<label for="newServiceAmount">Amount or coverage:</label><input class='ui-widget-content ui-corner-all' name="amount" id="newServiceAmount" placeholder="" size="25" type="text" maxlength="255" /><br>
					<label for="newServiceAssembly">Assembly:</label><input class='ui-widget-content ui-corner-all' name="assembly" id="newServiceAssembly" placeholder="" size="25" type="text" maxlength="255" /><br>
					<label for="newServiceAnalysis">Analysis:</label><input class='ui-widget-content ui-corner-all' name="analysis" id="newServiceAnalysis" placeholder="" size="25" type="text" maxlength="255" /><br>
					<label for="newServiceNote">Other:</label><br><textarea class='ui-widget-content ui-corner-all' name="serviceNote" id="newServiceNote" cols="40" rows="3" placeholder=""></textarea>
				</td>
			</tr>
			</table>
		</div>
		<div id="newServiceTabs-3">
			<table>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newServiceStatus"><b>Status</b></label></td>
				<td><div id="newServiceStatus">
					<input type="radio" id="newServiceStatusRadio1" name="status" value="1" checked="checked"><label for="newServiceStatusRadio1">Active</label>
					<input type="radio" id="newServiceStatusRadio2" name="status" value="0"><label for="newServiceStatusRadio2">Inactive</label>
					</div>
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newServiceStartDate"><b>Start Date</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="startDate" id="newServiceStartDate" value="$autoStartDate" placeholder="YYYY-MM-DD" size="20" type="text" maxlength="19" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newServiceEndDate"><b>End Date</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="endDate" id="newServiceEndDate" placeholder="YYYY-MM-DD" size="20" type="text" maxlength="19" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><b>Payment type</b></td>
				<td>
					<select class='ui-widget-content ui-corner-all' name="paymentType" id="newServicePaymentType" onchange="hideShowOther();">$paymentType</select><input class='ui-widget-content ui-corner-all' name="poNumber" id="newServicePoNumber" type="text" value="" placeholder="PO#" size="15" type="text" maxlength="25" />
					<br><input type="checkbox" id="newServicePaymentStatus" name="paymentStatus" value="1" onchange="hideShowOther();"><label for="newServicePaymentStatus">Fully paid</label>
					<input class='ui-widget-content ui-corner-all' name="paymentDate" id="newServicePaymentDate" placeholder="Paid on YYYY-MM-DD" size="20" type="text" maxlength="19" />
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newServiceInvoice"><b>Invoice number</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="invoice" id="newServiceInvoice" placeholder="Invoice#" size="20" type="text" maxlength="256" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newServiceComments"><b>Comments</b></label></td>
				<td><textarea class='ui-widget-content ui-corner-all' name="comments" id="newServiceComments" cols="40" rows="5" placeholder="Please write anything relevant to the service that may not be included on this form"></textarea></td>
			</tr>
			</table>
		</div>
	</div>
</form>
<script>
hideShowOther();
$( "#newServiceTabs" ).tabs();
$( "#newServiceStatus" ).buttonset();
$( "#newServiceSubmitDate" ).datepicker({dateFormat:"yy-mm-dd"});
$( "#newServiceStartDate" ).datepicker({
	dateFormat:"yy-mm-dd",
	onClose: function( selectedDate ) {
        $( "#newServiceEndDate" ).datepicker( "option", "minDate", selectedDate );
    }
});
$( "#newServiceEndDate" ).datepicker({
	dateFormat:"yy-mm-dd",
	onClose: function( selectedDate ) {
        $( "#newServiceStartDate" ).datepicker( "option", "maxDate", selectedDate );
    }
});
$( "#newServicePaymentDate" ).datepicker({dateFormat:"yy-mm-dd"});
$('#dialog').dialog("option", "title", "New Service");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('newService'); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
function hideShowOther()
{
	if($('#newServiceType').val() == 0)
	{
		$('#newServiceTypeOther').show();
	}
	else
	{
		$('#newServiceTypeOther').hide();
	}

	if($('#newServicePaymentType').val() == 3)
	{
		$('#newServicePoNumber').show();
	}
	else
	{
		$('#newServicePoNumber').hide();
	}

	if($('#newServicePaymentStatus').prop('checked'))
	{
		$('#newServicePaymentDate').show();
	}
	else
	{
		$('#newServicePaymentDate').hide();
	}
}
</script>
