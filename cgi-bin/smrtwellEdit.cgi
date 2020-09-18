#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use JSON::XS;
use DBI;
use POSIX qw(strftime);
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

my $smrtwellId = param ('smrtwellId') || '';
my %wellx = ( 96 => ['01' .. '12'], 384 =>['01' .. '24']);
my %welly = ( 96 => ['A' .. 'H'], 384 =>['A' .. 'P']);
my %movieLength = (
	0=>'Customized:',
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
my %smrtrunStatus = (
	0=>'Initialized',
	1=>'Started',
	2=>'Completed'
	);
my $smrtwell = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$smrtwell->execute($smrtwellId);
my @smrtwell=$smrtwell->fetchrow_array();
my $smrtwellDetails = decode_json $smrtwell[8];
$smrtwellDetails->{'loadingName'} = '' unless (exists $smrtwellDetails->{'loadingName'});
$smrtwellDetails->{'concentration'} = '' unless (exists $smrtwellDetails->{'concentration'});
$smrtwellDetails->{'polRatio'} = '' unless (exists $smrtwellDetails->{'polRatio'});
$smrtwellDetails->{'chemistry'} = '' unless (exists $smrtwellDetails->{'chemistry'});
$smrtwellDetails->{'customizedCondition'} = '' unless (exists $smrtwellDetails->{'customizedCondition'});
$smrtwellDetails->{'customizedMovieLength'} = '' unless (exists $smrtwellDetails->{'customizedMovieLength'});
$smrtwellDetails->{'comments'} = '' unless (exists $smrtwellDetails->{'comments'});

my $editNotice ='';
my $smrtrunIdList ='';
my $smrtwellWellList = '';
my $smrtwellCellNumber = '';
my $smrtwellMovieLength = '';
my $customizedMovieLength = '';
my $concentration = '';
my $polRatio = '';
my $smrtwellChemistry = '';
my $smrtwellCondition = '';
my $customizedCondition = '';

my $smrtrun=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$smrtrun->execute($smrtwell[6]);
my @smrtrun = $smrtrun->fetchrow_array();
if($smrtrun[7] > 0)
{
	$editNotice = "<div class='ui-state-error-text'>You are unable to change any settings of this $smrtrunStatus{$smrtrun[7]} Run Well!</div>";
	$smrtrunIdList = "$smrtrun[2]<input name='smrtrunId' id='editSmrtwellRunId' type='hidden' value='$smrtrun[0]' />";
	$smrtwellWellList = "$smrtwell[2] <input name='name' id='editSmrtwellName' type='hidden' value='$smrtwell[2]' />";
	$smrtwellCellNumber = "$smrtwell[4] <input name='cellNumber' id='editSmrtwellCellNumber' type='hidden' value='$smrtwell[4]' />";
	$smrtwellMovieLength = "$movieLength{$smrtwell[5]}<input name='movieLength' id='editSmrtwellMovieLength' type='hidden' value='$smrtwell[5]' />";
	$customizedMovieLength = "<span id='editSmrtwellCustomizedMovieLength'>$smrtwellDetails->{'customizedMovieLength'} <input name='customizedMovieLength' id='editSmrtwellCustomizedMovieLengthInput' type='hidden' value='$smrtwellDetails->{'customizedMovieLength'}' /></span>";
	$concentration = "$smrtwellDetails->{'concentration'} <input name='concentration' id='editSmrtwellConcentration' type='hidden' value='$smrtwellDetails->{'concentration'}' />";
	$polRatio = "$smrtwellDetails->{'polRatio'} <input name='polRatio' id='editSmrtwellPolRatio' type='hidden' value='$smrtwellDetails->{'polRatio'}' />";
	$smrtwellChemistry = "$chemistry{$smrtwellDetails->{'chemistry'}} <input name='chemistry' id='editSmrtwellChemistry' type='hidden' value='$smrtwellDetails->{'chemistry'}' />";
	$smrtwellCondition = "$condition{$smrtwell[7]}<input name='condition' id='editSmrtwellCondition' type='hidden' value='$smrtwell[7]' />";
	$customizedCondition = "<span id='editSmrtwellCustomizedCondition'>$smrtwellDetails->{'customizedCondition'} <input name='customizedCondition' id='editSmrtwellCustomizedConditionInput' type='hidden' value='$smrtwellDetails->{'customizedCondition'}' /></span>";
}
else
{
	my $allSmrtrun=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'smrtrun' AND barcode = 0");
	$allSmrtrun->execute();
	while (my @allSmrtrun = $allSmrtrun->fetchrow_array())
	{
		my $totalCellNumber = 0;
		my $checkCellNumber = $dbh->prepare("SELECT SUM(x) FROM matrix WHERE container LIKE 'smrtwell' AND z = ? AND id != ?");
		$checkCellNumber->execute($allSmrtrun[0],$smrtwellId);
		my @checkCellNumber=$checkCellNumber->fetchrow_array();
		$totalCellNumber = $checkCellNumber[0] if($checkCellNumber[0]);
		my $availableCell = $allSmrtrun[3] - $totalCellNumber;
		next if ($availableCell < 1);
		$smrtrunIdList .= ($smrtwell[6] eq $allSmrtrun[0]) ? "<option value='$allSmrtrun[0]' selected>$allSmrtrun[2] ($availableCell/$allSmrtrun[3] cells available)</option>" : "<option value='$allSmrtrun[0]'>$allSmrtrun[2] ($availableCell/$allSmrtrun[3] cells available)</option>";
	}
	$smrtrunIdList = "<select class='ui-widget-content ui-corner-all' name='smrtrunId' id='editSmrtwellRunId' onchange='hideShowCustomized();'>$smrtrunIdList<option value='0'>Create New Run</option></select>";

	$smrtwellWellList = "<table class='ui-widget-content ui-corner-all'><tr><td></td>";
	for(sort @{$wellx{96}})
	{
		$smrtwellWellList .= "<td style='text-align:center;' class='ui-state-highlight ui-corner-all'>$_</td>";
	}
	$smrtwellWellList .= "</tr>";

	for(sort @{$welly{96}})
	{
		my $y = $_;
		$smrtwellWellList .= "<tr><td style='text-align:center;' class='ui-state-highlight ui-corner-all'>$y</td>";
		for(sort @{$wellx{96}})
		{
			$smrtwellWellList .= ($smrtwell[2] eq "$y$_") ? "<td style='text-align:center;'><input type='radio' name='name' value='$y$_' title='$y$_' id='nameRadio$y$_' class='radioTootip' checked></td>" : "<td style='text-align:center;'><input type='radio' name='name' value='$y$_' title='$y$_' id='nameRadio$y$_' class='radioTootip'></td>";
		}
		$smrtwellWellList .= "</tr>";
	}
	$smrtwellWellList .= "</table>";

	$smrtwellCellNumber = "<input name='cellNumber' id='editSmrtwellCellNumber' size='2' type='text' maxlength='2' value='$smrtwell[4]' />";

	foreach (sort {$a <=> $b} keys %movieLength)
	{
		next unless ($_);
		$smrtwellMovieLength .= ($smrtwell[5] eq $_) ? "<option value='$_' selected>$movieLength{$_}</option>" : "<option value='$_'>$movieLength{$_}</option>";
	}
	$smrtwellMovieLength = ($smrtwell[5] eq 0) ? "<select class='ui-widget-content ui-corner-all' name='movieLength' id='editSmrtwellMovieLength' onchange='hideShowCustomized();'>$smrtwellMovieLength<option value='0' selected>$movieLength{0}</option></select>"
		: "<select class='ui-widget-content ui-corner-all' name='movieLength' id='editSmrtwellMovieLength' onchange='hideShowCustomized();'>$smrtwellMovieLength<option value='0'>$movieLength{0}</option></select>";
	$customizedMovieLength = "<input class='ui-widget-content ui-corner-all' name='customizedMovieLength' id='editSmrtwellCustomizedMovieLength' type='text' value='$smrtwellDetails->{'customizedMovieLength'}' placeholder='Customized Movie Length' size='15' type='text' maxlength='25'/>";
	$concentration = "<input class='ui-widget-content ui-corner-all' name='concentration' id='editSmrtwellConcentration' size='5' type='text' maxlength='10' value='$smrtwellDetails->{'concentration'}' />";
	$polRatio = "<input class='ui-widget-content ui-corner-all' name='polRatio' id='editSmrtwellPolRatio' size='5' type='text' maxlength='10' value='$smrtwellDetails->{'polRatio'}' />";
	foreach (sort keys %chemistry)
	{
		$smrtwellChemistry .= ($smrtwellDetails->{'chemistry'} eq $_) ? "<option value='$_' selected>$chemistry{$_}</option>" : "<option value='$_'>$chemistry{$_}</option>";
	}
	$smrtwellChemistry = "<select class='ui-widget-content ui-corner-all' name='chemistry' id='editSmrtwellChemistry'>$smrtwellChemistry</select>";

	foreach (sort {$a <=> $b} keys %condition)
	{
		next unless ($_);
		$smrtwellCondition .= ($smrtwell[7] eq $_) ? "<option value='$_' selected>$condition{$_}</option>" : "<option value='$_'>$condition{$_}</option>";
	}
	$smrtwellCondition = ($smrtwell[7] eq 0) ? "<select class='ui-widget-content ui-corner-all' name='condition' id='editSmrtwellCondition' onchange='hideShowCustomized();'>$smrtwellCondition<option value='0' selected>$condition{0}</option></select>"
		: "<select class='ui-widget-content ui-corner-all' name='condition' id='editSmrtwellCondition' onchange='hideShowCustomized();'>$smrtwellCondition<option value='0'>$condition{0}</option></select>";
	$customizedCondition = "<input class='ui-widget-content ui-corner-all' name='customizedCondition' id='editSmrtwellCustomizedCondition' type='text' value='$smrtwellDetails->{'customizedCondition'}' placeholder='Customized Condition' size='15' type='text' maxlength='25'/>";
}

my $paclib = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$paclib->execute($smrtwell[3]);
my @paclib=$paclib->fetchrow_array();
my $paclibToSample=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$paclibToSample->execute($paclib[6]);
my @paclibToSample = $paclibToSample->fetchrow_array();
my $sampleToService=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$sampleToService->execute($paclibToSample[6]);
my @sampleToService = $sampleToService->fetchrow_array();
my $serviceToProject=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$serviceToProject->execute($sampleToService[6]);
my @serviceToProject = $serviceToProject->fetchrow_array();

my $config = new config;
my $lastSmrtrun = $config->getFieldValueWithFieldName("smrtrun");
$lastSmrtrun = sprintf "%0*d", 3, $lastSmrtrun + 1;
my $datetimeString = strftime "%y %Y %m %d %H:%M:%S", localtime;
my ($tdYear,$year,$month,$dayofmonth,$time)=split(/\s+/,$datetimeString);
my $autoSmrtrunName = "$tdYear$month$dayofmonth\R$lastSmrtrun";
my $maxSmrtcell = 16;
my $smrtrunSmrtcell = '';
for (my $i = 1;$i <= $maxSmrtcell;$i++)
{
	$smrtrunSmrtcell .= ($i == 8) ? "<option value='$i' selected>$i</option>" : "<option value='$i'>$i</option>";
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$editNotice/$editNotice/g;
$html =~ s/\$smrtwellId/$smrtwellId/g;
$html =~ s/\$smrtrunId/$smrtrunIdList/g;
$html =~ s/\$autoSmrtrunName/$autoSmrtrunName/g;
$html =~ s/\$smrtrunSmrtcell/$smrtrunSmrtcell/g;
$html =~ s/\$smrtwellName/$smrtwell[2]/g;
$html =~ s/\$smrtwellWellList/$smrtwellWellList/g;
$html =~ s/\$smrtrunId/$smrtrunIdList/g;
$html =~ s/\$projectName/$serviceToProject[2]/g;
$html =~ s/\$serviceName/$sampleToService[2]/g;
$html =~ s/\$sampleName/$paclibToSample[2]/g;
$html =~ s/\$sampleId/$paclib[6]/g;
$html =~ s/\$paclibName/$paclib[2]/g;
$html =~ s/\$paclibId/$smrtwell[3]/g;
$html =~ s/\$smrtwellCellNumber/$smrtwellCellNumber/g;
$html =~ s/\$smrtwellMovieLength/$smrtwellMovieLength/g;
$html =~ s/\$customizedMovieLength/$customizedMovieLength/g;
$html =~ s/\$concentration/$concentration/g;
$html =~ s/\$polRatio/$polRatio/g;
$html =~ s/\$smrtwellChemistry/$smrtwellChemistry/g;
$html =~ s/\$smrtwellCondition/$smrtwellCondition/g;
$html =~ s/\$customizedCondition/$customizedCondition/g;
$html =~ s/\$loadingName/$smrtwellDetails->{'loadingName'}/g;
$html =~ s/\$comments/$smrtwellDetails->{'comments'}/g;
$html =~ s/\$smrtwellCreator/$smrtwell[9]/g;
$html =~ s/\$smrtwellEnteredDate/$smrtwell[10]/g;

print header;
print $html;

__DATA__
<form id="editSmrtwell" name="editSmrtwell" action="smrtwellSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="smrtwellId" id="editSmrtwellId" type="hidden" value="$smrtwellId" />
	$editNotice
	<table>
		<tr><td style='text-align:right;white-space: nowrap;'><label for="editSmrtwellRunId"><b>SMRT Run</b></label></td>
			<td>$smrtrunId
				<div id='editSmrtwellNewSmrtrun'>
					<table>
						<tr><td style='text-align:right;white-space: nowrap;'><label for="newSmrtrunName">Name</label></td>
							<td><input class='ui-widget-content ui-corner-all' name="smrtrunName" id="newSmrtrunName" value="$autoSmrtrunName" placeholder="Run Name" size="15" type="text" maxlength="32" /></td>
						</tr>
						<tr><td style='text-align:right;white-space: nowrap;'><label for="newSmrtrunSmrtcell"><b>SMRTCells</b></label></td>
							<td><select class='ui-widget-content ui-corner-all' name="smrtcell" id="newSmrtrunSmrtcell">$smrtrunSmrtcell</select></td>
						</tr>
						<tr><td style='text-align:right;white-space: nowrap;'><label for="newSmrtrunComments">Comments</label></td>
							<td><textarea class='ui-widget-content ui-corner-all' name="smrtrunComments" id="newSmrtrunComments" cols="40" rows="2" placeholder="Please write anything relevant to the run that may not be included on this form"></textarea></td>
						</tr>
					</table>
				</div>
			</td>
		</tr>
		<tr><td style='text-align:right;white-space: nowrap;'><label for="editSmrtwellName"><b>Sample Plate Well</b></label></td>
			<td>$smrtwellWellList<sup class='ui-state-disabled'>Last changed by $smrtwellCreator on $smrtwellEnteredDate</sup></td>
		</tr>
		<tr><td style='text-align:right;white-space: nowrap;'><label for="editSmrtwellPaclibId"><b>PacBio library</b></label></td>
			<td>$projectName > $serviceName > $sampleName > $paclibName<input name="paclibId" id="editSmrtwellPaclibId" type="hidden" value="$paclibId" /></td>
		</tr>
		<tr><td style='text-align:right;white-space: nowrap;'></td>
			<td><input class='ui-widget-content ui-corner-all' name="loadingName" id="editSmrtwellLoadingName" title="SMRT Sample Name" placeholder="Loading Name" size="20" type="text" value="$loadingName" readonly /></td>
		</tr>
		<tr><td style='text-align:right;white-space: nowrap;'><label for="editSmrtwellCellNumber"><b>Cell number</b></label></td><td>$smrtwellCellNumber</td></tr>
		<tr><td style='text-align:right;white-space: nowrap;'><label for="editSmrtwellMovieLength"><b>Movie length</b></label></td>
			<td>$smrtwellMovieLength $customizedMovieLength
			</td>
		</tr>
		<tr><td style='text-align:right;white-space: nowrap;'><label for="editSmrtwellConcentration"><b>DNA concentration</b></label></td><td>$concentration nM</td></tr>
		<tr><td style='text-align:right;white-space: nowrap;'><label for="editSmrtwellPolRatio"><b>Polymerase ratio</b></label></td><td>$polRatio</td></tr>
		<tr><td style='text-align:right;white-space: nowrap;'><label for="editSmrtwellChemistry"><b>Chemistry</b></label></td>
			<td>$smrtwellChemistry</td>
		</tr>
		<tr><td style='text-align:right;white-space: nowrap;'><label for="editSmrtwellCondition"><b>Condition</b></label></td>
			<td>$smrtwellCondition $customizedCondition
			</td>
		</tr>
		<tr><td style='text-align:right;white-space: nowrap;'><label for="editSmrtwellComments"><b>Comments</b></label></td>
			<td><textarea class='ui-widget-content ui-corner-all' name="comments" id="editSmrtwellComments" cols="40" rows="5" placeholder="Please write anything relevant to the well that may not be included on this form">$comments</textarea></td>
		</tr>
	</table>
</form>
<script>
hideShowCustomized();
$('#dialog').dialog("option", "title", "Edit Well $smrtwellName");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('editSmrtwell'); } }, { text: "Delete", click: function() { deleteItem($smrtwellId); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
$( "#editSmrtwellCellNumber" ).spinner({ min: 1, max: 96});
$( '.radioTootip' ).tooltip();
function hideShowCustomized()
{
	if($('#editSmrtwellRunId').val() == 0)
	{
		$('#editSmrtwellNewSmrtrun').show();
	}
	else
	{
		$('#editSmrtwellNewSmrtrun').hide();
	}
	if($('#editSmrtwellMovieLength').val() == 0)
	{
		$('#editSmrtwellCustomizedMovieLength').show();
	}
	else
	{
		$('#editSmrtwellCustomizedMovieLength').hide();
	}
	if($('#editSmrtwellCondition').val() == 0)
	{
		$('#editSmrtwellCustomizedCondition').show();
	}
	else
	{
		$('#editSmrtwellCustomizedCondition').hide();
	}
}
</script>
