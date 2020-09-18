#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
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

my $paclibId = param ('paclibId') || '';

my %wellx = ( 96 => ['01' .. '12'], 384 =>['01' .. '24']);
my %welly = ( 96 => ['A' .. 'H'], 384 =>['A' .. 'P']);

my $smrtwellWellList = "<table class='ui-widget-content ui-corner-all'><tr><td></td>";
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
		$smrtwellWellList .= "<td style='text-align:center;'><input type='radio' name='name' value='$y$_' title='$y$_' id='nameRadio$y$_' class='radioTootip'></td>";
	}
	$smrtwellWellList .= "</tr>";
}
$smrtwellWellList .= "</table>";

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
my $defaultMovieLength = 240;
my $smrtwellMovieLength = '';
foreach (sort {$a <=> $b} keys %movieLength)
{
	next unless ($_);
	$smrtwellMovieLength .= ($defaultMovieLength eq $_) ?  "<option value='$_' selected>$movieLength{$_}</option>" : "<option value='$_'>$movieLength{$_}</option>";
}
$smrtwellMovieLength .= "<option value='0'>$movieLength{0}</option>";

my %chemistry = (
	P4=>'P4',
	P5=>'P5',
	P6=>'P6'
	);
my $defaultChemistry = 'P5';
my $smrtwellChemistry = '';
foreach (sort keys %chemistry)
{
	$smrtwellChemistry .= ($defaultChemistry eq $_) ?  "<option value='$_' selected>$chemistry{$_}</option>" : "<option value='$_'>$chemistry{$_}</option>";
}
my %condition = (
	0=>'Customized:',
	1=>'100%',
	2=>'75%'
	);
my $smrtwellCondition = '';
foreach (sort {$a <=> $b} keys %condition)
{
	next unless ($_);
	$smrtwellCondition .= "<option value='$_'>$condition{$_}</option>";
}
$smrtwellCondition .= "<option value='0'>$condition{0}</option>";

my $paclib = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$paclib->execute($paclibId);
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

my $smrtrunIdList ='';
my $smrtrun=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'smrtrun' AND barcode = 0");
$smrtrun->execute();
while (my @smrtrun = $smrtrun->fetchrow_array())
{
	my $totalCellNumber = 0;
	my $checkCellNumber = $dbh->prepare("SELECT SUM(x) FROM matrix WHERE container LIKE 'smrtwell' AND z = ?");
	$checkCellNumber->execute($smrtrun[0]);
	my @checkCellNumber=$checkCellNumber->fetchrow_array();
	$totalCellNumber = $checkCellNumber[0] if($checkCellNumber[0]);
	my $availableCell = $smrtrun[3] - $totalCellNumber;
	next if ($availableCell < 1);
	$smrtrunIdList .= "<option value='$smrtrun[0]'>$smrtrun[2] ($availableCell/$smrtrun[3] cells available)</option>";
}
$smrtrunIdList .= "<option value='0'>Create New Run</option>";

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
$html =~ s/\$smrtrunId/$smrtrunIdList/g;
$html =~ s/\$autoSmrtrunName/$autoSmrtrunName/g;
$html =~ s/\$smrtrunSmrtcell/$smrtrunSmrtcell/g;
$html =~ s/\$projectName/$serviceToProject[2]/g;
$html =~ s/\$serviceName/$sampleToService[2]/g;
$html =~ s/\$sampleName/$paclibToSample[2]/g;
$html =~ s/\$sampleId/$paclib[6]/g;
$html =~ s/\$paclibName/$paclib[2]/g;
$html =~ s/\$paclibId/$paclibId/g;
$html =~ s/\$smrtwellWellList/$smrtwellWellList/g;
$html =~ s/\$smrtwellMovieLength/$smrtwellMovieLength/g;
$html =~ s/\$smrtwellChemistry/$smrtwellChemistry/g;
$html =~ s/\$smrtwellCondition/$smrtwellCondition/g;

print header;
print $html;

__DATA__
<form id="newSmrtwell" name="newSmrtwell" action="smrtwellSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<table>
		<tr><td style='text-align:right;white-space: nowrap;'><label for="newSmrtwellRunId"><b>SMRT Run</b></label></td>
			<td><select class='ui-widget-content ui-corner-all' name="smrtrunId" id="newSmrtwellRunId" onchange="hideShowCustomized();">$smrtrunId</select>
				<div id='newSmrtwellNewSmrtrun'>
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
		<tr><td style='text-align:right;white-space: nowrap;'><label for="newSmrtwellName"><b>Sample Plate Well</b></label></td>
			<td>$smrtwellWellList</td>
		</tr>
		<tr><td style='text-align:right;white-space: nowrap;'><label for="newSmrtwellPaclibId"><b>PacBio library</b></label></td>
			<td>$projectName > $serviceName > $sampleName > $paclibName<input name="paclibId" id="newSmrtwellPaclibId" type="hidden" value="$paclibId" /></td>
		</tr>
		<tr><td style='text-align:right;white-space: nowrap;'><label for="newSmrtwellCellNumber"><b>Cell number</b></label></td><td><input name="cellNumber" id="newSmrtwellCellNumber" size="2" type="text" maxlength="2" value="1" /></td></tr>
		<tr><td style='text-align:right;white-space: nowrap;'><label for="newSmrtwellMovieLength"><b>Movie length</b></label></td>
			<td><select class='ui-widget-content ui-corner-all' name="movieLength" id="newSmrtwellMovieLength" onchange="hideShowCustomized();">$smrtwellMovieLength</select>
			<input class='ui-widget-content ui-corner-all' name="customizedMovieLength" id="newSmrtwellCustomizedMovieLength" type="text" value="" placeholder="Customized Movie Length" size="15" type="text" maxlength="25" />
			</td>
		</tr>
		<tr><td style='text-align:right;white-space: nowrap;'><label for="newSmrtwellConcentration"><b>DNA concentration</b></label></td><td><input class='ui-widget-content ui-corner-all' name="concentration" id="newSmrtwellConcentration" size="5" type="text" maxlength="10" value="" />nM</td></tr>
		<tr><td style='text-align:right;white-space: nowrap;'><label for="newSmrtwellPolRatio"><b>Polymerase ratio</b></label></td><td><input class='ui-widget-content ui-corner-all' name="polRatio" id="newSmrtwellPolRatio" size="5" type="text" maxlength="10" value="" /></td></tr>
		<tr><td style='text-align:right;white-space: nowrap;'><label for="newSmrtwellChemistry"><b>Chemistry</b></label></td>
			<td><select class='ui-widget-content ui-corner-all' name="chemistry" id="newSmrtwellChemistry">$smrtwellChemistry</select></td>
		</tr>
		<tr><td style='text-align:right;white-space: nowrap;'><label for="newSmrtwellCondition"><b>Condition</b></label></td>
			<td><select class='ui-widget-content ui-corner-all' name="condition" id="newSmrtwellCondition" onchange="hideShowCustomized();">$smrtwellCondition</select>
			<input class='ui-widget-content ui-corner-all' name="customizedCondition" id="newSmrtwellCustomizedCondition" type="text" value="" placeholder="Customized Condition" size="15" type="text" maxlength="25" />
			</td>
		</tr>
		<tr><td style='text-align:right;white-space: nowrap;'><label for="newSmrtwellComments"><b>Comments</b></label></td>
			<td><textarea class='ui-widget-content ui-corner-all' name="comments" id="newSmrtwellComments" cols="40" rows="5" placeholder="Please write anything relevant to the well that may not be included on this form"></textarea></td>
		</tr>
	</table>
</form>
<script>
hideShowCustomized();
$('#dialog').dialog("option", "title", "New Well");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('newSmrtwell'); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
$( "#newSmrtwellCellNumber" ).spinner({ min: 1, max: 96});
$( '.radioTootip' ).tooltip();
function hideShowCustomized()
{
	if($('#newSmrtwellRunId').val() == 0)
	{
		$('#newSmrtwellNewSmrtrun').show();
	}
	else
	{
		$('#newSmrtwellNewSmrtrun').hide();
	}
	if($('#newSmrtwellMovieLength').val() == 0)
	{
		$('#newSmrtwellCustomizedMovieLength').show();
	}
	else
	{
		$('#newSmrtwellCustomizedMovieLength').hide();
	}
	if($('#newSmrtwellCondition').val() == 0)
	{
		$('#newSmrtwellCustomizedCondition').show();
	}
	else
	{
		$('#newSmrtwellCustomizedCondition').hide();
	}
}
</script>
