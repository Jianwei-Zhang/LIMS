#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use DBI;
use lib "lib/";
use lib "lib/pangu";
use pangu;
use userConfig;
use userCookie;

my $userCookie = new userCookie;
my $userId = (cookie('cid')) ? $userCookie->checkCookie(cookie('cid')) : 0;
exit if (!$userId);

my $commoncfg = readConfig("main.conf");
my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});
my $userConfig = new userConfig;
my $SEQTOGNMMINOVERLAP = $userConfig->getFieldValueWithUserIdAndFieldName($userId,"SEQTOGNMMINOVERLAP");
my $ENDTOENDMINOVERLAP = $userConfig->getFieldValueWithUserIdAndFieldName($userId,"ENDTOENDMINOVERLAP");
my $ENDTOENDIDENTITY = $userConfig->getFieldValueWithUserIdAndFieldName($userId,"ENDTOENDIDENTITY");

my $objectComponent;
$objectComponent->{0} = "Unknown";
$objectComponent->{1} = "Chr-Seq";
$objectComponent->{2} = "Ctg-Seq";

print header;

my $assemblyId = param ('assemblyId') || '';
my $assembly=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$assembly->execute($assemblyId);
my @assembly = $assembly->fetchrow_array();

unless ($assembly[7] == 1 || $assembly[7] == 0) # exit if for frozen or running assembly
{
	print <<END;
<script>
	closeDialog();
	errorPop("This assembly is running or frozen.");
</script>	
END
	exit;
}

my $checkRunningAssembly=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assembly' AND barcode < 0 ");
$checkRunningAssembly->execute();
if ($checkRunningAssembly->rows > 0)
{
	print <<END;
<script>
	closeDialog();
	errorPop("There is another running assembly. Please try again later!");
</script>	
END
	exit;
}

my $target=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$target->execute($assembly[4]);
my @target = $target->fetchrow_array();
my $fpcOrAgpId = '';
if($target[1] eq 'library')
{
	my $fpcList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'fpc' AND barcode = ?");
	$fpcList->execute($assembly[4]);
	while (my @fpcList = $fpcList->fetchrow_array())
	{
		$fpcOrAgpId .= ($fpcList[0] eq $assembly[6] ) ?
			"<option value='$fpcList[0]' title='$fpcList[8]' selected>FPC: $fpcList[2] v.$fpcList[3]</option>" :
			"<option value='$fpcList[0]' title='$fpcList[8]'>FPC: $fpcList[2] v.$fpcList[3]</option>";
	}
}
elsif($target[1] eq 'genome')
{
	my $agpList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'agp' AND x = ?");
	$agpList->execute($assembly[4]);
	while (my @agpList = $agpList->fetchrow_array())
	{
		$fpcOrAgpId .= ($agpList[0] eq $assembly[6] ) ?
			"<option value='$agpList[0]' title='$objectComponent->{$agpList[5]} AGP' selected>AGP: $agpList[2] v.$agpList[3]</option>" :
			"<option value='$agpList[0]' title='$objectComponent->{$agpList[5]} AGP'>AGP: $agpList[2] v.$agpList[3]</option>";
	}
}
if($fpcOrAgpId)
{
	$fpcOrAgpId = "<option class='ui-state-error-text' value='0'>Please select a reference</option>".$fpcOrAgpId;
}
else
{
	$fpcOrAgpId = "<option class='ui-state-error-text' value='0'>No reference available</option>".$fpcOrAgpId;
}
my $refGenomeId = '';
my $genomeList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'genome'");
$genomeList->execute();
while (my @genomeList = $genomeList->fetchrow_array())
{
	next if ($genomeList[0] eq $assembly[4]);
	next if ($genomeList[5] < 1); #remove not for reference
	$refGenomeId .= ($genomeList[0] eq $assembly[5] ) ?
		"<option value='$genomeList[0]' title='$genomeList[8]' selected>$genomeList[2]</option>" :
		"<option value='$genomeList[0]' title='$genomeList[8]'>$genomeList[2]</option>";
}
if($refGenomeId)
{
	$refGenomeId = "<option class='ui-state-error-text' value='0'>Please select a genome</option>".$refGenomeId;
}
else
{
	$refGenomeId = "<option class='ui-state-error-text' value='0'>No reference genome available</option>".$refGenomeId;
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$assemblyId/$assemblyId/g;
$html =~ s/\$assemblyName/$assembly[2]/g;
$html =~ s/\$fpcOrAgpId/$fpcOrAgpId/g;
$html =~ s/\$refGenomeId/$refGenomeId/g;
$html =~ s/\$SEQTOGNMMINOVERLAP/$SEQTOGNMMINOVERLAP/g;
$html =~ s/\$ENDTOENDMINOVERLAP/$ENDTOENDMINOVERLAP/g;
$html =~ s/\$ENDTOENDIDENTITY/$ENDTOENDIDENTITY/g;

print $html;

__DATA__
<form id="newAssembly" name="newAssembly" action="assemblyRun.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="assemblyId" id="newAssemblyId" type="hidden" value="$assemblyId" />
	<table>
	<tr><td style='text-align:left'><label for="newAssemblyReplace"><b>Start Over a New Assembly?</b></label></td>
		<td>
		<div id="newAssemblyReplace">
			<input type="radio" id="newAssemblyReplaceRadio2" name="replace" value="1"><label for="newAssemblyReplaceRadio2">Yes</label>
			<input type="radio" id="newAssemblyReplaceRadio1" name="replace" value="0" checked="checked"><label for="newAssemblyReplaceRadio1">No</label>
		</div>
	</td></tr>
	<tr><td style='text-align:left'><label for="newAssemblyWithAllSequences"><b>With All Sequences available?</b></label></td>
		<td>
		<div id="newAssemblyWithAllSequences">
			<input type="radio" id="newAssemblyWithAllSequencesRadio2" name="withAllSequences" value="1" checked="checked"><label for="newAssemblyWithAllSequencesRadio2">Yes</label>
			<input type="radio" id="newAssemblyWithAllSequencesRadio1" name="withAllSequences" value="0"><label for="newAssemblyWithAllSequencesRadio1">No</label>
		</div>
	</td></tr>
	<tr><td style='text-align:left'><label for="newAssemblySeqMinLength"><b>AssemblySeq Minimum Length</b></label></td>
		<td><input name="assemblySeqMinLength" id="newAssemblySeqMinLength" size="4" type="text" maxlength="6" VALUE="0" />(bp)</td>
	</tr>
	</table>
	<table>
	<tr><td style='text-align:right'><label for="newAssemblyFpcOrAgp"><b>Physical Reference</b></label></td>
	<td><select class='ui-widget-content ui-corner-all' name="fpcOrAgpId" id="newAssemblyFpcOrAgp">$fpcOrAgpId</select></td>
	</tr>
	</table>
	<table>
	<tr><td style='text-align:right'><label for="newAssemblyGenome"><b>Reference Genome</b></label></td>
	<td><select class='ui-widget-content ui-corner-all' name="refGenomeId" id="newAssemblyGenome">$refGenomeId</select><br>
			<label for="newAssemblyAlignmentBlockSize">Alignment Block Size</label> <input name="alignmentBlockSize" id="newAssemblyAlignmentBlockSize" size="4" type="text" maxlength="6" VALUE="$SEQTOGNMMINOVERLAP" />(bp)<br>
			<label for="newAssemblyAlignmentCoverage">Alignment Coverage Percent</label> <input name="alignmentCoverage" id="newAssemblyAlignmentCoverage" size="2" type="text" maxlength="3" VALUE="25" />(%)<br>
			<input type="checkbox" id="newAssemblyAssignChrCheckbox" name="assignChr" value="1" checked="checked"> <label for="newAssemblyAssignChrCheckbox">Assign chromosome number for unplaced contigs</label><br>
			<input type="checkbox" id="newAssemblyRepositionCheckbox" name="reposition" value="1"> <label for="newAssemblyRepositionCheckbox">Re-estimate position for anchored contigs</label><br>
			<input type="checkbox" id="newAssemblyOrientContigsCheckbox" name="orientContigs" value="1" checked="checked"> <label for="newAssemblyOrientContigsCheckbox">Orient contigs based-on reference genome</label>
	</td>
	</tr>
	</table>
	<hr>
	<table>
	<tr><td style='text-align:left' class='ui-state-highlight ui-corner-all' rowspan="2"><label for="newAssemblyEndToEnd"><b>End-to-End Merge</b></label>
		<div id="newAssemblyEndToEnd">
			<input type="radio" id="newAssemblyEndToEndRadio2" name="endToEnd" value="1"><label for="newAssemblyEndToEndRadio2">Yes</label>
			<input type="radio" id="newAssemblyEndToEndRadio1" name="endToEnd" value="0" checked="checked"><label for="newAssemblyEndToEndRadio1">No</label>
		</div>
	</td>
	<td style='text-align:right'><label for="newAssemblyEndToEndMinOverlap"><b>Minimum Overlap</b></label></td><td><input name="minOverlapEndToEnd" id="newAssemblyEndToEndMinOverlap" size="4" type="text" maxlength="6" VALUE="$ENDTOENDMINOVERLAP" />(bp)</td>
	</tr>
	<tr><td style='text-align:right'><label for="newAssemblyEndToEndIdentity"><b>Minimum Identity</b></label></td><td><input name="identityEndToEnd" id="newAssemblyEndToEndIdentity" size="4" type="text" maxlength="4" VALUE="$ENDTOENDIDENTITY" />(%)</td></tr>
	</table>
	<hr>
	<table>
	<tr><td><b>Post-Assembly Options</b></td><td style='text-align:right'><b>Filter</b></td>
		<td>
		<div id="newAssemblyRedundancyFilterSeq">
			<input type="checkbox" id="newAssemblyRedundancyFilterSeqCheckbox" name="redundancyFilterSeq" value="1" checked="checked"><label for="newAssemblyRedundancyFilterSeqCheckbox">Redundant Sequences</label>
		</div>
	</td>
	</tr>
	<tr><td></td><td style='text-align:right'></td>
		<td>
		<div id="newAssemblyRedundancyFilterOverlap">
			<input type="checkbox" id="newAssemblyRedundancyFilterOverlapCheckbox" name="redundancyFilterOverlap" value="1" checked="checked"><label for="newAssemblyRedundancyFilterOverlapCheckbox">Redundant Overlaps</label>
		</div>
	</td>
	</tr>
	<tr><td></td><td style='text-align:right'><label for="newAssemblyOrientSeqs"><b>Auto-Orient Sequences</b></label></td>
		<td>
		<div id="newAssemblyOrientSeqs">
			<input type="radio" id="newAssemblyOrientSeqsRadio2" name="orientSeqs" value="1" checked="checked"><label for="newAssemblyOrientSeqsRadio2">Yes</label>
			<input type="radio" id="newAssemblyOrientSeqsRadio1" name="orientSeqs" value="0"><label for="newAssemblyOrientSeqsRadio1">No</label>
		</div>
	</td>
	</tr>
	<tr><td></td><td style='text-align:right'><label for="newAssemblyRenumber"><b>Renumber Contigs</b></label></td>
		<td>
		<div id="newAssemblyRenumber">
			<input type="radio" id="newAssemblyRenumberRadio2" name="renumber" value="1" checked="checked"><label for="newAssemblyRenumberRadio2">Yes</label>
			<input type="radio" id="newAssemblyRenumberRadio1" name="renumber" value="0"><label for="newAssemblyRenumberRadio1">No</label>
		</div>
	</td>
	</tr>
	<tr><td></td><td style='text-align:left' colspan='2'><input type="checkbox" id="newEmailNotificationCheckbox" name="emailNotification" value="1" checked="checked"><label for="newEmailNotificationCheckbox">Email a notification after assembly is done.</label></td></tr>
	</table>
</form>
<script>
$('#dialog').dialog("option", "title", "Run Assembly ($assemblyName)");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Run Assembly", click: function() { submitForm('newAssembly'); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
$( "#newAssemblyReplace" ).buttonset();
$( "#newAssemblyWithAllSequences" ).buttonset();
$( "#newAssemblySeqMinLength" ).spinner({
	min: 0,
	max: 10000,
	step: 100,
	start: function( event, ui ) {
		var current = $( "#newAssemblySeqMinLength" ).spinner( "value" );
		if(current >= 2000)
		{
			$( "#newAssemblySeqMinLength" ).spinner( "option", "step", 500 );
		}
		else
		{
			$( "#newAssemblySeqMinLength" ).spinner( "option", "step", 100 );
		}
	}
});
$( "#newAssemblyAlignmentBlockSize" ).spinner({ min: 0, max: 30000, step: 1000});
$( "#newAssemblyAlignmentCoverage" ).spinner({ min: 10, max: 100, step: 5});
$( "#newAssemblyEndToEnd" ).buttonset();
$( "#newAssemblyEndToEndMinOverlap" ).spinner({
	min: 500,
	max: 99999,
	step: 100,
	start: function( event, ui ) {
		var current = $( "#newAssemblyEndToEndMinOverlap" ).spinner( "value" );
		if(current >= 2000)
		{
			$( "#newAssemblyEndToEndMinOverlap" ).spinner( "option", "step", 500 );
		}
		else
		{
			$( "#newAssemblyEndToEndMinOverlap" ).spinner( "option", "step", 100 );
		}
	}
});
$( "#newAssemblyEndToEndIdentity" ).spinner({ min: 90, max: 100});
$( "#newAssemblyOrientSeqs" ).buttonset();
$( "#newAssemblyRenumber" ).buttonset();
</script>