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

my $genomeId = param ('genomeId') || '';
my $genome = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$genome->execute($genomeId);
my @genome=$genome->fetchrow_array();

my $genomeStatus;
$genomeStatus->{0} = "Sequences not loaded";
$genomeStatus->{-1} = "Sequences are being loaded";
$genomeStatus->{1} = "$genome[3] sequences loaded";

my $objectComponent;
$objectComponent->{0} = "Unknown";
$objectComponent->{1} = "Chr-Seq";
$objectComponent->{2} = "Ctg-Seq";

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$genomeId/$genomeId/g;
$html =~ s/\$genomeName/$genome[2]/g;

if($genome[4])
{
	$html =~ s/\$editGenomeForAssembly/ checked="checked"/g;
}
else
{
	$html =~ s/\$editGenomeForAssembly//g;
}
if($genome[5])
{
	$html =~ s/\$editGenomeAsReference/ checked="checked"/g;
}
else
{
	$html =~ s/\$editGenomeAsReference//g;
}

my $libraryId = "<option class='ui-state-error-text' value='0'>None</option>";
my $libraryList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'library'");
$libraryList->execute();
if ($libraryList->rows > 0)
{
	my $libraryListResult;
	while (my @libraryList = $libraryList->fetchrow_array())
	{
		@{$libraryListResult->{$libraryList[2]}} = @libraryList;
	}
	foreach (sort {uc ($a) cmp uc($b)} keys %$libraryListResult)
	{
		my @libraryList = @{$libraryListResult->{$_}};
		my $libraryDetails = decode_json $libraryList[8];
		$libraryDetails->{'comments'} = escapeHTML($libraryDetails->{'comments'});
		$libraryId .= ($libraryList[0] eq $genome[6] ) ?
			"<option value='$libraryList[0]' title='$libraryDetails->{'comments'}' selected>Library: $libraryList[2]</option>" :
			"<option value='$libraryList[0]' title='$libraryDetails->{'comments'}'>Library: $libraryList[2]</option>";
	}
}

my $agpAvailable = "";
my $agpList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'agp' AND x = ? ORDER BY o");
$agpList->execute($genome[0]);
while (my @agpList = $agpList->fetchrow_array())
{
	$agpAvailable .= ($agpAvailable) ? "<li><a href='download.cgi?agpId=$agpList[0]' target='hiddenFrame' title='$objectComponent->{$agpList[5]} AGP v.$agpList[3]'>$agpList[2]</a><span class='ui-icon ui-icon-trash' onclick='deleteItem($agpList[0])' title='Delete this AGP' style='display:inline-block;'></span></li>" : "<br>Available AGP:<ul><li><a href='download.cgi?agpId=$agpList[0]' target='hiddenFrame' title='$objectComponent->{$agpList[5]} AGP v.$agpList[3]'>$agpList[2]</a><span class='ui-icon ui-icon-trash'  onclick='deleteItem($agpList[0])' title='Delete this AGP' style='display:inline-block;'></span></li>";
}
$agpAvailable .= ($agpAvailable) ? '</ul>' : '';
my $agpObjectComponent;
for (sort keys %{$objectComponent})
{
	$agpObjectComponent .= ($_ == 1) ? "<option value='$_' title='$objectComponent->{$_}' selected>$objectComponent->{$_}</option>" : "<option value='$_' title='$objectComponent->{$_}'>$objectComponent->{$_}</option>";
}

$html =~ s/\$libraryId/$libraryId/g;
$html =~ s/\$agpObjectComponent/$agpObjectComponent/g;
$html =~ s/\$agpAvailable/$agpAvailable/g;
$html =~ s/\$genomeStatus/$genomeStatus->{$genome[7]}/g;
$html =~ s/\$genomeDescription/$genome[8]/g;
$html =~ s/\$genomeCreator/$genome[9]/g;
$html =~ s/\$genomeCreationDate/$genome[10]/g;
print header;
print $html;

__DATA__
<form id="editGenome" name="editGenome" action="genomeSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="genomeId" id="editGenomeId" type="hidden" value="$genomeId" />
	<table>
	<tr><td style='text-align:right'><label for="editGenomeName"><b>Genome Name</b></label></td><td><input class='ui-widget-content ui-corner-all' name="name" id="editGenomeName" size="40" type="text" maxlength="32" value="$genomeName"/><br><sup class='ui-state-disabled'>$genomeStatus by $genomeCreator on $genomeCreationDate</sup></td></tr>
	<tr><td style='text-align:right' rowspan='3'><label for="editGenomeFile"><b>Sequence File</b></label></td><td><sup class="ui-state-error-text">Note: Sequence Id should contain no more than 32 characters.</sup><br><input name="genomeFile" id="editGenomeFile" type="file" />(in FASTA format)</td></tr>
	<tr><td>or <input name="genomeFilePath" id="editGenomeFilePath" type="text" />(On-server file name with full path)</td></tr>
	<tr>
		<td><div class="ui-state-error-text">The New Sequence File will be applied to</div>
			<div id="editGenomeReplace">
			<input type="radio" id="editGenomeReplaceRadio1" name="replace" value="1"><label for="editGenomeReplaceRadio1">replace the existing data set</label>
			<input type="radio" id="editGenomeReplaceRadio2" name="replace" value="0" checked="checked"><label for="editGenomeReplaceRadio2">append to the existing data set</label>
			</div>
		</td>
	</tr>
	<tr><td style='text-align:right'><b>Import Options</b></td>
		<td>
			<hr>
			<input type="checkbox" id="editGenomeAssignChr" name="assignChr" value="1"><label for="editGenomeAssignChr">Assign chromosome number based on sequence name</label><br><sub>Manual assignment for sequences on unknown chromosomes is required after uploading.</sub><hr width="80%">
			<input type="checkbox" id="editGenomeSplit" name="split" value="1"> <label for="editGenomeSplit">Split gapped sequences</label><hr width="80%">
			<input type="checkbox" id="editGenomeForAssembly" name="forAssembly" value="1"$editGenomeForAssembly><label for="editGenomeForAssembly">Enable this genome to be reassembled with GPM</label><br>
			<sub><label for="editAgpFile"><b>New <select class='ui-widget-content ui-corner-all' name='agpObjectComponent' id='editAgpObjectComponent'>$agpObjectComponent</select> AGP file for guiding re-assembling</b></label></sub>  <input name="agpFile" id="editAgpFile" type="file" />(<a title="You may upload an AGP with a long file name, but only the first 32 characters will be saved.">Maximum 32 characters</a>) $agpAvailable<hr width="80%">
			<input type="checkbox" id="editGenomeAsReference" name="asReference" value="1"$editGenomeAsReference><label for="editGenomeAsReference">Enable this genome to be used as a reference</label>
			<hr>
		</td>
	</tr>
	<tr><td style='text-align:right'><label for="editGenomeLibraryId"><b>Link to</b></label></td><td><select class='ui-widget-content ui-corner-all' name='libraryId' id='editGenomeLibraryId'>$libraryId</select></td></tr>
	<tr><td style='text-align:right'><label for="editGenomeDescription"><b>Description</b></label></td><td><textarea class='ui-widget-content ui-corner-all' name="description" id="editGenomeDescription" cols="50" rows="10">$genomeDescription</textarea></td></tr>
	</table>
</form>
<script>
$( "#editGenomeReplace" ).buttonset();
$('#dialog').dialog("option", "title", "Edit Genome");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('editGenome'); } }, { text: "Delete", click: function() { deleteItem($genomeId); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
</script>