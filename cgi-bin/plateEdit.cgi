#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
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

my $plateId = param ('plateId') || '';
my $plate = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$plate->execute($plateId);
my @plate=$plate->fetchrow_array();
my $noun = 'word';
my $plateToLibrary=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$plateToLibrary->execute($plate[6]);
my @plateToLibrary = $plateToLibrary->fetchrow_array();
my $libraryToProject=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$libraryToProject->execute($plateToLibrary[5]);
my @libraryToProject = $libraryToProject->fetchrow_array();

undef $/;# enable slurp mode
my $html = <DATA>;
if($plateToLibrary[5])
{
	my $sourceLibrary=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$sourceLibrary->execute($plateToLibrary[5]);
	my @sourceLibrary = $sourceLibrary->fetchrow_array();
	$html =~ s/\$rearrayingLibrary/<br><br>(Please provide source clone list from library '$sourceLibrary[2]', 1 clone\/line; no non-word characters allowed)/g;
	$noun = 'clone';
}
else
{
	$html =~ s/\$rearrayingLibrary//g;
}

$html =~ s/\$plateId/$plateId/g;
$html =~ s/\$project/$libraryToProject[2]/g;
$html =~ s/\$library/$plateToLibrary[2]/g;
$html =~ s/\$plateName/$plate[2]/g;
$html =~ s/\$plateCopyNumber/$plate[4]/g;
$html =~ s/\$plateCopiedFrom/$plate[5]/g;
$html =~ s/\$plateBarcode/$plate[7]/g;
$html =~ s/\$plateDescription/$plate[8]/g;
$html =~ s/\$plateCreator/$plate[9]/g;
$html =~ s/\$plateCreationDate/$plate[10]/g;
$html =~ s/\$noun/$noun/g;

print header;
print $html;

__DATA__
<form id="editPlate" name="editPlate" action="plateSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="plateId" id="editPlateId" type="hidden" value="$plateId" />
	<table>
	<tr><td style='text-align:right'><label for="editPlateName"><b>Plate number</b></label></td><td>$project > $library<input class='ui-widget-content ui-corner-all' name="name" id="editPlateName" size="10" type="text" maxlength="4" value="$plateName" readonly /><br><img alt='$plateBarcode' src='barcode.cgi?code=$plateBarcode'/><sup class='ui-state-disabled'>by $plateCreator on $plateCreationDate</sup></td></tr>
	<tr><td style='text-align:right'><label for="editPlateCopyNumber"><b>Copy</b></label></td><td><input class='ui-widget-content ui-corner-all' name="copyNumber" id="editPlateCopyNumber" size="4" type="text" maxlength="2" VALUE="$plateCopyNumber" readonly /> <label for="editPlateCopiedFrom"><b>from</b></label> <input class='ui-widget-content ui-corner-all' name="copiedFrom" id="editPlateCopiedFrom" size="4" type="text" maxlength="2" VALUE="$plateCopiedFrom" readonly /></td></tr>
	<tr><td style='text-align:right'><label for="editPlateDescription"><b>Description</b></label><br><sub id='editPlateDescription_count' style='display:none'></sub>$rearrayingLibrary</td><td><textarea class='ui-widget-content ui-corner-all word_count' name="description" id="editPlateDescription" cols="50" rows="10">$plateDescription</textarea></td></tr>
	</table>
</form>
<script>
$('#dialog').dialog("option", "title", "Edit Plate");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('editPlate'); } }, { text: "Delete", click: function() { deleteItem($plateId); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
wordCount('$noun');
</script>