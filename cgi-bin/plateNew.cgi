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

my $libraryId = param ('libraryId') || '';
my $noun = 'word';
my $plateId;
my $plateInLibrary=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'plate' AND z = $libraryId ORDER BY o");
$plateInLibrary->execute();
while (my @plateInLibrary = $plateInLibrary->fetchrow_array())
{
	$plateId->{$plateInLibrary[2]}->{$plateInLibrary[4]} = $plateInLibrary[0];
}
my $lastPlate;
for (sort keys %{$plateId})
{
	$lastPlate = $_;
}
my $startPlate = sprintf "%0*d", 4, $lastPlate + 1;
my $maxNofPlates = 9999 - $lastPlate;

my $library=$dbh->prepare("SELECT * FROM matrix WHERE id = $libraryId");
$library->execute();
my @library = $library->fetchrow_array();
my $libraryToProject=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$libraryToProject->execute($library[6]);
my @libraryToProject = $libraryToProject->fetchrow_array();
$maxNofPlates = 1 if ($library[5]);

undef $/;# enable slurp mode
my $html = <DATA>;
if($library[5])
{
	my $sourceLibrary=$dbh->prepare("SELECT * FROM matrix WHERE id = $library[5]");
	$sourceLibrary->execute();
	my @sourceLibrary = $sourceLibrary->fetchrow_array();
	$html =~ s/\$rearrayingLibrary/<br><br>(Please provide source clone list from library '$sourceLibrary[2]', 1 clone\/line; no non-word characters allowed)/g;
	$noun = 'clone';
}
else
{
	$html =~ s/\$rearrayingLibrary//g;
}

$html =~ s/\$libraryId/$libraryId/g;
$html =~ s/\$project/$libraryToProject[2]/g;
$html =~ s/\$libraryName/$library[2]/g;
$html =~ s/\$startPlate/$startPlate/g;
$html =~ s/\$maxNofPlates/$maxNofPlates/g;
$html =~ s/\$noun/$noun/g;

print header;
print $html;

__DATA__
<form id="newPlate" name="newPlate" action="plateSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="libraryId" id="libraryId" type="hidden" value="$libraryId" />
	<table>
	<tr><td style='text-align:right'><label for="newPlateName"><b>Plate start at</b></label></td><td>$project > $libraryName<input class='ui-widget-content ui-corner-all' name="name" id="newPlateName" size="10" type="text" maxlength="4" value="$startPlate" readonly /> (Copy 0)</td></tr>
	<tr><td style='text-align:right'><label for="newPlateNofPlates"><b>Number of plates</b></label></td><td><input name="nofPlates" id="newPlateNofPlates" size="4" type="text" maxlength="4" VALUE="1" /></td></tr>
	<tr><td style='text-align:right'><label for="newPlateDescription"><b>Description</b></label><br><sub id='newPlateDescription_count' style='display:none'></sub>$rearrayingLibrary</td><td><textarea class='ui-widget-content ui-corner-all word_count' name="description" id="newPlateDescription" cols="50" rows="10"></textarea></td></tr>
	</table>
</form>
<script>
$( "#newPlateNofPlates" ).spinner({ min: 1, max: $maxNofPlates});
$('#dialog').dialog("option", "title", "New Plate");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('newPlate'); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
wordCount('$noun');
</script>