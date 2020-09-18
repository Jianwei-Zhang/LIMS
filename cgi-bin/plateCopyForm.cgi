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

my $library = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$library->execute($libraryId);
my @library=$library->fetchrow_array();

my $allCopyNumber;
my $copyNumberInLibrary=$dbh->prepare("SELECT x FROM matrix WHERE container LIKE 'plate' AND z = $libraryId GROUP BY x");
$copyNumberInLibrary->execute();
while (my @copyNumberInLibrary = $copyNumberInLibrary->fetchrow_array())
{
	$allCopyNumber->{$copyNumberInLibrary[0]} = 1;
}
my $newCopy;
my $plateOriginalCopy;
for (sort {$a <=> $b} keys %{$allCopyNumber})
{
	$plateOriginalCopy .= "<option value='$_'>Copy $_</option>";
	$newCopy = $_ + 1;
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$libraryId/$libraryId/g;
$html =~ s/\$newCopy/$newCopy/g;
$html =~ s/\$libraryName/$library[2]/g;
$html =~ s/\$plateOriginalCopy/$plateOriginalCopy/g;

print header;
print $html;

__DATA__
<form id="copyPlate" name="copyPlate" action="plateCopy.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<h3 class='ui-state-error-text'>Are you sure to copy the below plates?</h3>
	<input name="libraryId" id="libraryId" type="hidden" value="$libraryId" />
	<table>
	<tr><td style='text-align:right'><label for="copyPlateOriginalCopy"><b>Orignal Plates</b></label></td><td><select class='ui-widget-content ui-corner-all' name="originalCopy" id="copyPlateOriginalCopy">$plateOriginalCopy</select> of '$libraryName'</td></tr>
	<tr><td style='text-align:right'><label for="copyPlateNewCopy"><b>New Copy Number</b></label></td><td><input class='ui-widget-content ui-corner-all' name="newCopy" id="copyPlateNewCopy" size="4" type="text" maxlength="2" VALUE="$newCopy" readonly /></td></tr>
	</table>
</form>
<script>
$('#dialog').dialog("option", "title", "Copy Plate");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Copy", click: function() { submitForm('copyPlate'); } }, { text: "Cancel", click: function() { closeDialog(); } } ] );
</script>