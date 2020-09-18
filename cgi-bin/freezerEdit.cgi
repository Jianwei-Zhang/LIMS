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

my $freezerId = param ('freezerId') || '';

my $freezer = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$freezer->execute($freezerId);
my @freezer=$freezer->fetchrow_array();
my $freezerRoomId;
my $allRoom=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'room' ORDER BY o");
$allRoom->execute();
while (my @allRoom = $allRoom->fetchrow_array())
{
	$freezerRoomId .= ($freezer[3] == $allRoom[0]) ? "<option value='$allRoom[0]' selected>Room $allRoom[2]</option>"
		: "<option value='$allRoom[0]'>Room $allRoom[2]</option>";
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$freezerId/$freezerId/g;
$html =~ s/\$freezerName/$freezer[2]/g;
$html =~ s/\$freezerRoomId/$freezerRoomId/g;
$html =~ s/\$freezerX/$freezer[4]/g;
$html =~ s/\$freezerY/$freezer[5]/g;
$html =~ s/\$freezerZ/$freezer[6]/g;
$html =~ s/\$freezerBarcode/$freezer[7]/g;
$html =~ s/\$freezerDescription/$freezer[8]/g;

print header;
print $html;

__DATA__
<form id="editFreezer" name="editFreezer" action="freezerSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="freezerId" id="editFreezerId" type="hidden" value="$freezerId" />
	<table>
	<tr><td style='text-align:right'><label for="editFreezerName"><b>Freezer Name</b></label></td><td><input class='ui-widget-content ui-corner-all' name="name" id="editFreezerName" size="20" type="text" maxlength="32" value="$freezerName"/> in <select class='ui-widget-content ui-corner-all' name="roomId" id="editFreezerRoomId">$freezerRoomId</select><br><img alt='$freezerBarcode' src='barcode.cgi?code=$freezerBarcode'/></td></tr>
	<tr><td style='text-align:right'><label for="editFreezerX"><b>Space</b></label></td><td><label for="editFreezerX">X:</label><input class='ui-widget-content ui-corner-all' name="x" id="editFreezerX" size="4" type="text" maxlength="2" VALUE="$freezerX" /> <label for="editFreezerY">Y:</label><input class='ui-widget-content ui-corner-all' name="y" id="editFreezerY" size="4" type="text" maxlength="2" VALUE="$freezerY" /> <label for="editFreezerZ">Z:</label><input class='ui-widget-content ui-corner-all' name="z" id="editFreezerZ" size="4" type="text" maxlength="2" VALUE="$freezerZ" /></td></tr>
	<tr><td style='text-align:right'><label for="editFreezerDescription"><b>Description</b></label></td><td><textarea class='ui-widget-content ui-corner-all' name="description" id="editFreezerDescription" cols="50" rows="10">$freezerDescription</textarea></td></tr>
	</table>
</form>
<script>
$('#dialog').dialog("option", "title", "Edit Freezer");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('editFreezer'); } }, { text: "Delete", click: function() { deleteItem($freezerId); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
</script>