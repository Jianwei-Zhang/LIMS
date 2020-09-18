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

my $boxId = param ('boxId') || '';

my $box = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$box->execute($boxId);
my @box=$box->fetchrow_array();

my $boxFreezerId;
my $room=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'room' ORDER BY o");
$room->execute();
while (my @room = $room->fetchrow_array())
{
	$boxFreezerId .= "<optgroup label='Room $room[2]'>";
	my $freezerInRoom=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'freezer' AND o = ? ORDER BY name");
	$freezerInRoom->execute($room[0]);
	while (my @freezerInRoom = $freezerInRoom->fetchrow_array())
	{
		$boxFreezerId .= ($box[3] == $freezerInRoom[0]) ? "<option value='$freezerInRoom[0]' selected>Freezer $freezerInRoom[2]</option>"
			: "<option value='$freezerInRoom[0]'>Freezer $freezerInRoom[2]</option>";
	}
	$boxFreezerId .= "</optgroup>";
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$boxId/$boxId/g;
$html =~ s/\$boxName/$box[2]/g;
$html =~ s/\$boxFreezerId/$boxFreezerId/g;
$html =~ s/\$boxX/$box[4]/g;
$html =~ s/\$boxY/$box[5]/g;
$html =~ s/\$boxZ/$box[6]/g;
$html =~ s/\$boxBarcode/$box[7]/g;
$html =~ s/\$boxDescription/$box[8]/g;
$html =~ s/\$boxCreator/$box[9]/g;
$html =~ s/\$boxCreationDate/$box[10]/g;

print header;
print $html;

__DATA__
<form id="editBox" name="editBox" action="boxSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="boxId" id="editBoxId" type="hidden" value="$boxId" />
	<table>
	<tr><td style='text-align:right'><label for="editBoxName"><b>Box Name</b></label></td><td><input class='ui-widget-content ui-corner-all' name="name" id="editBoxName" size="20" type="text" maxlength="32" value="$boxName"/> in <select class='ui-widget-content ui-corner-all' name="freezerId" id="editBoxFreezerId">$boxFreezerId</select><br><img alt='$boxBarcode' src='barcode.cgi?code=$boxBarcode'/><sup class='ui-state-disabled'>by $boxCreator on $boxCreationDate</sup></td></tr>
	<tr><td style='text-align:right'><label for="editBoxX"><b>Space</b></label></td><td><label for="editBoxX">X:</label><input class='ui-widget-content ui-corner-all' name="x" id="editBoxX" size="4" type="text" maxlength="2" VALUE="$boxX" /> <label for="editBoxY">Y:</label><input class='ui-widget-content ui-corner-all' name="y" id="editBoxY" size="4" type="text" maxlength="2" VALUE="$boxY" /> <label for="editBoxZ">Z:</label><input class='ui-widget-content ui-corner-all' name="z" id="editBoxZ" size="4" type="text" maxlength="2" VALUE="$boxZ" /></td></tr>
	<tr><td style='text-align:right'><label for="editBoxDescription"><b>Description</b></label></td><td><textarea class='ui-widget-content ui-corner-all' name="description" id="editBoxDescription" cols="50" rows="10">$boxDescription</textarea></td></tr>
	</table>
</form>
<script>
$('#dialog').dialog("option", "title", "Edit Box");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('editBox'); } }, { text: "Delete", click: function() { deleteItem($boxId); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
</script>