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

my $roomId = param ('roomId') || '';

my $room = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$room->execute($roomId);
my @room=$room->fetchrow_array();

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$roomId/$roomId/g;
$html =~ s/\$roomName/$room[2]/g;
$html =~ s/\$roomX/$room[4]/g;
$html =~ s/\$roomY/$room[5]/g;
$html =~ s/\$roomZ/$room[6]/g;
$html =~ s/\$roomDescription/$room[8]/g;

print header;
print $html;

__DATA__
<form id="editRoom" name="editRoom" action="roomSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="roomId" id="editRoomId" type="hidden" value="$roomId" />
	<table>
	<tr><td style='text-align:right'><label for="editRoomName"><b>Room Name</b></label></td><td><input class='ui-widget-content ui-corner-all' name="name" id="editRoomName" size="40" type="text" maxlength="32" value="$roomName"/></td></tr>
	<tr><td style='text-align:right'><label for="editRoomX"><b>Space</b></label></td><td><label for="editRoomX">X:</label><input class='ui-widget-content ui-corner-all' name="x" id="editRoomX" size="4" type="text" maxlength="2" VALUE="$roomX" /> <label for="editRoomY">Y:</label><input class='ui-widget-content ui-corner-all' name="y" id="editRoomY" size="4" type="text" maxlength="2" VALUE="$roomY" /> <label for="editRoomZ">Z:</label><input class='ui-widget-content ui-corner-all' name="z" id="editRoomZ" size="4" type="text" maxlength="2" VALUE="$roomZ" /></td></tr>
	<tr><td style='text-align:right'><label for="editRoomDescription"><b>Description</b></label></td><td><textarea class='ui-widget-content ui-corner-all' name="description" id="editRoomDescription" cols="50" rows="10">$roomDescription</textarea></td></tr>
	</table>
</form>
<script>
$('#dialog').dialog("option", "title", "Edit Room");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('editRoom'); } }, { text: "Delete", click: function() { deleteItem($roomId); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
</script>