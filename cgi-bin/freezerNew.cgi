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
my $freezerRoomId = '';
my $room=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'room' ORDER BY o");
$room->execute();
while (my @room = $room->fetchrow_array())
{
	if($roomId == $room[0])
	{
		$freezerRoomId .= "<option value='$room[0]' selected>Room $room[2]</option>";
	}
	else
	{
		$freezerRoomId .= "<option value='$room[0]'>Room $room[2]</option>";
	}
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$freezerRoomId/$freezerRoomId/g;

print header;
print $html;

__DATA__
<form id="newFreezer" name="newFreezer" action="freezerSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<table>
	<tr><td style='text-align:right'><label for="newFreezerName"><b>Freezer Name</b></label></td><td><input class='ui-widget-content ui-corner-all' name="name" id="newFreezerName" size="20" type="text" maxlength="32" /> in <select class='ui-widget-content ui-corner-all' name="roomId" id="newFreezerRoomId">$freezerRoomId</select></td></tr>
	<tr><td style='text-align:right'><label for="newFreezerX"><b>Space</b></label></td><td><label for="newFreezerX">X:</label><input class='ui-widget-content ui-corner-all' name="x" id="newFreezerX" size="4" type="text" maxlength="2" VALUE="10" /> <label for="newFreezerY">Y:</label><input class='ui-widget-content ui-corner-all' name="y" id="newFreezerY" size="4" type="text" maxlength="2" VALUE="1" /> <label for="newFreezerZ">Z:</label><input class='ui-widget-content ui-corner-all' name="z" id="newFreezerZ" size="4" type="text" maxlength="2" VALUE="1" /></td></tr>
	<tr><td style='text-align:right'><label for="newFreezerDescription"><b>Description</b></label></td><td><textarea class='ui-widget-content ui-corner-all' name="description" id="newFreezerDescription" cols="50" rows="10"></textarea></td></tr>
	</table>
</form>
<script>
$('#dialog').dialog("option", "title", "New Freezer");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('newFreezer'); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
</script>