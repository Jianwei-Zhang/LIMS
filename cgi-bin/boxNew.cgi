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
		$boxFreezerId .= ($freezerId == $freezerInRoom[0]) ? "<option value='$freezerInRoom[0]' selected>Freezer $freezerInRoom[2]</option>"
			: "<option value='$freezerInRoom[0]'>Freezer $freezerInRoom[2]</option>";
	}
	$boxFreezerId .= "</optgroup>";
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$boxFreezerId/$boxFreezerId/g;
print header;
print $html;

__DATA__
<form id="newBox" name="newBox" action="boxSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<table>
	<tr><td style='text-align:right'><label for="newBoxName"><b>Box Name</b></label></td><td><input class='ui-widget-content ui-corner-all' name="name" id="newBoxName" size="20" type="text" maxlength="32" /> in <select class='ui-widget-content ui-corner-all' name="freezerId" id="newBoxFreezerId">$boxFreezerId</select></td></tr>
	<tr><td style='text-align:right'><label for="newBoxX"><b>Space</b></label></td><td><label for="newBoxX">X:</label><input class='ui-widget-content ui-corner-all' name="x" id="newBoxX" size="4" type="text" maxlength="2" VALUE="10" /> <label for="newBoxY">Y:</label><input class='ui-widget-content ui-corner-all' name="y" id="newBoxY" size="4" type="text" maxlength="2" VALUE="1" /> <label for="newBoxZ">Z:</label><input class='ui-widget-content ui-corner-all' name="z" id="newBoxZ" size="4" type="text" maxlength="2" VALUE="1" /></td></tr>
	<tr><td style='text-align:right'><label for="newBoxDescription"><b>Description</b></label></td><td><textarea class='ui-widget-content ui-corner-all' name="description" id="newBoxDescription" cols="50" rows="10"></textarea></td></tr>
	</table>
</form>
<script>
$('#dialog').dialog("option", "title", "New Box");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('newBox'); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
</script>