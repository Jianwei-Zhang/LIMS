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

my @itemId = param ('itemId');
my $items;
for(@itemId)
{
	my $item = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$item->execute($_);
	my @item=$item->fetchrow_array();
	if($item[1] eq "plate")
	{
		my $plateToLibrary=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$plateToLibrary->execute($item[6]);
		my @plateToLibrary = $plateToLibrary->fetchrow_array();
		$items->{$item[0]} = "$plateToLibrary[2]-$item[2].$item[4]";
	}
	elsif($item[1] eq "sample")
	{
		my $sampleToService=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$sampleToService->execute($item[6]);
		my @sampleToService = $sampleToService->fetchrow_array();
		$items->{$item[0]} = "$sampleToService[2] Sample-$item[2]";
	}
	elsif($item[1] eq "paclib")
	{
		my $paclibToSample=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$paclibToSample->execute($item[6]);
		my @paclibToSample = $paclibToSample->fetchrow_array();
		my $sampleToService=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$sampleToService->execute($paclibToSample[6]);
		my @sampleToService = $sampleToService->fetchrow_array();
		$items->{$item[0]} = "$sampleToService[2] > $paclibToSample[2] > Paclib-$item[2]";
	}
	else
	{
		$items->{$item[0]} = "Misc-$item[2]";
	}
}
my $itemList = 'Please check items first!';
$itemList = checkbox_group(-name=>'items',
-values=>[sort {$a <=> $b} keys %{$items}],
-default=>[sort keys %{$items}],
-labels=>\%{$items},
-columns=>3) if (keys %{$items});

my $plateLocation;
my $room=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'room' ORDER BY o");
$room->execute();
while (my @room = $room->fetchrow_array())
{
	$plateLocation .= "<optgroup label='Room $room[2]'>";
	my $freezerInRoom=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'freezer' AND o = ? ORDER BY name");
	$freezerInRoom->execute($room[0]);
	while (my @freezerInRoom = $freezerInRoom->fetchrow_array())
	{
		$plateLocation .= "<optgroup label='--Freezer $freezerInRoom[2]'>";
		my $boxInFreezer=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'box' AND o = ? ORDER BY name");
		$boxInFreezer->execute($freezerInRoom[0]);
		while (my @boxInFreezer = $boxInFreezer->fetchrow_array())
		{
			$plateLocation .= "<option value='$boxInFreezer[0]'>Box $boxInFreezer[2]</option>";
		}
		$plateLocation .= "</optgroup>";
	}
	$plateLocation .= "</optgroup>";
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$itemList/$itemList/g;
$html =~ s/\$plateLocation/$plateLocation/g;

print header;
print $html;

__DATA__
<form id="depotItem" name="depotItem" action="itemDepot.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<h3 class='ui-state-error-text'>Where are you going to depot the below items to?</h3>
	<table width='100%'>
	<tr><td style='text-align:right'><label for="depotItemList"><b>Items</b></label></td><td>$itemList</td></tr>
	<tr><td style='text-align:right'><label for="depotItemLocation"><b>Container</b></label></td><td><select class='ui-widget-content ui-corner-all' name="location" id="depotItemLocation">$plateLocation</select></td></tr>
	</table>
</form>
<script>
$('#dialog').dialog("option", "title", "Depot Item");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Depot", click: function() { submitForm('depotItem'); } }, { text: "Cancel", click: function() { closeDialog(); } } ] );
</script>