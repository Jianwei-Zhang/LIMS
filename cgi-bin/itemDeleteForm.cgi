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
	elsif($item[1] eq "library")
	{
		$items->{$item[0]} = "Library-$item[2]";
	}
	elsif($item[1] eq "vector")
	{
		$items->{$item[0]} = "Vector-$item[2]";
	}
	elsif($item[1] eq "sample")
	{
		my $sampleToService=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$sampleToService->execute($item[6]);
		my @sampleToService = $sampleToService->fetchrow_array();
		$items->{$item[0]} = "$sampleToService[2] Sample-$item[2]";
	}
	elsif($item[1] eq "room")
	{
		$items->{$item[0]} = "Room-$item[2]";
	}
	elsif($item[1] eq "freezer")
	{
		$items->{$item[0]} = "Freezer-$item[2]";
	}
	elsif($item[1] eq "box")
	{
		$items->{$item[0]} = "Box-$item[2]";
	}
	elsif($item[1] eq "pool")
	{
		$items->{$item[0]} = "Pool-$item[2]";
	}
	elsif($item[1] eq "project")
	{
		$items->{$item[0]} = "Project-$item[2]";
	}
	elsif($item[1] eq "asbProject")
	{
		$items->{$item[0]} = "Assembly Project-$item[2]";
	}
	elsif($item[1] eq "assembly")
	{
		$items->{$item[0]} = "Assembly-$item[2]";
	}
	elsif($item[1] eq "genome")
	{
		$items->{$item[0]} = "Genome-$item[2]";
	}
	elsif($item[1] eq "sequence")
	{
		$items->{$item[0]} = "Sequence-$item[2]";
	}
	elsif($item[1] eq "dataset")
	{
		$items->{$item[0]} = "Dataset-$item[2]";
	}
	elsif($item[1] eq "record")
	{
		$items->{$item[0]} = "Record-$item[2]";
	}
	elsif($item[1] eq "fpc")
	{
		$items->{$item[0]} = "FPC-$item[2]";
	}
	elsif($item[1] eq "agp")
	{
		$items->{$item[0]} = "AGP-$item[2]";
	}
	elsif($item[1] eq "paclib")
	{
		my $paclibToSample=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$paclibToSample->execute($item[6]);
		my @paclibToSample = $paclibToSample->fetchrow_array();
		my $sampleToService=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$sampleToService->execute($paclibToSample[6]);
		my @sampleToService = $sampleToService->fetchrow_array();
		$items->{$item[0]} = "$sampleToService[2] > Sample $paclibToSample[2] > Paclib-$item[2]";
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

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$itemList/$itemList/g;

print header;
print $html;

__DATA__
<form id="deleteItem" name="deleteItem" action="itemDelete.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<h3 class='ui-state-error-text'>Are you sure to delete the below items?</h3>
	$itemList
</form>
<script>
$('#dialog').dialog("option", "title", "Delete Item");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Delete", click: function() { submitForm('deleteItem'); } }, { text: "Cancel", click: function() { closeDialog(); } } ] );
</script>