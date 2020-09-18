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

my $vectorId = param ('vectorId') || '';

my $vector = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$vector->execute($vectorId);
my @vector=$vector->fetchrow_array();

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$vectorId/$vectorId/g;
$html =~ s/\$vectorName/$vector[2]/g;
$html =~ s/\$vectorBarcode/$vector[7]/g;
$html =~ s/\$vectorDescription/$vector[8]/g;
$html =~ s/\$vectorCreator/$vector[9]/g;
$html =~ s/\$vectorCreationDate/$vector[10]/g;
print header;
print $html;

__DATA__
<form id="editVector" name="editVector" action="vectorSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="vectorId" id="editVectorId" type="hidden" value="$vectorId" />
	<table>
	<tr><td style='text-align:right'><label for="editVectorName"><b>Vector Name</b></label></td><td><input class='ui-widget-content ui-corner-all' name="name" id="editVectorName" size="40" type="text" maxlength="32" value="$vectorName"/><br><img alt='$vectorBarcode' src='barcode.cgi?code=$vectorBarcode'/><sup class='ui-state-disabled'>by $vectorCreator on $vectorCreationDate</sup></td></tr>
	<tr><td style='text-align:right'><label for="editVectorDescription"><b>Description</b></label><br>(Sequence)</td><td><textarea class='ui-widget-content ui-corner-all' name="description" id="editVectorDescription" cols="50" rows="10">$vectorDescription</textarea></td></tr>
	</table>
</form>
<script>
$('#dialog').dialog("option", "title", "Edit Vector");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('editVector'); } }, { text: "Delete", click: function() { deleteItem($vectorId); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
</script>