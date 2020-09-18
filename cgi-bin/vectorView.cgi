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
	<table>
	<tr><td style='text-align:right'><label for="viewVectorName"><b>Vector Name</b></label></td><td>$vectorName<br><img alt='$vectorBarcode' src='barcode.cgi?code=$vectorBarcode'/><sup class='ui-state-disabled'>by $vectorCreator on $vectorCreationDate</sup></td></tr>
	<tr><td style='text-align:right'><label for="viewVectorDescription"><b>Description</b></label><br>(Sequence)</td><td><textarea class='ui-widget-content ui-corner-all' name="description" id="viewVectorDescription" cols="50" rows="10" readonly="readonly">$vectorDescription</textarea></td></tr>
	</table>
<script>
$('#dialog').dialog("option", "title", "View Vector $vectorName");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Edit", click: function() { closeDialog();openDialog("vectorEdit.cgi?vectorId=$vectorId"); } }, { text: "OK", click: function() {closeDialog(); } } ] );
</script>