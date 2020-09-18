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
my $library=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$library->execute($libraryId);
my @library = $library->fetchrow_array();

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$libraryName/$library[2]/g;
$html =~ s/\$libraryId/$libraryId/g;

print header;
print $html;

__DATA__
<form id="newBes" name="newBes" action="besSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
<h3>Load BAC End Sequences for Library $libraryName<sup class='ui-state-disabled'>(Library id: $libraryId)</sup></h3>
	<input name="libraryId" id="newBesLibraryId" type="hidden" value="$libraryId" />
	<table>
	<tr><td style='text-align:right'><label for="newBesFile"><b>BES File</b></label></td><td><input name="besFile" id="newBesFile" type="file" />(in FASTA format)</td></tr>
	<tr><td></td><td>or <input name="besFilePath" id="newBesFilePath" type="text" />(On-server file name with full path)</td></tr>
	<tr><td style='text-align:right'></td><td>Sequence title must be formated like 'cloneName.direction'<br>e.g., $libraryName0001A01.f, $libraryName0001A01.r, ... </td></tr>
	<tr><td style='text-align:left' colspan='2'><label for="newBesReplace"><b>Are you going to update any existing BES in database?</b></label></td></tr>
	<tr><td style='text-align:right'></td>
		<td>
		<div id="newBesReplace">
			<input type="radio" id="newBesReplaceRadio2" name="replace" value="1" checked="checked"><label for="newBesReplaceRadio2">Yes</label>
			<input type="radio" id="newBesReplaceRadio1" name="replace" value="0"><label for="newBesReplaceRadio1">No</label>
		</div>
	</td></tr>
	</table>
</form>
<script>
$('#dialog').dialog("option", "title", "Load BES");
$( "#newBesReplace" ).buttonset();
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Upload", click: function() { submitForm('newBes'); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
</script>