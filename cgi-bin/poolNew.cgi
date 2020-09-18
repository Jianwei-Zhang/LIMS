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

my @cloneName = param ('cloneName');
my $libraryId = param ('libraryId') || '';
my $cloneList = join "\n", @cloneName;

my $library = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$library->execute($libraryId);
my @library=$library->fetchrow_array();

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$seqLibraryName/$library[2]/g;
$html =~ s/\$cloneList/$cloneList/g;
$html =~ s/\$libraryId/$libraryId/g;

print header;
print $html;

__DATA__
<form id="newPool" name="newPool" action="poolSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="libraryId" id="libraryId" type="hidden" value="$libraryId" />
	<table>
	<tr><td style='text-align:right'><label for="newPoolName"><b>Pool Name</b></label></td><td colspan="2">$seqLibraryName<input class='ui-widget-content ui-corner-all' name="name" id="newPoolName" size="40" type="text" maxlength="32" /></td></tr>
	<tr><td style='text-align:right'></td>
	<td><label for="newPoolClones"><b>Clone List</b></label> (1 clone/line) <sub id="newPoolClones_count" style="display:none"></sub><br><textarea class='ui-widget-content ui-corner-all word_count' name="poolClones" id="newPoolClones" cols="30" rows="10">$cloneList</textarea></td>
	<td><label for="newPoolJobs"><b>Assembly</b></label><br><textarea class='ui-widget-content ui-corner-all' name="poolJobs" id="newPoolJobs" cols="20" rows="10"></textarea></td></tr>
	<tr><td style='text-align:right'><label for="newPoolComments"><b>Comments</b></label></td><td colspan="2"><textarea class='ui-widget-content ui-corner-all' name="comments" id="newPoolComments" cols="50" rows="4"></textarea></td></tr>
	</table>
</form>
<script>
$('#dialog').dialog("option", "title", "New Pool");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('newPool'); } }, { text: "Cancel", click: function() { closeDialog(); } } ] );
wordCount('clone');
</script>