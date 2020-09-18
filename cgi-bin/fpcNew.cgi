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
<form id="newFpc" name="newFpc" action="fpcLoad.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
<h3>Load FPC for Library $libraryName<sup class='ui-state-disabled'>(Library id: $libraryId)</sup></h3>
	<input name="libraryId" id="newFpclibraryId" type="hidden" value="$libraryId" />
	<table>
	<tr><td style='text-align:right'><label for="newFpcFile"><b>FPC File</b></label></td><td><input name="fpcFile" id="newFpcFile" type="file" /></td></tr>
	<tr><td style='text-align:left' colspan='2'><label for="newFpcReplace"><b>Are you going to update any existing FPC in database?</b></label></td></tr>
	<tr><td style='text-align:right'></td>
		<td>
		<div id="newFpcReplace">
			<input type="radio" id="newFpcReplaceRadio2" name="replace" value="1" checked="checked"><label for="newFpcReplaceRadio2">Yes</label>
			<input type="radio" id="newFpcReplaceRadio1" name="replace" value="0"><label for="newFpcReplaceRadio1">No</label>
		</div>
	</td></tr>
	</table>
</form>
<script>
$('#dialog').dialog("option", "title", "Load FPC");
$( "#newFpcReplace" ).buttonset();
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Upload", click: function() { submitForm('newFpc'); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
</script>