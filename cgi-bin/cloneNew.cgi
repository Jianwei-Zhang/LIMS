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
$html =~ s/\$libraryId/$libraryId/g;
$html =~ s/\$libraryName/$library[2]/g;
$html =~ s/\$nofPlates/$library[3]/g;

print header;
print $html;

__DATA__
<form id="newClone" name="newClone" action="cloneSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="libraryId" id="libraryId" type="hidden" value="$libraryId" />
	This operation will generate a new clone list for $nofPlates plate(s) in $libraryName and may take several minutes.<br>
	<input type="checkbox" id="reformatSourceCloneName" name="reformat" value="1" checked="checked"><label for="reformatSourceCloneName">Reformat Source Clone Names (for rearraying library only).</label>
</form>
<script>
$('#dialog').dialog("option", "title", "Generate Clone List");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Generate", click: function() { submitForm('newClone'); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
</script>