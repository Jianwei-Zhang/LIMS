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

my $fpcId = param ('fpcId') || '';
if ($fpcId)
{
	my $fpc=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$fpc->execute($fpcId);
	my @fpc = $fpc->fetchrow_array();
	
	undef $/;# enable slurp mode
	my $html = <DATA>;
	$html =~ s/\$fpcId/$fpcId/g;
	$html =~ s/\$fpcName/$fpc[2]/g;
	$html =~ s/\$fpcVersion/$fpc[3]/g;
}
else
{
	print header(-type=>'text/html',-status=>'402 Invalid operation');
	print <<END;
<script>
	errorPop("Not a valid operation!");
</script>	
END
	exit;
}


print header;
print $html;

__DATA__
<form id="markClone" name="markClone" action="fpcMarkClone.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
<h3>Mark clones in '$fpcName v.$fpcVersion' as</h3>
	<input name="fpcId" id="fpcId" type="hidden" value="$fpcId" />
	<table>
	<tr><td style='text-align:right'></td>
		<td>
		<div id="markCloneAs">
			<input type="radio" id="markCloneAsRadio2" name="markAs" value="1" checked="checked"><label for="markCloneAsRadio2"><b>MTP Clones</b></label>
			<input type="radio" id="markCloneAsRadio1" name="markAs" value="2"><label for="markCloneAsRadio1"><b>Highlighted Clones</b></label>
		</div>
	</td></tr>
	<tr><td style='text-align:right'><sub id="markClones_count" style="display:none"></sub></td><td><textarea class='ui-widget-content ui-corner-all word_count' name="markClones" id="markClones" cols="50" rows="8"></textarea></td></tr>
	<tr><td style='text-align:left' colspan='2'><label for="markCloneReplace">Are you going to clear all marked clones before marking the above clones?</label></td></tr>
	<tr><td style='text-align:right'></td>
		<td>
		<div id="markCloneReplace">
			<input type="radio" id="markCloneReplaceRadio2" name="replace" value="1"><label for="markCloneReplaceRadio2">Yes</label>
			<input type="radio" id="markCloneReplaceRadio1" name="replace" value="0" checked="checked"><label for="markCloneReplaceRadio1">No</label>
		</div>
	</td></tr>
	</table>
</form>
<script>
$('#dialog').dialog("option", "title", "Mark Clones");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Submit", click: function() { submitForm('markClone'); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
$( "#markCloneAs" ).buttonset();
$( "#markCloneReplace" ).buttonset();
wordCount('clone');
</script>
