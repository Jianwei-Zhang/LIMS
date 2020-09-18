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

my $poolId = param ('poolId') || '';

my $pool = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$pool->execute($poolId);
my @pool=$pool->fetchrow_array();

my $poolToLibrary = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$poolToLibrary->execute($pool[4]);
my @poolToLibrary = $poolToLibrary->fetchrow_array();

my $poolJobs = '';
my $jobId=$dbh->prepare("SELECT * FROM link WHERE parent = ? AND type LIKE 'poolJob' ORDER BY child DESC");
$jobId->execute($pool[0]);
while(my @jobId = $jobId->fetchrow_array())
{
	$poolJobs .= "$jobId[1]\n";
}
my $poolClones = '';
my $poolClone=$dbh->prepare("SELECT * FROM link WHERE parent = ? AND type LIKE 'poolClone' ORDER BY child");
$poolClone->execute($pool[0]);
while(my @poolClone = $poolClone->fetchrow_array())
{
	$poolClones .= "$poolClone[1]\n";
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$poolId/$poolId/g;
$html =~ s/\$poolName/$pool[2]/g;
$html =~ s/\$seqLibraryName/$poolToLibrary[2]/g;
$html =~ s/\$seqLibrary/$pool[4]/g;
$html =~ s/\$poolBarcode/$pool[7]/g;
$html =~ s/\$poolClones/$poolClones/g;
$html =~ s/\$poolComments/$pool[8]/g;
$html =~ s/\$poolJobs/$poolJobs/g;
$html =~ s/\$poolCreator/$pool[9]/g;
$html =~ s/\$poolCreationDate/$pool[10]/g;
print header;
print $html;

__DATA__
<form id="editPool" name="editPool" action="poolSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="poolId" id="editPoolId" type="hidden" value="$poolId" />
	<input name="seqLibrary" id="seqLibrary" type="hidden" value="$seqLibrary" />
	<table>
	<tr><td style='text-align:right'><label for="editPoolName"><b>Pool Name</b></label></td><td colspan="2">$seqLibraryName<input class='ui-widget-content ui-corner-all' name="name" id="editPoolName" size="40" type="text" maxlength="32" value="$poolName"/><br><img alt='$poolBarcode' src='barcode.cgi?code=$poolBarcode'/><sup class='ui-state-disabled'>by $poolCreator on $poolCreationDate</sup></td></tr>
	<tr><td style='text-align:right'></td>
	<td><label for="editPoolClones"><b>Clone List</b></label> (1 clone/line) <sub id="editPoolClones_count" style="display:none"></sub><br><textarea class='ui-widget-content ui-corner-all word_count' name="poolClones" id="editPoolClones" cols="30" rows="10">$poolClones</textarea></td>
	<td><label for="editPoolJobs"><b>Assembly</b></label><br><textarea class='ui-widget-content ui-corner-all' name="poolJobs" id="editPoolJobs" cols="20" rows="10">$poolJobs</textarea></td></tr>
	<tr><td style='text-align:right'><label for="editPoolComments"><b>Comments</b></label></td><td colspan="2"><textarea class='ui-widget-content ui-corner-all' name="comments" id="editPoolComments" cols="50" rows="4">$poolComments</textarea></td></tr>
	</table>
</form>
<script>
$('#dialog').dialog("option", "title", "Edit Pool");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('editPool'); } }, { text: "Delete", click: function() { deleteItem($poolId); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
wordCount('clone');
</script>
