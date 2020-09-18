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
$pool[8] = 'none.' unless($pool[8]);

my $poolToLibrary = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$poolToLibrary->execute($pool[4]);
my @poolToLibrary = $poolToLibrary->fetchrow_array();

my $poolJobs = '';
my $jobId=$dbh->prepare("SELECT * FROM link WHERE parent = ? AND type LIKE 'poolJob' ORDER BY child DESC");
$jobId->execute($pool[0]);
while(my @jobId = $jobId->fetchrow_array())
{
	$poolJobs .= "<a onclick='closeDialog();openDialog(\"jobView.cgi?jobId=$jobId[1]\")'>$jobId[1]</a>\n";
}
$poolJobs = 'none.' unless($poolJobs);
my $poolClones = '';
my $poolClone=$dbh->prepare("SELECT clones.* FROM link,clones WHERE link.type LIKE 'poolClone' AND link.parent = ? AND clones.name LIKE link.child ORDER BY link.child");
$poolClone->execute($pool[0]);
while(my @poolClone = $poolClone->fetchrow_array())
{
	$poolClones .= ($poolClone[6]) ? 
		"<div style='float: left; margin-right: .7em;' title='Sequence available'><span class='ui-icon ui-icon-check' style='float: left; margin-right: 0;'></span><a onclick='closeDialog();openDialog(\"cloneView.cgi?cloneName=$poolClone[1]\")'>$poolClone[1]</a></div>" : 
		"<div class='ui-state-disabled' style='float: left; margin-right: .7em;' title='Sequence unavailable'><span class='ui-icon ui-icon-cancel' style='float: left; margin-right: 0;'></span><a onclick='closeDialog();openDialog(\"cloneView.cgi?cloneName=$poolClone[1]\")'>$poolClone[1]</a></div>";
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
<table>
	<tr><td style='text-align:right'><b>Pool Name</b></td><td>$seqLibraryName<a class='ui-state-content ui-corner-all' name="name" id="viewPoolName">$poolName</a><br><img alt='$poolBarcode' src='barcode.cgi?code=$poolBarcode'/><sup class='ui-state-disabled'>by $poolCreator on $poolCreationDate</sup></td></tr>
	<tr><td style='text-align:right'><b>Clone List</b></td><td>$poolClones</td></tr>
	<tr><td style='text-align:right'><b>Assembly</b></td><td>$poolJobs</td></tr>
	<tr><td style='text-align:right'><b>Comments</b></td><td>$poolComments</td></tr>
</table>
<script>
$('#dialog').dialog("option", "title", "View Pool");
$( "#dialog" ).dialog( "option", "buttons", [ { text: "Edit", click: function() {closeDialog(); openDialog("poolEdit.cgi?poolId=$poolId")} }, { text: "OK", click: function() {closeDialog(); } } ] );
</script>
