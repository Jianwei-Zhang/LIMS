#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use JSON::XS;
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

my %status = (
	0=>'Initialized',
	1=>'Started',
	2=>'Completed'
	);
my $maxSmrtcell = 16;

my $smrtrunId = param ('smrtrunId') || '';
my $smrtrun = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$smrtrun->execute($smrtrunId);
my @smrtrun=$smrtrun->fetchrow_array();
my $smrtrunDetails = decode_json $smrtrun[8];
$smrtrunDetails->{'comments'} = '' unless (exists $smrtrunDetails->{'comments'});

my $smrtrunStatus = '';
foreach (sort {$a <=> $b} keys %status)
{
	$smrtrunStatus .= ($smrtrun[7] eq $_) ? "<option value='$_' selected>$status{$_}</option>" : "<option value='$_'>$status{$_}</option>";
}

my $smrtrunSmrtcell = '';
for (my $i = 1;$i <= $maxSmrtcell;$i++)
{
	$smrtrunSmrtcell .= ($smrtrun[3] eq $i) ? "<option value='$i' selected>$i</option>" : "<option value='$i'>$i</option>";
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$smrtrunId/$smrtrunId/g;
$html =~ s/\$smrtrunName/$smrtrun[2]/g;
$html =~ s/\$smrtrunSmrtcell/$smrtrunSmrtcell/g;
$html =~ s/\$smrtrunStatus/$smrtrunStatus/g;
$html =~ s/\$comments/$smrtrunDetails->{'comments'}/g;
$html =~ s/\$smrtrunCreator/$smrtrun[9]/g;
$html =~ s/\$smrtrunEnteredDate/$smrtrun[10]/g;

print header;
print $html;

__DATA__
<form id="editSmrtrun" name="editSmrtrun" action="smrtrunSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="smrtrunId" id="editSmrtrunId" type="hidden" value="$smrtrunId" />
	<table>
		<tr><td style='text-align:right;white-space: nowrap;'><label for="editSmrtrunName"><b>SMRT Run Name</b></label></td>
			<td><input class='ui-widget-content ui-corner-all' name="name" id="editSmrtrunName" placeholder="Run Name" size="15" type="text" maxlength="32" value="$smrtrunName" /><sup class='ui-state-disabled'>Last changed by $smrtrunCreator on $smrtrunEnteredDate</sup></td>
		</tr>
		<tr><td style='text-align:right;white-space: nowrap;'><label for="editSmrtrunSmrtcell"><b>SMRTCells</b></label></td>
			<td><select class='ui-widget-content ui-corner-all' name="smrtcell" id="editSmrtrunSmrtcell">$smrtrunSmrtcell</select></td>
		</tr>
		<tr><td style='text-align:right;white-space: nowrap;'><label for="editSmrtrunStatus"><b>Status</b></label></td>
			<td><select class='ui-widget-content ui-corner-all' name="status" id="editSmrtrunStatus">$smrtrunStatus</select></td>
		</tr>
		<tr><td style='text-align:right;white-space: nowrap;'><label for="editSmrtrunComments"><b>Comments</b></label></td>
			<td><textarea class='ui-widget-content ui-corner-all' name="comments" id="editSmrtrunComments" cols="40" rows="5" placeholder="Please write anything relevant to the run that may not be included on this form">$comments</textarea></td>
		</tr>
	</table>
</form>
<script>
$('#dialog').dialog("option", "title", "Edit SMRT Run $smrtrunName");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('editSmrtrun'); } }, { text: "Delete", click: function() { deleteItem($smrtrunId); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
</script>
