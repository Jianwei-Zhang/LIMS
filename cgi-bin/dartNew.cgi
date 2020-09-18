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

my $givenParentId = param ('parentId') || '';
my $parentId = "<option class='ui-state-error-text' value='0'>None</option>";
my $parentList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'genebank'");# ORDER BY name
$parentList->execute();
if($parentList->rows > 0)
{
	my $parentListResult;
	while (my @parentList = $parentList->fetchrow_array())
	{
		@{$parentListResult->{$parentList[2]}} = @parentList;
	}
	foreach (sort {uc ($a) cmp uc($b)} keys %$parentListResult)
	{
		my @parentList = @{$parentListResult->{$_}};
		my $genebankDetails = decode_json $parentList[8];
		$genebankDetails->{'comments'} = escapeHTML($genebankDetails->{'comments'});
		$parentId .= ($givenParentId eq $parentList[0]) ? "<option value='$parentList[0]' title='$genebankDetails->{'comments'}' selected>$parentList[1]: $parentList[2]</option>" : "<option value='$parentList[0]' title='$genebankDetails->{'comments'}'>$parentList[1]: $parentList[2]</option>";
	}
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$parentId/$parentId/g;

print header;
print $html;

__DATA__
<form id="newDart" name="newDart" action="dartSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<table>
	<tr><td style='text-align:right'><label for="newDartName"><b>DArTseq Name</b></label></td><td><input class='ui-widget-content ui-corner-all' name="name" id="newDartName" size="40" type="text" maxlength="32" /></td></tr>
	<tr><td style='text-align:right' rowspan='2'><label for="newDartFile"><b>DArTseq File</b></label></td><td><input name="dartFile" id="newDartFile" type="file" />(in Tab Delimited Text format)</td></tr>
	<tr><td>or <input name="dartFilePath" id="newDartFilePath" type="text" />(On-server file name with full path)</td></tr>
	<tr><td style='text-align:right'><label for="newDartGenebankId"><b>Link to</b></label></td><td><select class='ui-widget-content ui-corner-all' name='genebankId' id='newDartGenebankId'>$parentId</select></td></tr>
	<tr><td style='text-align:right'><label for="newDartDescription"><b>Description</b></label></td><td><textarea class='ui-widget-content ui-corner-all' name="description" id="newDartDescription" cols="50" rows="10"></textarea></td></tr>
	</table>
</form>
<script>
$('#dialog').dialog("option", "title", "New DArTseq");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('newDart'); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
</script>