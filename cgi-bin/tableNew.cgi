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

my $type = param ('type') || '';
my $parentId = param ('parentId') || '';
my $refresh = param ('refresh') || 'menu';

my $parent=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$parent->execute($parentId);
my @parent = $parent->fetchrow_array();	

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$type/$type/g;
$html =~ s/\$parentId/$parentId/g;
$html =~ s/\$parentName/$parent[2]/g;
$html =~ s/\$refresh/$refresh/g;
$html =~ s/\$\$/$$/g;

print header;
print $html;

__DATA__
<form id="newTable" name="newTable" action="tableSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="parentId" id="newTableParentId" type="hidden" value="$parentId" />
	<input name="type" id="newTableType" type="hidden" value="$type" />
	<input name="refresh" id="newTableRefresh" type="hidden" value="$refresh" />
	<table>
	<tr><td style='text-align:right' rowspan='2'><label for="newTableFile"><b>New $type file</b></label></td><td><input name="tableFile" id="newTableFile" type="file" /><br>(in Tab Delimited Text format with header)</td></tr>
	<tr><td>or <input name="tableFilePath" id="newTableFilePath" type="text" />(On-server file name with full path)</td></tr>
	<tr><td><label for="newTableIdColumn">Column for Name </label></td><td><input type="text" id="newTableIdColumn" name="idColumn" value="1"></td></tr>
	<tr>
		<td colspan='2'><div class="ui-state-error-text">The new $type file will be applied to</div>
			<div id="newTableReplace">
			<input type="radio" id="newTableReplaceRadio1" name="replace" value="1"><label for="newTableReplaceRadio1">replace the existing data set</label>
			<input type="radio" id="newTableReplaceRadio2" name="replace" value="0" checked="checked"><label for="newTableReplaceRadio2">append to the existing data set</label>
			</div>
		</td>
	</tr>
	</table>
</form>
<script>
$( "#newTableReplace" ).buttonset();
$('#dialog').dialog("option", "title", "New $type in $parentName");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('newTable'); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
</script>