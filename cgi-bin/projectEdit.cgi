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

my $projectId = param ('projectId') || '';
my $project = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$project->execute($projectId);
my @project=$project->fetchrow_array();

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$projectId/$projectId/g;
$html =~ s/\$projectName/$project[2]/g;
$html =~ s/\$projectDescription/$project[8]/g;
print header;
print $html;

__DATA__
<form id="editProject" name="editProject" action="projectSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="projectId" id="editProjectId" type="hidden" value="$projectId" />
	<table>
	<tr><td style='text-align:right'><label for="editProjectName"><b>Project Name</b></label></td><td><input class='ui-widget-content ui-corner-all' name="name" id="editProjectName" size="40" type="text" maxlength="32" value="$projectName"/></td></tr>
	<tr><td style='text-align:right'><label for="editProjectDescription"><b>Description</b></label></td><td><textarea class='ui-widget-content ui-corner-all' name="description" id="editProjectDescription" cols="50" rows="10">$projectDescription</textarea></td></tr>
	</table>
</form>
<script>
$('#dialog').dialog("option", "title", "Edit Project");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('editProject'); } }, { text: "Delete", click: function() { deleteItem($projectId); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
</script>