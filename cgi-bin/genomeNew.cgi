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

my $asbProjectId = param ('asbProjectId') || '';
$asbProjectId = "<input name='asbProjectId' id='newGenomeAsbProjectId' type='hidden' value='$asbProjectId' />" if ($asbProjectId);
my $forAssembly = ($asbProjectId) ? " checked=checked": "";

my $objectComponent;
$objectComponent->{0} = "Unknown";
$objectComponent->{1} = "Chr-Seq";
$objectComponent->{2} = "Ctg-Seq";

my $libraryId = "<option class='ui-state-error-text' value='0'>None</option>";
my $libraryList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'library'");# ORDER BY name
$libraryList->execute();
if ($libraryList->rows > 0)
{
	my $libraryListResult;
	while (my @libraryList = $libraryList->fetchrow_array())
	{
		@{$libraryListResult->{$libraryList[2]}} = @libraryList;
	}
	foreach (sort {uc ($a) cmp uc($b)} keys %$libraryListResult)
	{
		my @libraryList = @{$libraryListResult->{$_}};
		my $libraryDetails = decode_json $libraryList[8];
		$libraryDetails->{'comments'} = escapeHTML($libraryDetails->{'comments'});
		$libraryId .= "<option value='$libraryList[0]' title='$libraryDetails->{'comments'}'>Library: $libraryList[2]</option>";
	}
}
my $agpObjectComponent;
for (sort keys %{$objectComponent})
{
	$agpObjectComponent .= ($_ == 1) ? "<option value='$_' title='$objectComponent->{$_}' selected>$objectComponent->{$_}</option>" : "<option value='$_' title='$objectComponent->{$_}'>$objectComponent->{$_}</option>";
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$asbProjectId/$asbProjectId/g;
$html =~ s/\$forAssembly/$forAssembly/g;
$html =~ s/\$libraryId/$libraryId/g;
$html =~ s/\$agpObjectComponent/$agpObjectComponent/g;

print header;
print $html;

__DATA__
<form id="newGenome" name="newGenome" action="genomeSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	$asbProjectId
	<table>
	<tr><td style='text-align:right'><label for="newGenomeName"><b>Genome Name</b></label></td><td><input class='ui-widget-content ui-corner-all' name="name" id="newGenomeName" size="40" type="text" maxlength="32" /></td></tr>
	<tr><td style='text-align:right' rowspan='2'><label for="newGenomeFile"><b>Sequence File</b></label></td><td><sup class="ui-state-error-text">Note: Sequence Id should contain no more than 32 characters.</sup><br><input name="genomeFile" id="newGenomeFile" type="file" />(in FASTA format)</td></tr>
	<tr><td>or <input name="genomeFilePath" id="newGenomeFilePath" type="text" />(On-server file name with full path)</td></tr>
	<tr><td style='text-align:right'><label for="newGenomeLibraryId"><b>Link to</b></label></td><td><select class='ui-widget-content ui-corner-all' name='libraryId' id='newGenomeLibraryId'>$libraryId</select></td></tr>
	<tr><td style='text-align:right'><b>Import Options</b></td>
		<td>
			<hr>
			<input type="checkbox" id="newGenomeAssignChr" name="assignChr" value="1"><label for="newGenomeAssignChr">Assign chromosome number based on sequence name</label><br><sub>Manual assignment for sequences on unknown chromosomes is required after uploading.</sub><hr width="80%">
			<input type="checkbox" id="newGenomeSplit" name="split" value="1"><label for="newGenomeSplit">Split gapped sequences</label><hr width="80%">
			<input type="checkbox" id="newGenomeForAssembly" name="forAssembly" value="1"$forAssembly><label for="newGenomeForAssembly">Enable this genome to be reassembled with GPM</label><br>
			<sub><label for="newAgpFile"><b><select class='ui-widget-content ui-corner-all' name='agpObjectComponent' id='newAgpObjectComponent'>$agpObjectComponent</select> AGP file for guiding re-assembling</b></label>(Optional)</sub><input name="agpFile" id="newAgpFile" type="file" />(<a title="You may upload an AGP with a long file name, but only the first 32 characters will be saved.">Maximum 32 characters</a>)<hr width="80%">
			<input type="checkbox" id="newGenomeAsReference" name="asReference" value="1"><label for="newGenomeAsReference">Enable this genome to be used as a reference</label>
			<hr>
		</td>
	</tr>
	<tr><td style='text-align:right'><label for="newGenomeDescription"><b>Description</b></label></td><td><textarea class='ui-widget-content ui-corner-all' name="description" id="newGenomeDescription" cols="50" rows="10"></textarea></td></tr>
	</table>
</form>
<script>
$( "#newGenomeForAssembly" ).buttonset();
$( "#newGenomeAsReference" ).buttonset();
$( "#newGenomeSplit" ).buttonset();
$('#dialog').dialog("option", "title", "New Genome");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('newGenome'); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
</script>