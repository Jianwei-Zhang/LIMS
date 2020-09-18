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

my %datasetType = (
	0=>'Universal',
	1=>'Species',
	2=>'Picture'
	);

my $datasetId = param ('datasetId') || '';
my $dataset = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$dataset->execute($datasetId);
my @dataset=$dataset->fetchrow_array();

my $datasetStatus;
$datasetStatus->{0} = "not ";
$datasetStatus->{-1} = "is being ";
$datasetStatus->{1} = ($dataset[4] > 1) ? "$dataset[4] records " : "$dataset[4] record ";

my $datasetTypeList = '';
foreach (sort {$a <=> $b} keys %datasetType)
{
	$datasetTypeList .= ($_ == $dataset[3]) ? "<option value='$_' selected>$datasetType{$_}</option>" : "<option value='$_'>$datasetType{$_}</option>";
}
$html =~ s/\$datasetType/$datasetTypeList/g;

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
		$parentId .= ($parentList[0] eq $dataset[6] ) ?
			"<option value='$parentList[0]' title='$parentList[2]' selected>$parentList[1]: $parentList[2]</option>" :
			"<option value='$parentList[0]' title='$parentList[2]'>$parentList[1]: $parentList[2]</option>";
	}
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$datasetId/$datasetId/g;
$html =~ s/\$datasetName/$dataset[2]/g;
$html =~ s/\$parentId/$parentId/g;
$html =~ s/\$datasetStatus/$datasetStatus->{$dataset[7]}/g;
$html =~ s/\$datasetDescription/$dataset[8]/g;
$html =~ s/\$datasetCreator/$dataset[9]/g;
$html =~ s/\$datasetCreationDate/$dataset[10]/g;
print header;
print $html;

__DATA__
<form id="editDataset" name="editDataset" action="datasetSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="datasetId" id="editDatasetId" type="hidden" value="$datasetId" />
	<table>
	<tr><td style='text-align:right'><label for="editDatasetName"><b>Dataset Name</b></label></td><td><input class='ui-widget-content ui-corner-all' name="name" id="editDatasetName" size="40" type="text" maxlength="32" value="$datasetName"/><br><sup class='ui-state-disabled'>$datasetStatus loaded by $datasetCreator on $datasetCreationDate</sup></td></tr>
	<tr><td style='text-align:right'><label for="editDatasetType"><b>Type</b></label></td><td><select class='ui-widget-content ui-corner-all' name="datasetType" id="editDatasetType">$datasetType</select></td></tr>
	<tr><td style='text-align:right' rowspan='3'><label for="editDatasetFile"><b>File</b></label></td><td><input name="datasetFile" id="editDatasetFile" type="file" />(in Tab Delimited Text format)</td></tr>
	<tr><td>or <input name="datasetFilePath" id="editDatasetFilePath" type="text" />(On-server file name with full path)</td></tr>
	<tr>
		<td><div class="ui-state-error-text">The new file will replace the existing data set.</div>
		</td>
	</tr>
	<tr><td><label for="editDatasetIdColumn">Column for Name </label></td><td><input type="text" id="editDatasetIdColumn" name="idColumn" value="1"></td></tr>
	<tr><td style='text-align:right'><label for="editDatasetParentId"><b>Link to</b></label></td><td><select class='ui-widget-content ui-corner-all' name='parentId' id='editDatasetParentId'>$parentId</select></td></tr>
	<tr><td style='text-align:right'><label for="editDatasetDescription"><b>Description</b></label></td><td><textarea class='ui-widget-content ui-corner-all' name="description" id="editDatasetDescription" cols="50" rows="10">$datasetDescription</textarea></td></tr>
	</table>
</form>
<script>
$('#dialog').dialog("option", "title", "Edit Dataset");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('editDataset'); } }, { text: "Delete", click: function() { deleteItem($datasetId); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
</script>