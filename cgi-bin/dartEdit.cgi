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

my $dartId = param ('dartId') || '';
my $dart = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$dart->execute($dartId);
my @dart=$dart->fetchrow_array();

my $dartStatus;
$dartStatus->{0} = "not ";
$dartStatus->{-1} = "is being ";
$dartStatus->{1} = ($dart[3] > 1) ? ($dart[4] > 1) ? "$dart[3] SNPs and $dart[4] genotypes " : "$dart[3] SNPs and $dart[4] genotype " : ($dart[4] > 1) ? "$dart[3] SNP and $dart[4] genotypes " :"$dart[3] SNP and $dart[4] genotype ";

my $genebankId = "<option class='ui-state-error-text' value='0'>None</option>";
my $genebankList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'genebank'");# ORDER BY name
$genebankList->execute();
if ($genebankList->rows > 0)
{
	my $genebankListResult;
	while (my @genebankList = $genebankList->fetchrow_array())
	{
		@{$genebankListResult->{$genebankList[2]}} = @genebankList;
	}
	foreach (sort {uc ($a) cmp uc($b)} keys %$genebankListResult)
	{
		my @genebankList = @{$genebankListResult->{$_}};
		my $genebankDetails = decode_json $genebankList[8];
		$genebankDetails->{'comments'} = escapeHTML($genebankDetails->{'comments'});
		$genebankId .= ($genebankList[0] eq $dart[6] ) ?
			"<option value='$genebankList[0]' title='$genebankDetails->{'comments'}' selected>Genebank: $genebankList[2]</option>" :
			"<option value='$genebankList[0]' title='$genebankDetails->{'comments'}'>Genebank: $genebankList[2]</option>";
	}
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$dartId/$dartId/g;
$html =~ s/\$dartName/$dart[2]/g;
$html =~ s/\$genebankId/$genebankId/g;
$html =~ s/\$dartStatus/$dartStatus->{$dart[7]}/g;
$html =~ s/\$dartDescription/$dart[8]/g;
$html =~ s/\$dartCreator/$dart[9]/g;
$html =~ s/\$dartCreationDate/$dart[10]/g;
print header;
print $html;

__DATA__
<form id="editDart" name="editDart" action="dartSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="dartId" id="editDartId" type="hidden" value="$dartId" />
	<table>
	<tr><td style='text-align:right'><label for="editDartName"><b>DArTseq Name</b></label></td><td><input class='ui-widget-content ui-corner-all' name="name" id="editDartName" size="40" type="text" maxlength="32" value="$dartName"/><br><sup class='ui-state-disabled'>$dartStatus loaded by $dartCreator on $dartCreationDate</sup></td></tr>
	<tr><td style='text-align:right' rowspan='3'><label for="editDartFile"><b>DArTseq File</b></label></td><td><input name="dartFile" id="editDartFile" type="file" />(in Tab Delimited Text format)</td></tr>
	<tr><td>or <input name="dartFilePath" id="editDartFilePath" type="text" />(On-server file name with full path)</td></tr>
	<tr>
		<td><div class="ui-state-error-text">The new file will replace the existing data set.</div>
		</td>
	</tr>
	<tr><td style='text-align:right'><label for="editDartGenebankId"><b>Link to</b></label></td><td><select class='ui-widget-content ui-corner-all' name='genebankId' id='editDartGenebankId'>$genebankId</select></td></tr>
	<tr><td style='text-align:right'><label for="editDartDescription"><b>Description</b></label></td><td><textarea class='ui-widget-content ui-corner-all' name="description" id="editDartDescription" cols="50" rows="10">$dartDescription</textarea></td></tr>
	</table>
</form>
<script>
$('#dialog').dialog("option", "title", "Edit DArTseq");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('editDart'); } }, { text: "Delete", click: function() { deleteItem($dartId); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
</script>