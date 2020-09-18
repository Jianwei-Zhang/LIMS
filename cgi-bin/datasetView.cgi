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
my @dataset = $dataset->fetchrow_array();
my $datasetStatus;
$datasetStatus->{0} = "not ";
$datasetStatus->{-1} = "is being ";
$datasetStatus->{1} = ($dataset[4] > 1) ? "$dataset[4] records " : "$dataset[4] record ";

my $relatedParent = 'None.';
if ($dataset[6] > 0)
{
	my $relatedParentForDataset=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$relatedParentForDataset->execute($dataset[6]);
	my @relatedParentForDataset = $relatedParentForDataset->fetchrow_array();
	$relatedParent = "<a onclick='closeDialog();openDialog(\"itemView.cgi?itemId=$dataset[6]\")'>$relatedParentForDataset[2]</a> ";
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$datasetId/$datasetId/g;
$html =~ s/\$datasetName/$dataset[2]/g;
$html =~ s/\$datasetType/$datasetType{$dataset[3]}/g;
$html =~ s/\$relatedParent/$relatedParent/g;
$html =~ s/\$datasetStatus/$datasetStatus->{$dataset[7]}/g;
$dataset[8] = escapeHTML($dataset[8]);
$dataset[8] =~ s/\n/<br>/g;
$html =~ s/\$datasetDescription/$dataset[8]/g;
$html =~ s/\$datasetCreator/$dataset[9]/g;
$html =~ s/\$datasetCreationDate/$dataset[10]/g;
$html =~ s/\$\$/$$/g;
print header;
print $html;

__DATA__
<table>
	<tr><td style='text-align:right;white-space: nowrap;'><b>$datasetType Dataset</b></td><td>$datasetName<br>$datasetStatus <sup class='ui-state-disabled'>loaded by $datasetCreator on $datasetCreationDate</sup></td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'>Link to:</td><td>$relatedParent</td></tr>
	<tr><td style='text-align:right'><b>Description</b></td><td>$datasetDescription</td></tr>
</table>
<script>
$('#dialog').dialog("option", "title", "View Dataset");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Edit", click: function() { closeDialog();openDialog("datasetEdit.cgi?datasetId=$datasetId"); } }, { text: "OK", click: function() {closeDialog(); } } ] );
</script>