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
my @dart = $dart->fetchrow_array();
my $dartStatus;
$dartStatus->{0} = "not ";
$dartStatus->{-1} = "is being ";
$dartStatus->{1} = ($dart[3] > 1) ? ($dart[4] > 1) ? "$dart[3] SNPs and $dart[4] genotypes " : "$dart[3] SNPs and $dart[4] genotype " : ($dart[4] > 1) ? "$dart[3] SNP and $dart[4] genotypes " :"$dart[3] SNP and $dart[4] genotype ";

my $relatedGenebank = 'None.';
if ($dart[6] > 0)
{
	my $relatedGenebankForDart=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$relatedGenebankForDart->execute($dart[6]);
	my @relatedGenebankForDart = $relatedGenebankForDart->fetchrow_array();
	$relatedGenebank = "<a onclick='closeDialog();openDialog(\"itemView.cgi?itemId=$dart[6]\")'>$relatedGenebankForDart[2]</a> ";
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$dartId/$dartId/g;
$html =~ s/\$dartName/$dart[2]/g;
$html =~ s/\$relatedGenebank/$relatedGenebank/g;
$html =~ s/\$dartStatus/$dartStatus->{$dart[7]}/g;
$dart[8] = escapeHTML($dart[8]);
$dart[8] =~ s/\n/<br>/g;
$html =~ s/\$dartDescription/$dart[8]/g;
$html =~ s/\$dartCreator/$dart[9]/g;
$html =~ s/\$dartCreationDate/$dart[10]/g;
$html =~ s/\$\$/$$/g;
print header;
print $html;

__DATA__
<table>
	<tr><td style='text-align:right;white-space: nowrap;'><b>DArTseq</b></td><td>$dartName<br>$dartStatus <sup class='ui-state-disabled'>loaded by $dartCreator on $dartCreationDate</sup></td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'>Link to genebank:</td><td>$relatedGenebank</td></tr>
	<tr><td style='text-align:right'><b>Description</b></td><td>$dartDescription</td></tr>
</table>
<script>
$('#dialog').dialog("option", "title", "View DArTseq");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Edit", click: function() { closeDialog();openDialog("dartEdit.cgi?dartId=$dartId"); } }, { text: "OK", click: function() {closeDialog(); } } ] );
</script>