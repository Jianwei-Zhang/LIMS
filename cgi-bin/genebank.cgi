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

my $genebankId = param ('genebankId') || '';
my $button;
my $activeDart = 0;
my $acitveDartDetector = 0;
my $activeDataset = 0;
my $acitveDatasetDetector = 0;
my $cookieDart = cookie('dart') || '';
my $cookieDataset = cookie('dataset') || '';
if ($genebankId)
{
	my $genebank=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$genebank->execute($genebankId);
	my @genebank = $genebank->fetchrow_array();
	my $genebankDetails = decode_json $genebank[8];

	my $parent=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$parent->execute($genebank[6]);
	my @parent = $parent->fetchrow_array();	

	my $genebankLine = '';
	for (sort {$a <=> $b} keys %$genebankDetails)
	{
		$genebankDetails->{$_}->{'field'} = '' unless ($genebankDetails->{$_}->{'field'});
		$genebankDetails->{$_}->{'value'} = '' unless ($genebankDetails->{$_}->{'value'});
		$genebankDetails->{$_}->{'value'} = escapeHTML($genebankDetails->{$_}->{'value'});
		$genebankDetails->{$_}->{'value'} =~ s/\n/<br>/g;
		$genebankLine .= ($genebankDetails->{$_}->{'value'} =~ /:\/\//) ? "<tr><td style='text-align:right;width:200px;'><b>$genebankDetails->{$_}->{'field'}</b></td>
						<td><a href='$genebankDetails->{$_}->{'value'}' target='_blank'>$genebankDetails->{$_}->{'value'}</a></td>
					</tr>":
					"<tr><td style='text-align:right;width:200px;'><b>$genebankDetails->{$_}->{'field'}</b></td>
						<td>$genebankDetails->{$_}->{'value'}</td>
					</tr>";
	}

	my $relatedDart = 0;
	my $dartList=$dbh->prepare("SELECT * FROM matrix WHERE z = ? AND container LIKE 'dart'");
	$dartList->execute($genebankId);
	while (my @dartList = $dartList->fetchrow_array())
	{
		$relatedDart = "<div id='dartList$genebankId$$'>" unless ($relatedDart);
		$relatedDart .= "<h3 title='$dartList[8]'>
			<a href='genebankDart.cgi?dartId=$dartList[0]'>$dartList[2] ($dartList[3] SNPs & $dartList[4] Genotypes ) <sup class='ui-state-disabled'>loaded by $dartList[9] on $dartList[10]</sup></a>
			</h3>
			<div><img src='$commoncfg->{HTDOCS}/css/images/loading.gif'>Loading...</div>";
		$activeDart = $acitveDartDetector if ($dartList[0] eq $cookieDart);
		$acitveDartDetector++;
	}
	$relatedDart .= "</div>" if($relatedDart);

	my $relatedDataset = 0;
	my $datasetList=$dbh->prepare("SELECT * FROM matrix WHERE z = ? AND container LIKE 'dataset'");
	$datasetList->execute($genebankId);
	while (my @datasetList = $datasetList->fetchrow_array())
	{
		$relatedDataset = "<div id='datasetList$genebankId$$'>" unless ($relatedDataset);
		$relatedDataset .= "<h3 title='$datasetList[8]'>
			<a href='table.cgi?type=record&parentId=$datasetList[0]&refresh=viewGenebankTabs$genebankId'>$datasetType{$datasetList[3]} Dataset: $datasetList[2] ($datasetList[4] records ) <sup class='ui-state-disabled'>loaded by $datasetList[9] on $datasetList[10]</sup></a>
			</h3>
			<div><img src='$commoncfg->{HTDOCS}/css/images/loading.gif'>Loading...</div>";
		$activeDataset = $acitveDatasetDetector if ($datasetList[0] eq $cookieDataset);
		$acitveDatasetDetector++;
	}
	$relatedDataset .= "</div>" if($relatedDataset);

	$button = "<ul id='genebankInfoMenu$genebankId$$' style='margin-top: .3em;width: 250px;'><li><a><span class='ui-icon ui-icon-triangle-1-e'></span><b>Genebank '$genebank[2]'</b></a>
				<ul style='z-index: 1000;white-space: nowrap;'>
				<li><a onclick='openDialog(\"itemEdit.cgi?itemId=$genebankId\")' title='Edit/Delete $genebank[2]'><span class='ui-icon ui-icon-pencil'></span>Edit/Delete</a></li>
				<li><a onclick='openDialog(\"datasetNew.cgi?parentId=$genebankId\")' title='Upload New Dataset'>New Dataset</a></li>
				<li><a onclick='openDialog(\"dartNew.cgi?parentId=$genebankId\")' title='Upload New DArTseq'>New DArTseq</a></li>
				</ul></li></ul>
			<div id='viewGenebankTabs$genebankId'>
				<ul>
				<li><a href='#viewGenebankTabs-1'>General</a></li>
				";
	$button .= ($relatedDataset) ? "<li><a href='#datasetList$genebankId$$'>Datasets</a></li>" : "";				
	$button .= ($relatedDart) ? "<li><a href='#dartList$genebankId$$'>DArTseq</a></li>" : "";
	$button .= "</ul>
				<div id='viewGenebankTabs-1'>
					<table>
					$genebankLine
					</table>
				</div>
				";
	$button .= ($relatedDataset) ? $relatedDataset : "";
	$button .= ($relatedDart) ? $relatedDart : "";
	$button .= "</div>";

	undef $/;# enable slurp mode
	my $html = <DATA>;
	$html =~ s/\$button/$button/g;
	$html =~ s/\$genebankId/$genebankId/g;
	$html =~ s/\$\$/$$/g;
	$html =~ s/\$activeDataset/$activeDataset/g;
	$html =~ s/\$activeDart/$activeDart/g;
	$html =~ s/\$commoncfg->{HTDOCS}/$commoncfg->{HTDOCS}/g;

	print header(-cookie=>[cookie(-name=>'genebank',-value=>$genebankId),cookie(-name=>'dart',-value=>0),cookie(-name=>'dataset',-value=>0)]);
	print $html;
}
else
{
	print header(-type=>'text/html',-status=>'402 Invalid operation');
	exit;
}

__DATA__
$button
<script>
buttonInit();
$( "#viewGenebankTabs$genebankId" ).tabs({
	// loading spinner
	beforeLoad: function(event, ui) {
		ui.panel.html('<img src="$commoncfg->{HTDOCS}/css/images/loading.gif" style="vertical-align:middle;"> Loading...');
	}
});
$( "#genebankInfoMenu$genebankId$$" ).menu();
$( "#dartList$genebankId$$" ).accordion({
        active:$activeDart,
        create:function(event, ui) {
            ui.panel.load(ui.header.find('a').attr('href'));
        },
        activate:function(event, ui) {
            ui.newPanel.load(ui.newHeader.find('a').attr('href'));
        },
        heightStyle: "content"
});
$( "#datasetList$genebankId$$" ).accordion({
        active:$activeDataset,
        create:function(event, ui) {
            ui.panel.load(ui.header.find('a').attr('href'));
        },
        activate:function(event, ui) {
            ui.newPanel.load(ui.newHeader.find('a').attr('href'));
        },
        heightStyle: "content"
});
loadingHide();
</script>