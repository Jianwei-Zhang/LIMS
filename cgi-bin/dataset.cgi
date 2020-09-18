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
if(!$userId)
{
	print header;
	print <<eof;
	<h2>Error</h2>
	Please <a onclick="openDialog('login.cgi');">login</a> first!
	<script type="text/javascript">
		openDialog('login.cgi');
		loadingHide();
	</script>
eof
    exit;
}

my $commoncfg = readConfig("main.conf");
my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});

my %datasetType = (
	0=>'Universal',
	1=>'Species',
	2=>'Picture'
	);

my $datasets = '';
my $allDataset=$dbh->prepare("SELECT * FROM matrix WHERE container Like 'dataset'");
$allDataset->execute();
while (my @allDataset = $allDataset->fetchrow_array())
{
	$allDataset[8] = escapeHTML($allDataset[8]);
	$allDataset[8] =~ s/\n/<br>/g;

	$datasets = "<form id='datasetList$$' name='datasetList$$'>
		<table id='datasets$$' class='display' style='width: 100%;'>
			<thead>
				<tr>
					<th style='text-align:left'><b>Dataset</b></th>
					<th style='text-align:left'><b>Type</b></th>
					<th style='text-align:left'><b>Records</b></th>
					<th style='text-align:left'><b>Creator</b></th>
					<th style='text-align:left'><b>Creation Date</b></th>
				</tr>
			</thead>
			<tbody>" unless($datasets);
	$datasets .= "<tr>
		<td title='Dataset'><a id='datasetId$allDataset[0]$$' onclick='openDialog(\"datasetView.cgi?datasetId=$allDataset[0]\")' title='View'>$allDataset[2]</a></td>
		<td title='Type'>$datasetType{$allDataset[3]}</td>
		<td title='Click to download dataset'><a href='download.cgi?datasetId=$allDataset[0]' target='hiddenFrame'>$allDataset[4]</a></td>
		<td title='Creator'>$allDataset[9]</td>
		<td title='Creation date'>$allDataset[10]</td>
		</tr>";
}
$datasets .= "</tbody></table></form>\n" if($datasets);
my $button = "<div class='ui-state-highlight ui-corner-all' style='padding: 0 .7em;'>
	<button style='float: right; margin-top: .3em; margin-right: .3em;' onclick='openDialog(\"datasetNew.cgi\")'>New Dataset</button>
	<button style='float: right; margin-top: .3em; margin-right: .3em;' onclick='refresh(\"general\")'>Refresh</button>
	";
$button .= "<h2>Datasets</h2>";

unless($datasets)
{
	$button .= "<p class='ui-state-error ui-corner-all' style='padding: .7em;'><span class='ui-icon ui-icon-alert' style='float: left; margin-right: .3em;'></span>
			<strong>No dataset, please upload one!</strong></p>";
}
$button .= "</div>\n";

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$button/$button/g;
$html =~ s/\$datasets/$datasets/g;
$html =~ s/\$\$/$$/g;

print header(-cookie=>cookie(-name=>'general',-value=>5));
print $html;

__DATA__
$button
$datasets
<script>
buttonInit();
$( "#datasets$$" ).dataTable({
	"scrollY": "600px",
	"scrollCollapse": true,
	"paging": false
});
loadingHide();
</script>