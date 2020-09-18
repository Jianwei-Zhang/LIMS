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
my $objectComponent;
$objectComponent->{0} = "Unknown";
$objectComponent->{1} = "Chr-Seq";
$objectComponent->{2} = "Ctg-Seq";

my $genomes = '';
my %yesOrNo = ( 0=>'No', 1=>'Yes' );
my $allGenome=$dbh->prepare("SELECT * FROM matrix WHERE container Like 'genome'");
$allGenome->execute();
if($allGenome->rows > 0)
{
# 	my $allLibraryResult;
# 	my $allLibrary=$dbh->prepare("SELECT * FROM matrix WHERE container Like 'library'");
# 	$allLibrary->execute();
# 	while (my @allLibrary = $allLibrary->fetchrow_array())
# 	{
# 		@{$allLibraryResult->{$allLibrary[0]}} = @allLibrary;
# 	}
	my $allGenomeResult;
	while (my @allGenome = $allGenome->fetchrow_array())
	{
		@{$allGenomeResult->{$allGenome[2]}} = @allGenome;
	}
	foreach (sort {uc ($a) cmp uc($b)} keys %$allGenomeResult)
	{
		my @allGenome = @{$allGenomeResult->{$_}};
		$allGenome[8] = escapeHTML($allGenome[8]);
		$allGenome[8] =~ s/\n/<br>/g;
#		my $relatedLibraries = ($allGenome[6]) ? "<a onclick='openDialog(\"libraryView.cgi?libraryId=$allGenome[6]\")'>${$allLibraryResult->{$allGenome[6]}}[2]</a> " :'';
		my $relatedLibraries = ($allGenome[6]) ? "<a onclick='openDialog(\"libraryView.cgi?libraryId=$allGenome[6]\")'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-link' title='View'></span>View</a> " :'';
	
		my $agpAvailable = "";
		my $agpList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'agp' AND x = ?");# ORDER BY o
		$agpList->execute($allGenome[0]);
		while (my @agpList = $agpList->fetchrow_array())
		{
			$agpAvailable .= ($agpAvailable) ? "<br> <a href='download.cgi?agpId=$agpList[0]' target='hiddenFrame' title='$objectComponent->{$agpList[5]} AGP v.$agpList[3]'>$agpList[2]</a>" : " Available AGP:<br> <a href='download.cgi?agpId=$agpList[0]' target='hiddenFrame' title='$objectComponent->{$agpList[5]} AGP v.$agpList[3]'>$agpList[2]</a>";
		}

		$genomes = "<form id='genomeList$$' name='genomeList$$'>
			<table id='genomes$$' class='display' style='width: 100%;'>
				<thead>
					<tr>
						<th>
							<input type='checkbox' id='checkAllBox$$' name='checkAllBox$$' value='Check all' checked='checked' onClick='checkAll(\"itemId\");return false;' title='Check all'>
							<input type='checkbox' id='uncheckAllBox$$' name='uncheckAllBox$$' value='Uncheck all' onClick='uncheckAll(\"itemId\");return false;' title='Uncheck all'>
						</th>
						<th style='text-align:left'><b>Genome</b></th>
						<th style='text-align:left'><b>Sequences</b></th>
						<th style='text-align:left'><b>For Reassembly</b></th>
						<th style='text-align:left'><b>As Reference</b></th>
						<th style='text-align:left'><b>Related Library</b></th>
						<th style='text-align:left'><b>Related Assemblies</b></th>
						<th style='text-align:left'><b>Creator</b></th>
						<th style='text-align:left'><b>Creation Date</b></th>
					</tr>
				</thead>
				<tbody>" unless($genomes);
		$genomes .= "<tr onclick='checkOrUncheck(\"genomeList$allGenome[0]$$\");'>
			<td style='text-align:center;'><input type='checkbox' id='genomeList$allGenome[0]$$' name='itemId' value='$allGenome[0]' onClick='checkOrUncheck(\"genomeList$allGenome[0]$$\");'></td>
			<td title='Genome'><a id='genomeId$allGenome[0]$$' onclick='openDialog(\"genomeView.cgi?genomeId=$allGenome[0]\")' title='View'>$allGenome[2]</a></td>
			<td title='Click to download $allGenome[3] Sequences'><a href='download.cgi?genomeId=$allGenome[0]' target='hiddenFrame'>$allGenome[3]<span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-disk'></span></a></td>
			<td>$yesOrNo{$allGenome[4]}. $agpAvailable</td>
			<td title='As a reference for assembly'>$yesOrNo{$allGenome[5]}</td>
			<td title='Related Library'>$relatedLibraries</td>
			<td><div id='relatedAssembly$allGenome[0]' style='position: relative;'><script>loaddiv(\"relatedAssembly$allGenome[0]\",\"genomeToAssembly.cgi?genomeId=$allGenome[0]\");</script></div></td>
			<td title='Creator'>$allGenome[9]</td>
			<td title='Creation date'>$allGenome[10]</td>
			</tr>";
	}
}
$genomes .= "</tbody></table></form>\n" if($genomes);

my $button = "<div class='ui-state-highlight ui-corner-all' style='padding: 0 .7em;'>
	<button style='float: right; margin-top: .3em; margin-right: .3em;' onclick='openDialog(\"alignmentLoadForm.cgi\")'>Load Tabular Alignment</button>
	<button style='float: right; margin-top: .3em; margin-right: .3em;' onclick='openDialog(\"alignmentForm.cgi\")'>Run Alignment</button>";
$button .= "<button style='float: right; margin-top: .3em; margin-right: .3em;' onclick='openDialogForm(\"itemDeleteForm.cgi\",\"genomeList$$\")'>Delete</button>" if($allGenome->rows > 0);
$button .= "<button style='float: right; margin-top: .3em; margin-right: .3em;' onclick='openDialogForm(\"genomeMergeForm.cgi\",\"genomeList$$\")'>Merge</button>" if($allGenome->rows > 1);
$button .= "<button style='float: right; margin-top: .3em; margin-right: .3em;' onclick='openDialogForm(\"genomeSplitForm.cgi\",\"genomeList$$\")'>Split</button>
	<button style='float: right; margin-top: .3em; margin-right: .3em;' onclick='openDialog(\"genomeNew.cgi\")'>New Genome</button>
	<button style='float: right; margin-top: .3em; margin-right: .3em;' onclick='refresh(\"general\")'>Refresh</button>
	<input style='float: right;; margin-top: .3em; margin-right: .3em;' class='ui-widget-content ui-corner-all' name='seqName' id='searchSeqName$$' size='16' type='text' maxlength='32' VALUE='' placeholder='Search Seq' />
	";
$button .= "<h2>Genomes</h2>";

unless($genomes)
{
	$button .= "<p class='ui-state-error ui-corner-all' style='padding: .7em;'><span class='ui-icon ui-icon-alert' style='float: left; margin-right: .3em;'></span>
			<strong>No genome, please upload one!</strong></p>";
}
$button .= "</div>\n";

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$button/$button/g;
$html =~ s/\$genomes/$genomes/g;
$html =~ s/\$\$/$$/g;

print header(-cookie=>cookie(-name=>'general',-value=>4));
print $html;

__DATA__
$button
$genomes
<script>
buttonInit();
$( "#genomes$$" ).dataTable({
	"scrollY": "600px",
	"scrollCollapse": true,
	"paging": false,
	"columnDefs": [
    { "orderable": false, "targets": 0 }
  ]
});
$( "#searchSeqName$$" ).autocomplete({
	source: "autoSeqSearch.cgi",
	minLength: 4,
	select: function( event, ui ) {
		openDialog("seqView.cgi?seqId=" + ui.item.id);
	}
});
loadingHide();
</script>