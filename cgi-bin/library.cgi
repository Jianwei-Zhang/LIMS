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

my $libraryId = param ('libraryId') || '';
my %yesOrNo = ( 0=>'No', 1=>'Yes' );
my %status = ( active=>'Active', inactive=>'Inactive', deleted=>'Deleted' );
my $button;
my $tableHeader;
my $plates = '';
my $nonOrderableTargets;
my $activeFpc = 0;
my $acitveFpcDetector = 0;
my $cookieFpc = cookie('fpc') || '';
my $activeGenome = 0;
my $acitveGenomeDetector = 0;
my $cookieGenome = cookie('genome') || '';
if ($libraryId)
{
	my $library=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$library->execute($libraryId);
	my @library = $library->fetchrow_array();
	my $libraryDetails = decode_json $library[8];

	$libraryDetails->{'comments'} = escapeHTML($libraryDetails->{'comments'});
	$libraryDetails->{'orderPageComments'} = escapeHTML($libraryDetails->{'orderPageComments'});
	$libraryDetails->{'comments'} =~ s/\n/<br>/g;
	$libraryDetails->{'orderPageComments'} =~ s/\n/<br>/g;
	$libraryDetails->{'commonName'} = ($libraryDetails->{'commonName'}) ? "($libraryDetails->{'commonName'})" : "";
	if($library[4])
	{
		my $checkVector = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$checkVector->execute($library[4]);
		my @checkVector=$checkVector->fetchrow_array();
		$library[4] = "<a onclick='openDialog(\"vectorView.cgi?vectorId=$library[4]\");'>$checkVector[2]</a>";
	}
	if($library[5])
	{
		my $checkSource = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$checkSource->execute($library[5]);
		my @checkSource=$checkSource->fetchrow_array();
		$library[5] = "<tr><td style='text-align:right'><b>Rearraying source</b></td>
						<td><a onclick='openDialog(\"libraryView.cgi?libraryId=$library[5]\");'>$checkSource[2]</a></td>
					</tr>";
	}
	else
	{
		$library[5] = "";
	}
	my $cloneList=$dbh->prepare("SELECT * FROM clones WHERE libraryId = ?");
	$cloneList->execute($libraryId);
	my $tagList=$dbh->prepare("SELECT * FROM matrix WHERE x = ? AND container LIKE 'tag'");
	$tagList->execute($libraryId);
	my $besList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 98 AND x = ?");
	$besList->execute($libraryId);
	my $relatedFpc = 0;
	my $fpcList=$dbh->prepare("SELECT * FROM matrix WHERE barcode = ? AND container LIKE 'fpc'");
	$fpcList->execute($libraryId);
	while (my @fpcList = $fpcList->fetchrow_array())
	{
		$relatedFpc = "<div id='fpcList$libraryId$$'>" unless ($relatedFpc);
		$relatedFpc .= "<h3 title='$fpcList[8]'>
			<a href='fpcList.cgi?fpcId=$fpcList[0]'>$fpcList[2] v.$fpcList[3]</a>
			</h3>
			<div><img src='$commoncfg->{HTDOCS}/css/images/loading.gif'>Loading...</div>";
		$activeFpc = $acitveFpcDetector if ($fpcList[0] eq $cookieFpc);
		$acitveFpcDetector++;
	}
	$relatedFpc .= "</div>" if($relatedFpc);

	my $relatedGenome = 0;
	my $genomeList=$dbh->prepare("SELECT * FROM matrix WHERE z = ? AND container LIKE 'genome'");
	$genomeList->execute($libraryId);
	while (my @genomeList = $genomeList->fetchrow_array())
	{
		$relatedGenome = "<div id='genomeList$libraryId$$'>" unless ($relatedGenome);
		$relatedGenome .= "<h3 title='$genomeList[8]'>
			<a href='libraryGenome.cgi?genomeId=$genomeList[0]'>$genomeList[2] ($genomeList[3] Sequences) <sup class='ui-state-disabled'>loaded by $genomeList[9] on $genomeList[10]</sup></a>
			</h3>
			<div><img src='$commoncfg->{HTDOCS}/css/images/loading.gif'>Loading...</div>";
		$activeGenome = $acitveGenomeDetector if ($genomeList[0] eq $cookieGenome);
		$acitveGenomeDetector++;
	}
	$relatedGenome .= "</div>" if($relatedGenome);

	$button = "<ul id='libraryInfoMenu$libraryId$$' style='margin-top: .3em;width: 250px;'><li><a><span class='ui-icon ui-icon-triangle-1-e'></span><b>Library '$library[2]'</b></a>
				<ul style='z-index: 1000;white-space: nowrap;'>
				<li><a onclick='openDialog(\"libraryEdit.cgi?libraryId=$libraryId\")' title='Edit/Delete $library[2]'><span class='ui-icon ui-icon-pencil'></span>Edit/Delete</a></li>
				<li><a><span class='ui-icon ui-icon-document'></span>Clone List</a><ul style='z-index: 1000;white-space: nowrap;'>";
	$button .= ($library[3] > 0) ? "<li><a onclick='openDialog(\"cloneNew.cgi?libraryId=$libraryId\")'><span class='ui-icon ui-icon-bullet'></span>Generate</a></li>" : "<li class='ui-state-disabled'><a><span class='ui-icon ui-icon-bullet'></span>Generate</a></li>";
	$button .= ($cloneList->rows > 0) ? "<li><a href='download.cgi?cloneLibraryId=$libraryId' target='hiddenFrame'><span class='ui-icon ui-icon-bullet'></span>Download</a></li></ul></li>" : "<li class='ui-state-disabled'><a><span class='ui-icon ui-icon-bullet'></span>Download</a></li></ul></li>";
	$button .= ($tagList->rows > 0) ? "<li><a><span class='ui-icon ui-icon-tag'></span>WGP Tags</a><ul style='z-index: 1000;'><li><a href='download.cgi?tagLibraryId=$libraryId' target='hiddenFrame'><span class='ui-icon ui-icon-bullet'></span>Download</a></li><li><a onclick='openDialog(\"tagNew.cgi?libraryId=$libraryId\")'><span class='ui-icon ui-icon-bullet'></span>Load New</a></li></ul></li>" : "<li><a onclick='openDialog(\"tagNew.cgi?libraryId=$libraryId\")'><span class='ui-icon ui-icon-tag'></span>Load WGP Tags</a></li>";
	$button .= ($besList->rows > 0) ? "<li><a><span class='ui-icon ui-icon-seek-end'></span>BAC End Seqs</a><ul style='z-index: 1000;'><li><a onclick='openDialog(\"besNew.cgi?libraryId=$libraryId\")'><span class='ui-icon ui-icon-bullet'></span>Load New BES</a></li><li><a href='download.cgi?besLibraryId=$libraryId' target='hiddenFrame'><span class='ui-icon ui-icon-bullet'></span>Download</a></li><li><a onclick='openDialog(\"besToSeqForm.cgi?libraryId=$libraryId&targetId=\")'><span class='ui-icon ui-icon-bullet'></span>BES to Seq</a></li><li><a onclick='openDialog(\"besReport.cgi?libraryId=$libraryId\")'><span class='ui-icon ui-icon-bullet'></span>Report</a></li></ul></li>" : "<li><a onclick='openDialog(\"besNew.cgi?libraryId=$libraryId\")'><span class='ui-icon ui-icon-seek-end'></span>Load BAC End Seqs</a></li>";
	$button .= "<li><a onclick='openDialog(\"fpcNew.cgi?libraryId=$libraryId\")'><span class='ui-icon ui-icon-document-b'></span>Load FPC</a></li>";
	$button .= ($library[7]) ? "<li><a onclick='openDialog(\"seqLibraryReport.cgi?libraryId=$libraryId\")'><span class='ui-icon ui-icon-note'></span>Sequencing Report</a></li>
					<li><a href='download.cgi?libraryId=$libraryId' target='hiddenFrame'><span class='ui-icon ui-icon-disk'></span>Download Sequences</a></li>" : "";
	$button .= "</ul></li></ul>
			<div id='viewLibraryTabs$libraryId'>
				<ul>
				<li><a href='#viewLibraryTabs-1'>General</a></li>
				<li><a href='#viewLibraryTabs-2'>Specs</a></li>
				<li><a href='#viewLibraryTabs-3'>Status</a></li>
				<li><a href='libraryPlate.cgi?libraryId=$libraryId'>Plates</a></li>
				";
				
	$button .= ($relatedFpc) ? "<li><a href='#fpcList$libraryId$$'>FPC</a></li>" : "";
	$button .= ($relatedGenome) ? "<li><a href='#genomeList$libraryId$$'>Genome</a></li>" : "";
	$button .= ($library[7]) ? "<li><a href='seqLibrary.cgi?libraryId=$libraryId'>Sequencing</a></li>" : "";
	$button .= "</ul>
				<div id='viewLibraryTabs-1'>
					<table>
					<tr><td style='text-align:right;width:200px;'><b>Library Name</b></label>
						<td>$library[2]<sup class='ui-state-disabled'>(Library id: $libraryId)</sup>
						<br><sup class='ui-state-disabled'>Last changed by $library[9] on $library[10]</sup></td>
					</tr>
					<tr><td style='text-align:right'><b>Nickname</b></td>
						<td>$libraryDetails->{'nickname'}</td>
					</tr>
					<tr><td style='text-align:right'><b>Library type</b></td>
						<td>$libraryDetails->{'type'}</td>
					</tr>
					$library[5]
					<tr><td style='text-align:right'><b>Plate format</b></td>
						<td>$libraryDetails->{'format'}</td>
					</tr>
					<tr><td style='text-align:right'><b>Distributor institute</b></td>
						<td>$libraryDetails->{'distributorInstitute'}</td>
					</tr>
					<tr><td style='text-align:right'><b>Ref. to publication</b></td>
						<td>$libraryDetails->{'refToPublication'}</td>
					</tr>
					<tr><td style='text-align:right'><b>Provided by</b></td>
						<td>$libraryDetails->{'providedBy'}</td>
					</tr>
					<tr><td style='text-align:right'><b>Organization</b></td>
						<td>$libraryDetails->{'organization'}</td>
					</tr>
					<tr><td style='text-align:right;width:200px;'><b>Genus Species Subspecies</b><br>(Common Name)</td>
						<td>$libraryDetails->{'genus'} $libraryDetails->{'species'} $libraryDetails->{'subspecies'}<br>$libraryDetails->{'commonName'}</td>
					</tr>
					<tr><td style='text-align:right'><b>Cultivar</b></td>
						<td>$libraryDetails->{'cultivar'}</td>
					</tr>
					<tr><td style='text-align:right'><b>Genome equivalents</b></td>
						<td>$libraryDetails->{'genomeEquivalents'}</td>
					</tr>
					<tr><td style='text-align:right'><b>Genome size</b></td>
						<td>$libraryDetails->{'genomeSize'}</td>
					</tr>
					<tr><td style='text-align:right'><b>Genome type</b></td>
						<td>$libraryDetails->{'genomeType'}</td>
					</tr>
					<tr><td style='text-align:right'><b>Tissue</b></td>
						<td>$libraryDetails->{'tissue'}</td>
					</tr>
					<tr><td style='text-align:right'><b>Treatment</b></td>
						<td>$libraryDetails->{'treatment'}</td>
					</tr>
					<tr><td style='text-align:right'><b>Development stage</b></td>
						<td>$libraryDetails->{'developmentStage'}</td>
					</tr>
					</table>		
				</div>
				<div id='viewLibraryTabs-2'>
					<table>
					<tr><td style='text-align:right;width:200px;'><b># of plates</b></td>
						<td>$library[3]</td>
					</tr>
					<tr><td style='text-align:right'><b># of filters</b></td>
						<td>$libraryDetails->{'nofFilters'}</td>
					</tr>
					<tr><td style='text-align:right'><b>Average insert size</b></td>
						<td>$libraryDetails->{'averageInsertSize'}</td>
					</tr>
					<tr><td style='text-align:right'><b>Standard deviation</b></td>
						<td>$libraryDetails->{'standardDeviation'}</td>
					</tr>
					<tr><td style='text-align:right'><b>Vector</b></td>
						<td>$library[4]</td>
					</tr>
					<tr><td style='text-align:right'><b>Host</b></td>
						<td>$libraryDetails->{'host'}</td>
					</tr>
					<tr><td style='text-align:right'><b>Enzyme 5' prime</b></td>
						<td>$libraryDetails->{'enzymeFivePrime'}</td>
					</tr>
					<tr><td style='text-align:right'><b>Enzyme 3' prime</b></td>
						<td>$libraryDetails->{'enzymeThreePrime'}</td>
					</tr>
					<tr><td style='text-align:right'><b>End seq primer 5' prime</b></td>
						<td>$libraryDetails->{'endSeqFivePrime'}</td>
					</tr>
					<tr><td style='text-align:right'><b>End seq primer 3' prime</b></td>
						<td>$libraryDetails->{'endSeqThreePrime'}</td>
					</tr>
					<tr><td style='text-align:right'><b>Cloning linker 5' prime</b></td>
						<td>$libraryDetails->{'cloningFivePrime'}</td>
					</tr>
					<tr><td style='text-align:right'><b>Cloning linker 3' prime</b></td>
						<td>$libraryDetails->{'cloningThreePrime'}</td>
					</tr>
					<tr><td style='text-align:right'><b>Antibiotic</b></td>
						<td>$libraryDetails->{'antibiotic'}</td>
					</tr>
					</table>
				</div>
				<div id='viewLibraryTabs-3'>
					<table>
					<tr><td style='text-align:right;width:200px;'><b>Status</b></td>
						<td>$status{$libraryDetails->{'status'}}</td>
					</tr>
					<tr><td style='text-align:right'><b>Is external?</b></td>
						<td>$yesOrNo{$libraryDetails->{'isExternal'}}</td>
					</tr>
					<tr><td style='text-align:right'><b>Is public?</b></td>
						<td>$yesOrNo{$libraryDetails->{'isPublic'}}</td>
					</tr>
					<tr><td style='text-align:right'><b>Is finished?</b></td>
						<td>$yesOrNo{$libraryDetails->{'isFinished'}}</td>
					</tr>
					<tr><td style='text-align:right'><b>Is library for sale?</b></td>
						<td>$yesOrNo{$libraryDetails->{'isLibraryForSale'}}</td>
					</tr>
					<tr><td style='text-align:right'><b>Is filter for sale?</b></td>
						<td>$yesOrNo{$libraryDetails->{'isFilterForSale'}}</td>
					</tr>
					<tr><td style='text-align:right'><b>Is clone for sale?</b></td>
						<td>$yesOrNo{$libraryDetails->{'isCloneForSale'}}</td>
					</tr>
					<tr><td style='text-align:right'><b>Is library for sequencing?</b></td>
						<td>$yesOrNo{$library[7]}</td>
					</tr>
					<tr><td style='text-align:right'><b>Order page comments</b></td>
						<td>$libraryDetails->{'orderPageComments'}</td>
					</tr>
					<tr><td style='text-align:right;width:200px;'><b>Library was made on</b></td>
						<td>$libraryDetails->{'dateLibraryWasMade'}</td>
					</tr>
					<tr><td style='text-align:right'><b>Library was autoclaved on</b></td>
						<td>$libraryDetails->{'dateLibraryWasAutoclaved'}</td>
					</tr>
					<tr><td style='text-align:right'><b>Date entered the system</b></td>
						<td>$library[10]</td>
					</tr>
					<tr><td style='text-align:right'><b>Archive start date</b></td>
						<td>$libraryDetails->{'archiveStartDate'}</td>
					</tr>
					<tr><td style='text-align:right'><b>Archive end date</b></td>
						<td>$libraryDetails->{'archiveEndDate'}</td>
					</tr>
					<tr><td style='text-align:right'><b>Comments</b></td>
						<td>$libraryDetails->{'comments'}</td>
					</tr>
					</table>
				</div>
				";
	$button .= ($relatedFpc) ? $relatedFpc : "";
	$button .= ($relatedGenome) ? $relatedGenome : "";
	$button .= "</div>";
	
	undef $/;# enable slurp mode
	my $html = <DATA>;
	$html =~ s/\$button/$button/g;
	$html =~ s/\$libraryId/$libraryId/g;
	$html =~ s/\$\$/$$/g;
	$html =~ s/\$nonOrderableTargets/$nonOrderableTargets/g;
	$html =~ s/\$activeFpc/$activeFpc/g;
	$html =~ s/\$activeGenome/$activeGenome/g;
	$html =~ s/\$commoncfg->{HTDOCS}/$commoncfg->{HTDOCS}/g;

	print header(-cookie=>[cookie(-name=>'library',-value=>$libraryId),cookie(-name=>'sample',-value=>0)]);
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
$( "#viewLibraryTabs$libraryId" ).tabs({
	// loading spinner
	beforeLoad: function(event, ui) {
		ui.panel.html('<img src="$commoncfg->{HTDOCS}/css/images/loading.gif" style="vertical-align:middle;"> Loading...');
	}
});
$( "#libraryInfoMenu$libraryId$$" ).menu();
$( "#fpcList$libraryId$$" ).accordion({
        active:$activeFpc,
        create:function(event, ui) {
            ui.panel.load(ui.header.find('a').attr('href'));
        },
        activate:function(event, ui) {
            ui.newPanel.load(ui.newHeader.find('a').attr('href'));
        },
        heightStyle: "content"
});
$( "#genomeList$libraryId$$" ).accordion({
        active:$activeGenome,
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
