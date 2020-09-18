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
my $projectId = '';
my $library = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$library->execute($libraryId);
my @library=$library->fetchrow_array();
my $libraryDetails = decode_json $library[8];

my $libraryToProject=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$libraryToProject->execute($library[6]);
my @libraryToProject = $libraryToProject->fetchrow_array();

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$libraryId/$libraryId/g;
$html =~ s/\$libraryName/$library[2]/g;
$html =~ s/\$libraryNickname/$libraryDetails->{'nickname'}/g;
$html =~ s/\$libraryType/$libraryDetails->{'type'}/g;
if($libraryDetails->{'isExternal'})
{
	$html =~ s/\$editLibraryIsExternalRadio2/checked="checked"/g;
	$html =~ s/\$editLibraryIsExternalRadio1//g;
}
else
{
	$html =~ s/\$editLibraryIsExternalRadio2//g;
	$html =~ s/\$editLibraryIsExternalRadio1/checked="checked"/g;
}
$html =~ s/\$libraryNofPlates/$library[3]/g;
$html =~ s/\$libraryNofFilters/$libraryDetails->{'nofFilters'}/g;
$html =~ s/\$libraryGenomeEquivalents/$libraryDetails->{'genomeEquivalents'}/g;
$html =~ s/\$libraryGenomeSize/$libraryDetails->{'genomeSize'}/g;
$html =~ s/\$libraryAverageInsertSize/$libraryDetails->{'averageInsertSize'}/g;
$html =~ s/\$libraryStandardDeviation/$libraryDetails->{'standardDeviation'}/g;
if($libraryDetails->{'format'} == 96)
{
	$html =~ s/\$editLibraryFormatRadio1/checked="checked"/g;
	$html =~ s/\$editLibraryFormatRadio2//g;
}
else
{
	$html =~ s/\$editLibraryFormatRadio1//g;
	$html =~ s/\$editLibraryFormatRadio2/checked="checked"/g;
}
if($libraryDetails->{'status'} eq 'inactive')
{
	$html =~ s/\$editLibraryStatusRadio1//g;
	$html =~ s/\$editLibraryStatusRadio2/checked="checked"/g;
	$html =~ s/\$editLibraryStatusRadio3//g;
}
elsif($libraryDetails->{'status'} eq 'deleted')
{
	$html =~ s/\$editLibraryStatusRadio1//g;
	$html =~ s/\$editLibraryStatusRadio2//g;
	$html =~ s/\$editLibraryStatusRadio3/checked="checked"/g;
}
else
{
	$html =~ s/\$editLibraryStatusRadio1/checked="checked"/g;
	$html =~ s/\$editLibraryStatusRadio2//g;
	$html =~ s/\$editLibraryStatusRadio3//g;
}
$html =~ s/\$libraryDistributorInstitute/$libraryDetails->{'distributorInstitute'}/g;
if($library[4])
{
	my $checkVector = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$checkVector->execute($library[4]);
	my @checkVector=$checkVector->fetchrow_array();
	$library[4] = $checkVector[2];
}
$html =~ s/\$libraryVector/$library[4]/g;
$html =~ s/\$libraryHost/$libraryDetails->{'host'}/g;
$html =~ s/\$libraryEnzymeFivePrime/$libraryDetails->{'enzymeFivePrime'}/g;
$html =~ s/\$libraryEnzymeThreePrime/$libraryDetails->{'enzymeThreePrime'}/g;
$html =~ s/\$libraryEndSeqFivePrime/$libraryDetails->{'endSeqFivePrime'}/g;
$html =~ s/\$libraryEndSeqThreePrime/$libraryDetails->{'endSeqThreePrime'}/g;
$html =~ s/\$libraryCloningFivePrime/$libraryDetails->{'cloningFivePrime'}/g;
$html =~ s/\$libraryCloningThreePrime/$libraryDetails->{'cloningThreePrime'}/g;
$html =~ s/\$libraryAntibiotic/$libraryDetails->{'antibiotic'}/g;
$html =~ s/\$libraryGenus/$libraryDetails->{'genus'}/g;
$html =~ s/\$librarySpecies/$libraryDetails->{'species'}/g;
$html =~ s/\$librarySubspecies/$libraryDetails->{'subspecies'}/g;
$html =~ s/\$libraryCommonName/$libraryDetails->{'commonName'}/g;
$html =~ s/\$libraryGenomeType/$libraryDetails->{'genomeType'}/g;
$html =~ s/\$libraryTissue/$libraryDetails->{'tissue'}/g;
$html =~ s/\$libraryTreatment/$libraryDetails->{'treatment'}/g;
$html =~ s/\$libraryDevelopmentStage/$libraryDetails->{'developmentStage'}/g;
$html =~ s/\$libraryCultivar/$libraryDetails->{'cultivar'}/g;
$html =~ s/\$libraryRefToPublication/$libraryDetails->{'refToPublication'}/g;
if($libraryDetails->{'isPublic'})
{
	$html =~ s/\$editLibraryIsPublicRadio2/checked="checked"/g;
	$html =~ s/\$editLibraryIsPublicRadio1//g;
}
else
{
	$html =~ s/\$editLibraryIsPublicRadio2//g;
	$html =~ s/\$editLibraryIsPublicRadio1/checked="checked"/g;
}
if($libraryDetails->{'isFinished'})
{
	$html =~ s/\$editLibraryIsFinishedRadio2/checked="checked"/g;
	$html =~ s/\$editLibraryIsFinishedRadio1//g;
}
else
{
	$html =~ s/\$editLibraryIsFinishedRadio2//g;
	$html =~ s/\$editLibraryIsFinishedRadio1/checked="checked"/g;
}
if($libraryDetails->{'isLibraryForSale'})
{
	$html =~ s/\$editLibraryIsLibraryForSaleRadio2/checked="checked"/g;
	$html =~ s/\$editLibraryIsLibraryForSaleRadio1//g;
}
else
{
	$html =~ s/\$editLibraryIsLibraryForSaleRadio2//g;
	$html =~ s/\$editLibraryIsLibraryForSaleRadio1/checked="checked"/g;
}
if($libraryDetails->{'isFilterForSale'})
{
	$html =~ s/\$editLibraryIsFilterForSaleRadio2/checked="checked"/g;
	$html =~ s/\$editLibraryIsFilterForSaleRadio1//g;
}
else
{
	$html =~ s/\$editLibraryIsFilterForSaleRadio2//g;
	$html =~ s/\$editLibraryIsFilterForSaleRadio1/checked="checked"/g;
}
if($libraryDetails->{'isCloneForSale'})
{
	$html =~ s/\$editLibraryIsCloneForSaleRadio2/checked="checked"/g;
	$html =~ s/\$editLibraryIsCloneForSaleRadio1//g;
}
else
{
	$html =~ s/\$editLibraryIsCloneForSaleRadio2//g;
	$html =~ s/\$editLibraryIsCloneForSaleRadio1/checked="checked"/g;
}
if($library[7])
{
	$html =~ s/\$editLibraryIsSequencingRadio2/checked="checked"/g;
	$html =~ s/\$editLibraryIsSequencingRadio1//g;
}
else
{
	$html =~ s/\$editLibraryIsSequencingRadio2//g;
	$html =~ s/\$editLibraryIsSequencingRadio1/checked="checked"/g;
}
$html =~ s/\$libraryProvidedBy/$libraryDetails->{'providedBy'}/g;
$html =~ s/\$libraryOrganization/$libraryDetails->{'organization'}/g;
$html =~ s/\$libraryComments/$libraryDetails->{'comments'}/g;
$html =~ s/\$libraryOrderPageComments/$libraryDetails->{'orderPageComments'}/g;
$html =~ s/\$libraryCreator/$library[9]/g;
$html =~ s/\$libraryDateLibraryWasMade/$libraryDetails->{'dateLibraryWasMade'}/g;
$html =~ s/\$libraryDateLibraryWasAutoclaved/$libraryDetails->{'dateLibraryWasAutoclaved'}/g;
$html =~ s/\$libraryDateLibraryWasEntered/$library[10]/g;
if($library[5])
{
	my $checkSource = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$checkSource->execute($library[5]);
	my @checkSource=$checkSource->fetchrow_array();
	$library[5] = $checkSource[2];
}
$html =~ s/\$libraryRearrayingSource/$library[5]/g;
$html =~ s/\$projectId/$library[6]/g;
$html =~ s/\$projectName/$libraryToProject[2]/g;
$html =~ s/\$libraryArchiveStartDate/$libraryDetails->{'archiveStartDate'}/g;
$html =~ s/\$libraryArchiveEndDate/$libraryDetails->{'archiveEndDate'}/g;

print header;
print $html;

__DATA__
<form id="editLibrary" name="editLibrary" action="librarySave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="libraryId" id="editLibraryId" type="hidden" value="$libraryId" />
	<div id="editLibraryTabs">
		<ul>
		<li><a href="#editLibraryTabs-1">General</a></li>
		<li><a href="#editLibraryTabs-2">Materials</a></li>
		<li><a href="#editLibraryTabs-3">Specs</a></li>
		<li><a href="#editLibraryTabs-4">Status</a></li>
		<li><a href="#editLibraryTabs-5">Dates</a></li>
		</ul>
		<div id="editLibraryTabs-1">
			<table>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryName"><b>Library Name</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="name" id="editLibraryName" placeholder="Library Name" size="15" type="text" maxlength="9" value="$libraryName" readonly />
				in $projectName<input name="projectId" id="editLibraryProjectId" type="hidden" value="$projectId" />
				<br><sup class='ui-state-disabled'>Last changed by $libraryCreator on $libraryDateLibraryWasEntered</sup></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryNickname"><b>Nickname</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="nickname" id="editLibraryNickname" placeholder="Nickname" size="40" type="text" maxlength="135" value="$libraryNickname" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryType"><b>Library type</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="type" id="editLibraryType" placeholder="Type" size="40" type="text" maxlength="45" value="$libraryType" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryRearrayingSource"><b>Rearraying source</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="rearrayingSource" id="editLibraryRearrayingSource" placeholder="Rearraying Source" size="15" type="text" maxlength="9" value="$libraryRearrayingSource" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryFormat"><b>Plate format</b></label></td>
				<td><div id="editLibraryFormat">
					<input type="radio" id="editLibraryFormatRadio1" name="format" value="96" $editLibraryFormatRadio1><label for="editLibraryFormatRadio1">96 wells</label>
					<input type="radio" id="editLibraryFormatRadio2" name="format" value="384" $editLibraryFormatRadio2><label for="editLibraryFormatRadio2">384 wells</label>
					</div>
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryDistributorInstitute"><b>Distributor institute</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="distributorInstitute" id="editLibraryDistributorInstitute" placeholder="Distributor Institute" size="40" type="text" maxlength="25" value="$libraryDistributorInstitute" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryRefToPublication"><b>Ref. to publication</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="refToPublication" id="editLibraryRefToPublication" placeholder="Ref. To Publication" size="40" type="text" maxlength="255" value="$libraryRefToPublication" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryProvidedBy"><b>Provided by</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="providedBy" id="editLibraryProvidedBy" placeholder="Provider" size="40" type="text" maxlength="255" value="$libraryProvidedBy" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryOrganization"><b>Organization</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="organization" id="editLibraryOrganization" placeholder="Organization" size="40" type="text" maxlength="45" value="$libraryOrganization" /></td>
			</tr>
			</table>
		</div>
		<div id="editLibraryTabs-2">
			<table>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryGenus"><b>Genus Species Subspecies</b></label><br>(Common Name)</td>
				<td><INPUT class="ui-widget-content ui-corner-all" TYPE="text" name="genus" id="editLibraryGenus" placeholder="Genus" maxlength="45" size="12" value="$libraryGenus">
					<INPUT class="ui-widget-content ui-corner-all" TYPE="text" name="species" id="editLibrarySpecies" placeholder="Species" maxlength="45" size="12" value="$librarySpecies">
					<INPUT class="ui-widget-content ui-corner-all" TYPE="text" name="subspecies" id="editLibrarySubspecies" placeholder="Subspecies" maxlength="45" size="12" value="$librarySubspecies">
					<INPUT class="ui-widget-content ui-corner-all" TYPE="text" name="commonName" id="editLibraryCommonName" placeholder="Common Name" maxlength="45" size="40" value="$libraryCommonName">
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryCultivar"><b>Cultivar</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="cultivar" id="editLibraryCultivar" placeholder="Cultivar" size="40" type="text" maxlength="255" value="$libraryCultivar" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryGenomeEquivalents"><b>Genome equivalents</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="genomeEquivalents" id="editLibraryGenomeEquivalents" placeholder="Genome Equivalents" size="20" type="text" maxlength="13" value="$libraryGenomeEquivalents" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryGenomeSize"><b>Genome size</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="genomeSize" id="editLibraryGenomeSize" placeholder="Genome Size" size="20" type="text" maxlength="45" value="$libraryGenomeSize" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryGenomeType"><b>Genome type</b></label><br>(like CC, BBCC or AA)</td>
				<td><input class='ui-widget-content ui-corner-all' name="genomeType" id="editLibraryGenomeType" placeholder="Genome Type" size="20" type="text" maxlength="9" value="$libraryGenomeType" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryTissue"><b>Tissue</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="tissue" id="editLibraryTissue" placeholder="Tissue" size="40" type="text" maxlength="255" value="$libraryTissue" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryTreatment"><b>Treatment</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="treatment" id="editLibraryTreatment" placeholder="Treatment" size="40" type="text" maxlength="255" value="$libraryTreatment" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryDevelopmentStage"><b>Development stage</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="developmentStage" id="editLibraryDevelopmentStage" placeholder="Development Stage" size="40" type="text" maxlength="255" value="$libraryDevelopmentStage" /></td>
			</tr>
			</table>		
		</div>
		<div id="editLibraryTabs-3">
			<table>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryNofFilters"><b># of filters</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="nofFilters" id="editLibraryNofFilters" placeholder="Number of Filters" size="20" type="text" maxlength="4" value="$libraryNofFilters" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryAverageInsertSize"><b>Average insert size</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="averageInsertSize" id="editLibraryAverageInsertSize" placeholder="Average Insert Size" size="20" type="text" maxlength="45" value="$libraryAverageInsertSize" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryStandardDeviation"><b>Standard deviation</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="standardDeviation" id="editLibraryStandardDeviation" placeholder="Standard Deviation" size="10" type="text" maxlength="3" value="$libraryStandardDeviation" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryVector"><b>Vector</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="vector" id="editLibraryVector" placeholder="Vector" size="40" type="text" maxlength="45" value="$libraryVector" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryHost"><b>Host</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="host" id="editLibraryHost" placeholder="Host" size="40" type="text" maxlength="45" value="$libraryHost" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryEnzymeFivePrime"><b>Enzyme 5' prime</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="enzymeFivePrime" id="editLibraryEnzymeFivePrime" placeholder="Enzyme 5' Prime" size="40" type="text" maxlength="45" value="$libraryEnzymeFivePrime" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryEnzymeThreePrime"><b>Enzyme 3' prime</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="enzymeThreePrime" id="editLibraryEnzymeThreePrime" placeholder="Enzyme 3' Prime" size="40" type="text" maxlength="45" value="$libraryEnzymeThreePrime" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryEndSeqFivePrime"><b>End seq primer 5' prime</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="endSeqFivePrime" id="editLibraryEndSeqFivePrime" placeholder="End Seq 5' Prime" size="40" type="text" maxlength="45" value="$libraryEndSeqFivePrime" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryEndSeqThreePrime"><b>End seq primer 3' prime</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="endSeqThreePrime" id="editLibraryEndSeqThreePrime" placeholder="End Seq 3' Prime" size="40" type="text" maxlength="45" value="$libraryEndSeqThreePrime" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryCloningFivePrime"><b>Cloning linker 5' prime</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="cloningFivePrime" id="editLibraryCloningFivePrime" placeholder="Cloning linker 5' Prime" size="40" type="text" maxlength="45" value="$libraryCloningFivePrime" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryCloningThreePrime"><b>Cloning linker 3' prime</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="cloningThreePrime" id="editLibraryCloningThreePrime" placeholder="Cloning linker 3' Prime" size="40" type="text" maxlength="45" value="$libraryCloningThreePrime" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryAntibiotic"><b>Antibiotic</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="antibiotic" id="editLibraryAntibiotic" placeholder="Antibiotic" size="40" type="text" maxlength="45" value="$libraryAntibiotic" /></td>
			</tr>
			</table>
		</div>
		<div id="editLibraryTabs-4">
			<table>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryStatus"><b>Status</b></label></td>
				<td><div id="editLibraryStatus">
					<input type="radio" id="editLibraryStatusRadio1" name="status" value="active" $editLibraryStatusRadio1><label for="editLibraryStatusRadio1">Active</label>
					<input type="radio" id="editLibraryStatusRadio2" name="status" value="inactive" $editLibraryStatusRadio2><label for="editLibraryStatusRadio2">Inactive</label>
					<input type="radio" id="editLibraryStatusRadio3" name="status" value="deleted" $editLibraryStatusRadio3><label for="editLibraryStatusRadio3">Deleted</label></div>
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryIsExternal"><b>Is external?</b></label><br>(Will not store plates)</td>
				<td><div id="editLibraryIsExternal">
					<input type="radio" id="editLibraryIsExternalRadio2" name="isExternal" value="1" $editLibraryIsExternalRadio2><label for="editLibraryIsExternalRadio2">Yes</label>
					<input type="radio" id="editLibraryIsExternalRadio1" name="isExternal" value="0" $editLibraryIsExternalRadio1><label for="editLibraryIsExternalRadio1">No</label>
					</div>
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryIsPublic"><b>Is public?</b></label></td>
				<td><div id="editLibraryIsPublic">
					<input type="radio" id="editLibraryIsPublicRadio2" name="isPublic" value="1" $editLibraryIsPublicRadio2><label for="editLibraryIsPublicRadio2">Yes</label>
					<input type="radio" id="editLibraryIsPublicRadio1" name="isPublic" value="0" $editLibraryIsPublicRadio1><label for="editLibraryIsPublicRadio1">No</label>
					</div>
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryIsFinished"><b>Is finished?</b></label></td>
				<td><div id="editLibraryIsFinished">
					<input type="radio" id="editLibraryIsFinishedRadio2" name="isFinished" value="1" $editLibraryIsFinishedRadio2><label for="editLibraryIsFinishedRadio2">Yes</label>
					<input type="radio" id="editLibraryIsFinishedRadio1" name="isFinished" value="0" $editLibraryIsFinishedRadio1><label for="editLibraryIsFinishedRadio1">No</label>
					</div>
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryIsLibraryForSale"><b>Is library for sale?</b></label></td>
				<td><div id="editLibraryIsLibraryForSale">
					<input type="radio" id="editLibraryIsLibraryForSaleRadio2" name="isLibraryForSale" value="1" $editLibraryIsLibraryForSaleRadio2><label for="editLibraryIsLibraryForSaleRadio2">Yes</label>
					<input type="radio" id="editLibraryIsLibraryForSaleRadio1" name="isLibraryForSale" value="0" $editLibraryIsLibraryForSaleRadio1><label for="editLibraryIsLibraryForSaleRadio1">No</label>
					</div>
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryIsFilterForSale"><b>Is filter for sale?</b></label></td>
				<td><div id="editLibraryIsFilterForSale">
					<input type="radio" id="editLibraryIsFilterForSaleRadio2" name="isFilterForSale" value="1" $editLibraryIsFilterForSaleRadio2><label for="editLibraryIsFilterForSaleRadio2">Yes</label>
					<input type="radio" id="editLibraryIsFilterForSaleRadio1" name="isFilterForSale" value="0" $editLibraryIsFilterForSaleRadio1><label for="editLibraryIsFilterForSaleRadio1">No</label>
					</div>
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryIsCloneForSale"><b>Is clone for sale?</b></label></td>
				<td><div id="editLibraryIsCloneForSale">
					<input type="radio" id="editLibraryIsCloneForSaleRadio2" name="isCloneForSale" value="1" $editLibraryIsCloneForSaleRadio2><label for="editLibraryIsCloneForSaleRadio2">Yes</label>
					<input type="radio" id="editLibraryIsCloneForSaleRadio1" name="isCloneForSale" value="0" $editLibraryIsCloneForSaleRadio1><label for="editLibraryIsCloneForSaleRadio1">No</label>
					</div>
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryIsSequencing"><b>Is library for sequencing?</b></label></td>
				<td><div id="editLibraryIsSequencing">
					<input type="radio" id="editLibraryIsSequencingRadio2" name="isSequencing" value="1" $editLibraryIsSequencingRadio2><label for="editLibraryIsSequencingRadio2">Yes</label>
					<input type="radio" id="editLibraryIsSequencingRadio1" name="isSequencing" value="0" $editLibraryIsSequencingRadio1><label for="editLibraryIsSequencingRadio1">No</label>
					</div>
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryOrderPageComments"><b>Order page comments</b></label></td>
				<td><textarea class='ui-widget-content ui-corner-all' name="orderPageComments" id="editLibraryOrderPageComments" cols="40" rows="5">$libraryOrderPageComments</textarea></td>
			</tr>
			</table>		
		</div>
		<div id="editLibraryTabs-5">
			<table>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryDateLibraryWasMade"><b>Library was made on</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="dateLibraryWasMade" id="editLibraryDateLibraryWasMade" placeholder="YYYY-MM-DD" size="20" type="text" maxlength="19" value="$libraryDateLibraryWasMade" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryDateLibraryWasAutoclaved"><b>Library was autoclaved on</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="dateLibraryWasAutoclaved" id="editLibraryDateLibraryWasAutoclaved" placeholder="YYYY-MM-DD" size="20" type="text" maxlength="19" value="$libraryDateLibraryWasAutoclaved" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryArchiveStartDate"><b>Archive start date</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="archiveStartDate" id="editLibraryArchiveStartDate" placeholder="YYYY-MM-DD" size="20" type="text" maxlength="19" value="$libraryArchiveStartDate" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryArchiveEndDate"><b>Archive end date</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="archiveEndDate" id="editLibraryArchiveEndDate" placeholder="YYYY-MM-DD" size="20" type="text" maxlength="19" value="$libraryArchiveEndDate" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="editLibraryComments"><b>Comments</b></label></td>
				<td><textarea class='ui-widget-content ui-corner-all' name="comments" id="editLibraryComments" cols="40" rows="5">$libraryComments</textarea></td>
			</tr>
			</table>
		</div>
	</div>
</form>
<script>
$( "#editLibraryTabs" ).tabs();
$( "#editLibraryFormat" ).buttonset();
$( "#editLibraryIsExternal" ).buttonset();
$( "#editLibraryStatus" ).buttonset();
$( "#editLibraryIsPublic" ).buttonset();
$( "#editLibraryIsFinished" ).buttonset();
$( "#editLibraryIsLibraryForSale" ).buttonset();
$( "#editLibraryIsFilterForSale" ).buttonset();
$( "#editLibraryIsCloneForSale" ).buttonset();
$( "#editLibraryIsSequencing" ).buttonset();
$( "#editLibraryRearrayingSource" ).autocomplete({
	source: "autoNameSearch.cgi",
	minLength: 0,
	select: function( event, ui ) {
		if( !$("#editLibraryGenomeEquivalents").val() ) {$( "#editLibraryGenomeEquivalents" ).val( ui.item.genomeEquivalents );}
		if( !$("#editLibraryGenomeSize").val() ) {$( "#editLibraryGenomeSize" ).val( ui.item.genomeSize );}
		if( !$("#editLibraryAverageInsertSize").val() ) {$( "#editLibraryAverageInsertSize" ).val( ui.item.averageInsertSize );}
		if( !$("#editLibraryStandardDeviation").val() ) {$( "#editLibraryStandardDeviation" ).val( ui.item.standardDeviation );}
		if( !$("#editLibraryDistributorInstitute").val() ) {$( "#editLibraryDistributorInstitute" ).val( ui.item.distributorInstitute );}
		if( !$("#editLibraryVector").val() ) {$( "#editLibraryVector" ).val( ui.item.vector );}
		if( !$("#editLibraryHost").val() ) {$( "#editLibraryHost" ).val( ui.item.host );}
		if( !$("#editLibraryEnzymeFivePrime").val() ) {$( "#editLibraryEnzymeFivePrime" ).val( ui.item.enzymeFivePrime );}
		if( !$("#editLibraryEnzymeThreePrime").val() ) {$( "#editLibraryEnzymeThreePrime" ).val( ui.item.enzymeThreePrime );}
		if( !$("#editLibraryEndSeqFivePrime").val() ) {$( "#editLibraryEndSeqFivePrime" ).val( ui.item.endSeqFivePrime );}
		if( !$("#editLibraryEndSeqThreePrime").val() ) {$( "#editLibraryEndSeqThreePrime" ).val( ui.item.endSeqThreePrime );}
		if( !$("#editLibraryCloningFivePrime").val() ) {$( "#editLibraryCloningFivePrime" ).val( ui.item.cloningFivePrime );}
		if( !$("#editLibraryCloningThreePrime").val() ) {$( "#editLibraryCloningThreePrime" ).val( ui.item.cloningThreePrime );}
		if( !$("#editLibraryAntibiotic").val() ) {$( "#editLibraryAntibiotic" ).val( ui.item.antibiotic );}
		if( !$("#editLibraryGenus").val() ) {$( "#editLibraryGenus" ).val( ui.item.genus );}
		if( !$("#editLibrarySpecies").val() ) {$( "#editLibrarySpecies" ).val( ui.item.species );}
		if( !$("#editLibrarySubspecies").val() ) {$( "#editLibrarySubspecies" ).val( ui.item.subspecies );}
		if( !$("#editLibraryCommonName").val() ) {$( "#editLibraryCommonName" ).val( ui.item.commonName );}
		if( !$("#editLibraryGenomeType").val() ) {$( "#editLibraryGenomeType" ).val( ui.item.genomeType );}
		if( !$("#editLibraryTissue").val() ) {$( "#editLibraryTissue" ).val( ui.item.tissue );}
		if( !$("#editLibraryTreatment").val() ) {$( "#editLibraryTreatment" ).val( ui.item.treatment );}
		if( !$("#editLibraryDevelopmentStage").val() ) {$( "#editLibraryDevelopmentStage" ).val( ui.item.developmentStage );}
		if( !$("#editLibraryCultivar").val() ) {$( "#editLibraryCultivar" ).val( ui.item.cultivar );}
		if( !$("#editLibraryRefToPublication").val() ) {$( "#editLibraryRefToPublication" ).val( ui.item.refToPublication );}
		if( !$("#editLibraryProvidedBy").val() ) {$( "#editLibraryProvidedBy" ).val( ui.item.providedBy );}
		if( !$("#editLibraryOrganization").val() ) {$( "#editLibraryOrganization" ).val( ui.item.organization );}
		}
	});
$( "#editLibraryType" ).autocomplete({
	source: "autoTypeSearch.cgi",
	minLength: 0
	});
$( "#editLibraryDistributorInstitute" ).autocomplete({
	source: "autoInstituteSearch.cgi",
	minLength: 0
	});
$( "#editLibraryGenus" ).autocomplete({
	source: "autoGenusSearch.cgi",
	minLength: 0
	});
$( "#editLibrarySpecies" ).autocomplete({
	source: "autoSpeciesSearch.cgi",
	minLength: 0
	});
$( "#editLibrarySubspecies" ).autocomplete({
	source: "autoSubspeciesSearch.cgi",
	minLength: 0
	});
$( "#editLibraryCommonName" ).autocomplete({
	source: "autoCommonNameSearch.cgi",
	minLength: 0
	});
$( "#editLibraryVector" ).autocomplete({
	source: "autoVectorSearch.cgi",
	minLength: 0
	});
$( "#editLibraryHost" ).autocomplete({
	source: "autoHostSearch.cgi",
	minLength: 0
	});
$( "#editLibraryAntibiotic" ).autocomplete({
	source: "autoAntibioticSearch.cgi",
	minLength: 0
	});
$( "#editLibraryProvidedBy" ).autocomplete({
	source: "autoProviderSearch.cgi",
	minLength: 0
	});
$( "#editLibraryOrganization" ).autocomplete({
	source: "autoOrganizationSearch.cgi",
	minLength: 0
	});
$( "#editLibraryDateLibraryWasMade" ).datepicker({dateFormat:"yy-mm-dd"});
$( "#editLibraryDateLibraryWasAutoclaved" ).datepicker({dateFormat:"yy-mm-dd"});
$( "#editLibraryArchiveStartDate" ).datepicker({dateFormat:"yy-mm-dd"});
$( "#editLibraryArchiveEndDate" ).datepicker({dateFormat:"yy-mm-dd"});
$('#dialog').dialog("option", "title", "Edit Library");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('editLibrary'); } }, { text: "Delete", click: function() { deleteItem($libraryId); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
</script>


