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

my $projectId = param ('projectId') || '';
my $project=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$project->execute($projectId);
my @project = $project->fetchrow_array();

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$projectId/$projectId/g;
$html =~ s/\$projectName/$project[2]/g;

print header;
print $html;

__DATA__
<form id="newLibrary" name="newLibrary" action="librarySave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="projectId" id="newLibraryProjectId" type="hidden" value="$projectId" />
	<div id="newLibraryTabs">
		<ul>
		<li><a href="#newLibraryTabs-1">General</a></li>
		<li><a href="#newLibraryTabs-2">Materials</a></li>
		<li><a href="#newLibraryTabs-3">Specs</a></li>
		<li><a href="#newLibraryTabs-4">Status</a></li>
		<li><a href="#newLibraryTabs-5">Dates</a></li>
		</ul>
		<div id="newLibraryTabs-1">
			<table>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryName"><b>Library Name</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="name" id="newLibraryName" placeholder="Library Name" size="15" type="text" maxlength="9" />
				in <b>$projectName</b></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryNickname"><b>Nickname</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="nickname" id="newLibraryNickname" placeholder="Nickname" size="40" type="text" maxlength="135" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryType"><b>Library type</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="type" id="newLibraryType" placeholder="Type" size="40" type="text" maxlength="45" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryRearrayingSource"><b>Rearraying source</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="rearrayingSource" id="newLibraryRearrayingSource" placeholder="Rearraying Source" size="15" type="text" maxlength="9" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryFormat"><b>Plate format</b></label></td>
				<td><div id="newLibraryFormat">
					<input type="radio" id="newLibraryFormatRadio1" name="format" value="96"><label for="newLibraryFormatRadio1">96 wells</label>
					<input type="radio" id="newLibraryFormatRadio2" name="format" value="384" checked="checked"><label for="newLibraryFormatRadio2">384 wells</label>
					</div>
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryDistributorInstitute"><b>Distributor institute</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="distributorInstitute" id="newLibraryDistributorInstitute" placeholder="Distributor Institute" size="40" type="text" maxlength="25" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryRefToPublication"><b>Ref. to publication</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="refToPublication" id="newLibraryRefToPublication" placeholder="Ref. To Publication" size="40" type="text" maxlength="255" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryProvidedBy"><b>Provided by</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="providedBy" id="newLibraryProvidedBy" placeholder="Provider" size="40" type="text" maxlength="255" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryOrganization"><b>Organization</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="organization" id="newLibraryOrganization" placeholder="Organization" size="40" type="text" maxlength="45" /></td>
			</tr>
			</table>
		</div>
		<div id="newLibraryTabs-2">
			<table>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryGenus"><b>Genus Species Subspecies</b></label><br>(Common Name)</td>
				<td><INPUT class="ui-widget-content ui-corner-all" TYPE="text" name="genus" id="newLibraryGenus" placeholder="Genus" maxlength="45" size="12">
					<INPUT class="ui-widget-content ui-corner-all" TYPE="text" name="species" id="newLibrarySpecies" placeholder="Species" maxlength="45" size="12">
					<INPUT class="ui-widget-content ui-corner-all" TYPE="text" name="subspecies" id="newLibrarySubspecies" placeholder="Subspecies" maxlength="45" size="12">
					<INPUT class="ui-widget-content ui-corner-all" TYPE="text" name="commonName" id="newLibraryCommonName" placeholder="Common Name" maxlength="45" size="40">
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryCultivar"><b>Cultivar</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="cultivar" id="newLibraryCultivar" placeholder="Cultivar" size="40" type="text" maxlength="255" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryGenomeEquivalents"><b>Genome equivalents</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="genomeEquivalents" id="newLibraryGenomeEquivalents" placeholder="Genome Equivalents" size="20" type="text" maxlength="13" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryGenomeSize"><b>Genome size</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="genomeSize" id="newLibraryGenomeSize" placeholder="Genome Size" size="20" type="text" maxlength="45" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryGenomeType"><b>Genome type</b></label><br>(like CC, BBCC or AA)</td>
				<td><input class='ui-widget-content ui-corner-all' name="genomeType" id="newLibraryGenomeType" placeholder="Genome Type" size="20" type="text" maxlength="9" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryTissue"><b>Tissue</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="tissue" id="newLibraryTissue" placeholder="Tissue" size="40" type="text" maxlength="255" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryTreatment"><b>Treatment</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="treatment" id="newLibraryTreatment" placeholder="Treatment" size="40" type="text" maxlength="255" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryDevelopmentStage"><b>Development stage</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="developmentStage" id="newLibraryDevelopmentStage" placeholder="Development Stage" size="40" type="text" maxlength="255" /></td>
			</tr>
			</table>		
		</div>
		<div id="newLibraryTabs-3">
			<table>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryNofFilters"><b># of filters</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="nofFilters" id="newLibraryNofFilters" placeholder="Number of Filters" size="20" type="text" maxlength="4" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryAverageInsertSize"><b>Average insert size</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="averageInsertSize" id="newLibraryAverageInsertSize" placeholder="Average Insert Size" size="20" type="text" maxlength="45" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryStandardDeviation"><b>Standard deviation</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="standardDeviation" id="newLibraryStandardDeviation" placeholder="Standard Deviation" size="10" type="text" maxlength="3" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryVector"><b>Vector</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="vector" id="newLibraryVector" placeholder="Vector" size="40" type="text" maxlength="45" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryHost"><b>Host</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="host" id="newLibraryHost" placeholder="Host" size="40" type="text" maxlength="45" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryEnzymeFivePrime"><b>Enzyme 5' prime</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="enzymeFivePrime" id="newLibraryEnzymeFivePrime" placeholder="Enzyme 5' Prime" size="40" type="text" maxlength="45" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryEnzymeThreePrime"><b>Enzyme 3' prime</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="enzymeThreePrime" id="newLibraryEnzymeThreePrime" placeholder="Enzyme 3' Prime" size="40" type="text" maxlength="45" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryEndSeqFivePrime"><b>End seq primer 5' prime</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="endSeqFivePrime" id="newLibraryEndSeqFivePrime" placeholder="End Seq 5' Prime" size="40" type="text" maxlength="45" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryEndSeqThreePrime"><b>End seq primer 3' prime</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="endSeqThreePrime" id="newLibraryEndSeqThreePrime" placeholder="End Seq 3' Prime" size="40" type="text" maxlength="45" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryCloningFivePrime"><b>Cloning linker 5' prime</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="cloningFivePrime" id="newLibraryCloningFivePrime" placeholder="Cloning linker 5' Prime" size="40" type="text" maxlength="45" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryCloningThreePrime"><b>Cloning linker 3' prime</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="cloningThreePrime" id="newLibraryCloningThreePrime" placeholder="Cloning linker 3' Prime" size="40" type="text" maxlength="45" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryAntibiotic"><b>Antibiotic</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="antibiotic" id="newLibraryAntibiotic" placeholder="Antibiotic" size="40" type="text" maxlength="45" /></td>
			</tr>
			</table>
		</div>
		<div id="newLibraryTabs-4">
			<table>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryStatus"><b>Status</b></label></td>
				<td><div id="newLibraryStatus">
					<input type="radio" id="newLibraryStatusRadio1" name="status" value="active" checked="checked"><label for="newLibraryStatusRadio1">Active</label>
					<input type="radio" id="newLibraryStatusRadio2" name="status" value="inactive"><label for="newLibraryStatusRadio2">Inactive</label>
					</div>
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryIsExternal"><b>Is external?</b></label><br>(Will not store plates)</td>
				<td><div id="newLibraryIsExternal">
					<input type="radio" id="newLibraryIsExternalRadio2" name="isExternal" value="1"><label for="newLibraryIsExternalRadio2">Yes</label>
					<input type="radio" id="newLibraryIsExternalRadio1" name="isExternal" value="0" checked="checked"><label for="newLibraryIsExternalRadio1">No</label>
					</div>
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryIsPublic"><b>Is public?</b></label></td>
				<td><div id="newLibraryIsPublic">
					<input type="radio" id="newLibraryIsPublicRadio2" name="isPublic" value="1"><label for="newLibraryIsPublicRadio2">Yes</label>
					<input type="radio" id="newLibraryIsPublicRadio1" name="isPublic" value="0" checked="checked"><label for="newLibraryIsPublicRadio1">No</label>
					</div>
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryIsFinished"><b>Is finished?</b></label></td>
				<td><div id="newLibraryIsFinished">
					<input type="radio" id="newLibraryIsFinishedRadio2" name="isFinished" value="1"><label for="newLibraryIsFinishedRadio2">Yes</label>
					<input type="radio" id="newLibraryIsFinishedRadio1" name="isFinished" value="0" checked="checked"><label for="newLibraryIsFinishedRadio1">No</label>
					</div>
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryIsLibraryForSale"><b>Is library for sale?</b></label></td>
				<td><div id="newLibraryIsLibraryForSale">
					<input type="radio" id="newLibraryIsLibraryForSaleRadio2" name="isLibraryForSale" value="1"><label for="newLibraryIsLibraryForSaleRadio2">Yes</label>
					<input type="radio" id="newLibraryIsLibraryForSaleRadio1" name="isLibraryForSale" value="0" checked="checked"><label for="newLibraryIsLibraryForSaleRadio1">No</label>
					</div>
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryIsFilterForSale"><b>Is filter for sale?</b></label></td>
				<td><div id="newLibraryIsFilterForSale">
					<input type="radio" id="newLibraryIsFilterForSaleRadio2" name="isFilterForSale" value="1"><label for="newLibraryIsFilterForSaleRadio2">Yes</label>
					<input type="radio" id="newLibraryIsFilterForSaleRadio1" name="isFilterForSale" value="0" checked="checked"><label for="newLibraryIsFilterForSaleRadio1">No</label>
					</div>
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryIsCloneForSale"><b>Is clone for sale?</b></label></td>
				<td><div id="newLibraryIsCloneForSale">
					<input type="radio" id="newLibraryIsCloneForSaleRadio2" name="isCloneForSale" value="1"><label for="newLibraryIsCloneForSaleRadio2">Yes</label>
					<input type="radio" id="newLibraryIsCloneForSaleRadio1" name="isCloneForSale" value="0" checked="checked"><label for="newLibraryIsCloneForSaleRadio1">No</label>
					</div>
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryIsSequencing"><b>Is library for sequencing?</b></label></td>
				<td><div id="newLibraryIsSequencing">
					<input type="radio" id="newLibraryIsSequencingRadio2" name="isSequencing" value="1"><label for="newLibraryIsSequencingRadio2">Yes</label>
					<input type="radio" id="newLibraryIsSequencingRadio1" name="isSequencing" value="0" checked="checked"><label for="newLibraryIsSequencingRadio1">No</label>
					</div>
				</td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryOrderPageComments"><b>Order page comments</b></label></td>
				<td><textarea class='ui-widget-content ui-corner-all' name="orderPageComments" id="newLibraryOrderPageComments" cols="40" rows="5"></textarea></td>
			</tr>
			</table>		
		</div>
		<div id="newLibraryTabs-5">
			<table>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryDateLibraryWasMade"><b>Library was made on</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="dateLibraryWasMade" id="newLibraryDateLibraryWasMade" placeholder="YYYY-MM-DD" size="20" type="text" maxlength="19" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryDateLibraryWasAutoclaved"><b>Library was autoclaved on</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="dateLibraryWasAutoclaved" id="newLibraryDateLibraryWasAutoclaved" placeholder="YYYY-MM-DD" size="20" type="text" maxlength="19" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryArchiveStartDate"><b>Archive start date</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="archiveStartDate" id="newLibraryArchiveStartDate" placeholder="YYYY-MM-DD" size="20" type="text" maxlength="19" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryArchiveEndDate"><b>Archive end date</b></label></td>
				<td><input class='ui-widget-content ui-corner-all' name="archiveEndDate" id="newLibraryArchiveEndDate" placeholder="YYYY-MM-DD" size="20" type="text" maxlength="19" /></td>
			</tr>
			<tr><td style='text-align:right;white-space: nowrap;'><label for="newLibraryComments"><b>Comments</b></label></td>
				<td><textarea class='ui-widget-content ui-corner-all' name="comments" id="newLibraryComments" cols="40" rows="5"></textarea></td>
			</tr>
			</table>
		</div>
	</div>
</form>
<script>
$( "#newLibraryTabs" ).tabs();
$( "#newLibraryFormat" ).buttonset();
$( "#newLibraryIsExternal" ).buttonset();
$( "#newLibraryStatus" ).buttonset();
$( "#newLibraryIsPublic" ).buttonset();
$( "#newLibraryIsFinished" ).buttonset();
$( "#newLibraryIsLibraryForSale" ).buttonset();
$( "#newLibraryIsFilterForSale" ).buttonset();
$( "#newLibraryIsCloneForSale" ).buttonset();
$( "#newLibraryIsSequencing" ).buttonset();
$( "#newLibraryRearrayingSource" ).autocomplete({
	source: "autoNameSearch.cgi",
	minLength: 0,
	select: function( event, ui ) {
		if( !$("#newLibraryGenomeEquivalents").val() ) {$( "#newLibraryGenomeEquivalents" ).val( ui.item.genomeEquivalents );}
		if( !$("#newLibraryGenomeSize").val() ) {$( "#newLibraryGenomeSize" ).val( ui.item.genomeSize );}
		if( !$("#newLibraryAverageInsertSize").val() ) {$( "#newLibraryAverageInsertSize" ).val( ui.item.averageInsertSize );}
		if( !$("#newLibraryStandardDeviation").val() ) {$( "#newLibraryStandardDeviation" ).val( ui.item.standardDeviation );}
		if( !$("#newLibraryDistributorInstitute").val() ) {$( "#newLibraryDistributorInstitute" ).val( ui.item.distributorInstitute );}
		if( !$("#newLibraryVector").val() ) {$( "#newLibraryVector" ).val( ui.item.vector );}
		if( !$("#newLibraryHost").val() ) {$( "#newLibraryHost" ).val( ui.item.host );}
		if( !$("#newLibraryEnzymeFivePrime").val() ) {$( "#newLibraryEnzymeFivePrime" ).val( ui.item.enzymeFivePrime );}
		if( !$("#newLibraryEnzymeThreePrime").val() ) {$( "#newLibraryEnzymeThreePrime" ).val( ui.item.enzymeThreePrime );}
		if( !$("#newLibraryEndSeqFivePrime").val() ) {$( "#newLibraryEndSeqFivePrime" ).val( ui.item.endSeqFivePrime );}
		if( !$("#newLibraryEndSeqThreePrime").val() ) {$( "#newLibraryEndSeqThreePrime" ).val( ui.item.endSeqThreePrime );}
		if( !$("#newLibraryCloningFivePrime").val() ) {$( "#newLibraryCloningFivePrime" ).val( ui.item.cloningFivePrime );}
		if( !$("#newLibraryCloningThreePrime").val() ) {$( "#newLibraryCloningThreePrime" ).val( ui.item.cloningThreePrime );}
		if( !$("#newLibraryAntibiotic").val() ) {$( "#newLibraryAntibiotic" ).val( ui.item.antibiotic );}
		if( !$("#newLibraryGenus").val() ) {$( "#newLibraryGenus" ).val( ui.item.genus );}
		if( !$("#newLibrarySpecies").val() ) {$( "#newLibrarySpecies" ).val( ui.item.species );}
		if( !$("#newLibrarySubspecies").val() ) {$( "#newLibrarySubspecies" ).val( ui.item.subspecies );}
		if( !$("#newLibraryCommonName").val() ) {$( "#newLibraryCommonName" ).val( ui.item.commonName );}
		if( !$("#newLibraryGenomeType").val() ) {$( "#newLibraryGenomeType" ).val( ui.item.genomeType );}
		if( !$("#newLibraryTissue").val() ) {$( "#newLibraryTissue" ).val( ui.item.tissue );}
		if( !$("#newLibraryTreatment").val() ) {$( "#newLibraryTreatment" ).val( ui.item.treatment );}
		if( !$("#newLibraryDevelopmentStage").val() ) {$( "#newLibraryDevelopmentStage" ).val( ui.item.developmentStage );}
		if( !$("#newLibraryCultivar").val() ) {$( "#newLibraryCultivar" ).val( ui.item.cultivar );}
		if( !$("#newLibraryRefToPublication").val() ) {$( "#newLibraryRefToPublication" ).val( ui.item.refToPublication );}
		if( !$("#newLibraryProvidedBy").val() ) {$( "#newLibraryProvidedBy" ).val( ui.item.providedBy );}
		if( !$("#newLibraryOrganization").val() ) {$( "#newLibraryOrganization" ).val( ui.item.organization );}
		}
	});
$( "#newLibraryType" ).autocomplete({
	source: "autoTypeSearch.cgi",
	minLength: 0
	});
$( "#newLibraryDistributorInstitute" ).autocomplete({
	source: "autoInstituteSearch.cgi",
	minLength: 0
	});
$( "#newLibraryGenus" ).autocomplete({
	source: "autoGenusSearch.cgi",
	minLength: 0
	});
$( "#newLibrarySpecies" ).autocomplete({
	source: "autoSpeciesSearch.cgi",
	minLength: 0
	});
$( "#newLibrarySubspecies" ).autocomplete({
	source: "autoSubspeciesSearch.cgi",
	minLength: 0
	});
$( "#newLibraryCommonName" ).autocomplete({
	source: "autoCommonNameSearch.cgi",
	minLength: 0
	});
$( "#newLibraryVector" ).autocomplete({
	source: "autoVectorSearch.cgi",
	minLength: 0
	});
$( "#newLibraryHost" ).autocomplete({
	source: "autoHostSearch.cgi",
	minLength: 0
	});
$( "#newLibraryAntibiotic" ).autocomplete({
	source: "autoAntibioticSearch.cgi",
	minLength: 0
	});
$( "#newLibraryProvidedBy" ).autocomplete({
	source: "autoProviderSearch.cgi",
	minLength: 0
	});
$( "#newLibraryOrganization" ).autocomplete({
	source: "autoOrganizationSearch.cgi",
	minLength: 0
	});
$( "#newLibraryDateLibraryWasMade" ).datepicker({dateFormat:"yy-mm-dd"});
$( "#newLibraryDateLibraryWasAutoclaved" ).datepicker({dateFormat:"yy-mm-dd"});
$( "#newLibraryArchiveStartDate" ).datepicker({dateFormat:"yy-mm-dd"});
$( "#newLibraryArchiveEndDate" ).datepicker({dateFormat:"yy-mm-dd"});
$('#dialog').dialog("option", "title", "New Library");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('newLibrary'); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
</script>
