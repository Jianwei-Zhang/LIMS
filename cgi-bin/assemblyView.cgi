#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use JSON::XS; #JSON::XS is recommended to be installed for handling JSON string of big size 
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

my $assemblyId = param ('assemblyId') || '';
my $assembly=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$assembly->execute($assemblyId);
my @assembly = $assembly->fetchrow_array();
my $target=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$target->execute($assembly[4]);
my @target = $target->fetchrow_array();

my $fpcOrAgpId = '';
if($assembly[6])
{
	my $fpcOrAgpList=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$fpcOrAgpList->execute($assembly[6]);
	my @fpcOrAgpList = $fpcOrAgpList->fetchrow_array();
	if($fpcOrAgpList[1] eq 'fpc')
	{
		$fpcOrAgpId .= "<tr><td style='text-align:right'><b>Physical Reference</b></td><td title='$fpcOrAgpList[8]'>FPC: $fpcOrAgpList[2] v.$fpcOrAgpList[3]</td></tr>";
	}
	elsif($fpcOrAgpList[1]  eq 'agp')
	{
		$fpcOrAgpId .= "<tr><td style='text-align:right'><b>Physical Reference</b></td><td>AGP: $fpcOrAgpList[2] v.$fpcOrAgpList[3]</td></tr>";
	}
}

my $refGenomeId = '';
if($assembly[5])
{
	my $genomeList=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$genomeList->execute($assembly[5]);
	my @genomeList = $genomeList->fetchrow_array();
	$refGenomeId = "<tr><td style='text-align:right'><b>Reference Genome</b></td><td title='$genomeList[8]'>$genomeList[2]</td></tr>";
}

my $assemblyExtraIds = "<tr><td style='text-align:right'><b>Extra Genome</b></td><td>";
my $checkAsbGenome = $dbh->prepare("SELECT * FROM matrix,link WHERE link.type LIKE 'asbGenome' AND link.child = matrix.id AND link.parent = ?");
$checkAsbGenome->execute($assemblyId);
if ($checkAsbGenome->rows > 0)
{
	while(my @checkAsbGenome=$checkAsbGenome->fetchrow_array())
	{
		$assemblyExtraIds .= "<a onclick='closeDialog();openDialog(\"genomeView.cgi?genomeId=$checkAsbGenome[0]\")' title='$checkAsbGenome[8]'>$checkAsbGenome[2]<span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-arrow-1-ne'></span></a><br>";
	}
	$assemblyExtraIds .= "</td></tr>";
}
else
{
	$assemblyExtraIds .= "None.</td></tr>";
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$assemblyId/$assemblyId/g;
$html =~ s/\$assemblyName/$assembly[2]/g;
$html =~ s/\$assemblyVersion/$assembly[3]/g;
$target[1] = ucfirst ($target[1]);
$html =~ s/\$sourceType/$target[1]/g;
$html =~ s/\$sourceName/$target[2]/g;
$html =~ s/\$fpcOrAgpId/$fpcOrAgpId/g;
$html =~ s/\$refGenomeId/$refGenomeId/g;
$html =~ s/\$assemblyExtraIds/$assemblyExtraIds/g;

my $assemblyDetails = decode_json $assembly[8];
$assemblyDetails->{'description'} = '' if (!exists $assemblyDetails->{'description'});
$html =~ s/\$assemblyDescription/$assemblyDetails->{'description'}/g;

$html =~ s/\$assemblyCreator/$assembly[9]/g;
$html =~ s/\$assemblyCreationDate/$assembly[10]/g;
print header;
print $html;

__DATA__
<table>
	<tr><td style='text-align:right'><b>Assembly Name</b></td><td>$assemblyName<br>Version $assemblyVersion <sup class='ui-state-disabled'>by $assemblyCreator on $assemblyCreationDate</sup></td></tr>
	<tr><td style='text-align:right'><b>Source $sourceType</b></td><td>$sourceName</td></tr>
	$fpcOrAgpId
	$refGenomeId
	$assemblyExtraIds
	<tr><td style='text-align:right'><b>Description</b></td><td><textarea class='ui-widget-content ui-corner-all' name="description" id="viewAssemblyDescription" cols="50" rows="10" placeholder="Give some information about this assembly. Or you may do it later." readonly="readonly">$assemblyDescription</textarea></td></tr>
	<tr><td style='text-align:right'></td><td>
		<form id="copyAssembly" name="copyAssembly" action="assemblyCopy.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
			<input type="checkbox" name="assemblyId" id="copyAssemblyIdCheckbox" value="$assemblyId" /><b><label for="copyAssemblyIdCheckbox">Make a copy of this assembly</label></b>
			<input type="hidden" name="action" value="copy">
		</form>
	</td></tr>
</table>
<script>
$('#dialog').dialog("option", "title", "View Assembly $assemblyName");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Copy Assembly", click: function() { submitForm('copyAssembly'); } }, { text: "OK", click: function() {closeDialog(); } } ] );
</script>