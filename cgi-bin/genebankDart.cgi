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

my $dartId = param ('dartId') || '';
my $dart = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$dart->execute($dartId);
my @dart=$dart->fetchrow_array();

#my @headerRows = ("genotype name","Order number","DArT plate barcode","client plate barcode","well row position","well column position","sample comments");
my @headerRows = ("genotype name");
my $genotypes = "<form id='dartList$$' name='dartList$$'>
	<table id='dart$$' class='display'>
		<thead>
			<tr>
				<th>G1
					<input type='checkbox' id='checkAllBox$$' name='checkAllBox$$' value='Check all' checked='checked' onClick='checkAll(\"genotypeIdOne\");return false;' title='Check all'>
					<input type='checkbox' id='uncheckAllBox$$' name='uncheckAllBox$$' value='Uncheck all' onClick='uncheckAll(\"genotypeIdOne\");return false;' title='Uncheck all'>
				</th>
				<th>G2
					<input type='checkbox' id='checkAllBox$$' name='checkAllBox$$' value='Check all' checked='checked' onClick='checkAll(\"genotypeIdTwo\");return false;' title='Check all'>
					<input type='checkbox' id='uncheckAllBox$$' name='uncheckAllBox$$' value='Uncheck all' onClick='uncheckAll(\"genotypeIdTwo\");return false;' title='Uncheck all'>
				</th>
			";

foreach (@headerRows)
{
	$genotypes .= "<th><b>" . uc $_ . "</b></th>";
}

my @relatedDataset;
my $relatedDatasetRecord;
my $datasetList=$dbh->prepare("SELECT * FROM matrix WHERE z = ? AND container LIKE 'dataset'");
$datasetList->execute($dart[6]);
while (my @datasetList = $datasetList->fetchrow_array())
{
	$genotypes .= "<th><b>$datasetType{$datasetList[3]} Dataset: $datasetList[2]</b></th>";
	push @relatedDataset, $datasetList[0];
	my $allRecord=$dbh->prepare("SELECT * FROM matrix WHERE container Like 'record' AND z = ?");
	$allRecord->execute($datasetList[0]);
	while (my @allRecord = $allRecord->fetchrow_array())
	{
	 	my $recordDetails = decode_json $allRecord[8];
		my $viewDetails = '';
		my $viewDetailsShort = '';
		my $lineNumber = 0;
		for (sort {$a <=> $b} keys %$recordDetails)
		{
			$recordDetails->{$_}->{'field'} = '' unless ($recordDetails->{$_}->{'field'});
			$recordDetails->{$_}->{'value'} = '' unless ($recordDetails->{$_}->{'value'});
			$recordDetails->{$_}->{'value'} = escapeHTML($recordDetails->{$_}->{'value'});
			if($recordDetails->{$_}->{'value'} =~ /\.(jpg|jpeg|png|tif|tiff)$/i)
			{
				next;
			}
			else
			{
				$viewDetails .= ($recordDetails->{$_}->{'value'} =~ /:\/\//) ? "" : "$recordDetails->{$_}->{'field'}: $recordDetails->{$_}->{'value'}\n";
				unless ($lineNumber > 1)
				{
					$viewDetailsShort .= ($recordDetails->{$_}->{'value'} =~ /:\/\//) ? "" : "$recordDetails->{$_}->{'field'}: $recordDetails->{$_}->{'value'}<br>";
					$lineNumber++
				}
			}
		}
		
		if (exists $relatedDatasetRecord->{$allRecord[2]}->{$datasetList[0]})
		{
			$relatedDatasetRecord->{$allRecord[2]}->{$datasetList[0]} .= "<br><a id='moreId$allRecord[0]$$' onclick='openDialog(\"itemView.cgi?itemId=$allRecord[0]\")' title='$viewDetails'>$viewDetailsShort</a>";
		}
		else
		{
			$relatedDatasetRecord->{$allRecord[2]}->{$datasetList[0]} = "<a id='moreId$allRecord[0]$$' onclick='openDialog(\"itemView.cgi?itemId=$allRecord[0]\")' title='$viewDetails'>$viewDetailsShort</a>";
		}
	}
}

$genotypes .= "</tr>
	</thead>
	<tbody>";

my $getGenotypes = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'dartGenotype' AND z = ?");
$getGenotypes->execute($dartId);
while (my @getGenotypes = $getGenotypes->fetchrow_array())
{

	my $genotypeDetails = decode_json $getGenotypes[8];
	foreach (@headerRows)
	{
		$genotypeDetails->{$_} = '' unless (exists $genotypeDetails->{$_});
	}
	my $descriptionLength = length ($genotypeDetails->{$headerRows[6]});
	if ($descriptionLength > 35)
	{
		$genotypeDetails->{$headerRows[6]} = "<a title='$genotypeDetails->{$headerRows[6]}'>". substr($genotypeDetails->{$headerRows[6]},0,25). "...". substr($genotypeDetails->{$headerRows[6]},-5). "</a>";
	}

	$genotypes .= "<tr>";
	foreach (@headerRows)
	{
		$genotypes .= ($_ eq 'genotype name') ? "<td style='text-align:center;'><input type='checkbox' id='dartOne$getGenotypes[0]$$' name='genotypeIdOne' value='$getGenotypes[0]'></td><td style='text-align:center;'><input type='checkbox' id='dartTwo$getGenotypes[0]$$' name='genotypeIdTwo' value='$getGenotypes[0]'></td><td><a onclick='closeDialog();openDialog(\"genotypeView.cgi?genotypeId=$getGenotypes[0]\")' title='Id: $getGenotypes[0]'>$genotypeDetails->{$_}</a></td>" : "<td>$genotypeDetails->{$_}</td>";
	}

	foreach (@relatedDataset)
	{
		$genotypes .= (exists $relatedDatasetRecord->{$getGenotypes[2]}->{$_}) ? "<td style='text-align:center;'>$relatedDatasetRecord->{$getGenotypes[2]}->{$_}</td>" : "<td></td>";
		
	}
	$genotypes .= "</tr>";
}
$genotypes .= "</tbody></table></form>\n";
my $button = "<ul id='dartInfoMenu$$' style='margin-top: .3em; width: 250px;'>
		<li><a><span class='ui-icon ui-icon-triangle-1-e'></span><b>SNP Analysis Tools</b></a>
			<ul style='z-index: 1000;white-space: nowrap;'>
				<li><a onclick='openDialogForm(\"genotypeComparison.cgi\",\"dartList$$\")'>Genotype Comparison</a></li>
				<li><a onclick='openDialogForm(\"genotypeDistance.cgi\",\"dartList$$\")'>Compute Genetic Distance</a></li>
				<li><a onclick='openDialogForm(\"genotypeTree.cgi\",\"dartList$$\")'>Make Phylogenetic Tree</a></li>
				<li><a onclick='openDialogForm(\"genotypeDownload.cgi\",\"dartList$$\")'>Download Genotype Data</a></li>
			</ul>
		</li>
	</ul>";

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$dartId/$dartId/g;
$html =~ s/\$dartName/$dart[2]/g;
$html =~ s/\$snpNumber/$dart[3]/g;
$html =~ s/\$genotypeNumber/$dart[4]/g;
$dart[8] = escapeHTML($dart[8]);
$dart[8] =~ s/\n/<br>/g;
$html =~ s/\$dartDescription/$dart[8]/g;
$html =~ s/\$dartCreator/$dart[9]/g;
$html =~ s/\$dartCreationDate/$dart[10]/g;
$html =~ s/\$genotypes/$genotypes/g;
$html =~ s/\$button/$button/g;
$html =~ s/\$\$/$$/g;
print header;
print $html;

__DATA__
<h3><a onclick='openDialog("dartView.cgi?dartId=$dartId")' title='View'>$dartName</a> <sub class='ui-state-disabled'>$snpNumber SNPs and $genotypeNumber genotypes</sub></h3>
<b>$dartDescription</b><br>
$button
$genotypes
<script>
$( "#dart$$" ).dataTable({
	"scrollY": "400px",
	"scrollCollapse": true,
	"paging": false
});
$( "#dartInfoMenu$$" ).menu();
</script>