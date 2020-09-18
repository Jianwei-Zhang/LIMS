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

my @genotypeIdOne = param ('genotypeIdOne');
my @genotypeIdTwo = param ('genotypeIdTwo');
push @genotypeIdOne, @genotypeIdTwo;
my %genotypeIdAll = map { $_, 1 } @genotypeIdOne;
       # or a hash slice: @genotypeId{ @genotypeIdOne } = ();
       # or a foreach: $genotypeId{$_} = 1 foreach ( @genotypeIdOne );
my @genotypeIdAll = keys %genotypeIdAll;
my @headerRows = ("genotype name");
#my @headerRows = ("genotype name","Order number","DArT plate barcode","client plate barcode","well row position","well column position","sample comments");

print header;
my $genotypeNumber = scalar (@genotypeIdAll);
if($genotypeNumber < 2)
{
	print 'Select at least 2 genotypes first!';
	exit;
}

if($genotypeNumber > 100)
{
	print 'Too many (>100) genotypes selected!';
	exit;
}

my $genotypeName;
my $genotypeDetails;
my $dartId;
foreach my $genotypeId (@genotypeIdAll)
{
	my $getGenotype = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$getGenotype->execute($genotypeId);
	my @getGenotype = $getGenotype->fetchrow_array();
	$genotypeDetails->{$genotypeId} = decode_json $getGenotype[8];
	$dartId = $getGenotype[6];
	if (exists $genotypeName->{$genotypeDetails->{$genotypeId}->{'genotype name'}})
	{
		$genotypeName->{$genotypeDetails->{$genotypeId}->{'genotype name'}}++;
	}
	else
	{
		$genotypeName->{$genotypeDetails->{$genotypeId}->{'genotype name'}} = 1;
	}
}

my $distancePair;
my $getSNPs = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'dartSNP' AND z = ?");
$getSNPs->execute($dartId);
my $totalSNP = $getSNPs->rows;
while (my @getSNPs = $getSNPs->fetchrow_array())
{
	foreach my $genotypeIdA (@genotypeIdAll)
	{
		foreach my $genotypeIdB (@genotypeIdAll)
		{
			next if ($genotypeIdA eq $genotypeIdB);
			if ($genotypeDetails->{$genotypeIdA}->{$getSNPs[0]} eq "NN" || $genotypeDetails->{$genotypeIdB}->{$getSNPs[0]} eq "NN")
			{
				if (exists $distancePair->{$genotypeIdA}->{$genotypeIdB})
				{
					$distancePair->{$genotypeIdA}->{$genotypeIdB}++;
				}
				else
				{
					$distancePair->{$genotypeIdA}->{$genotypeIdB} = 1;
				}
			}
			else
			{
				if ($genotypeDetails->{$genotypeIdA}->{$getSNPs[0]} ne $genotypeDetails->{$genotypeIdB}->{$getSNPs[0]})
				{
					if (exists $distancePair->{$genotypeIdA}->{$genotypeIdB})
					{
						$distancePair->{$genotypeIdA}->{$genotypeIdB}++;
					}
					else
					{
						$distancePair->{$genotypeIdA}->{$genotypeIdB} = 1;
					}
				}
			}
		}
	}
}

open (DISTANCE,">$commoncfg->{TMPDIR}/genotypeDistance.$$.meg") or die "can't open file: $commoncfg->{TMPDIR}/genotypeDistance.$$.meg";
my $tableHead = "<thead><tr><th><b>$genotypeNumber genotypes</b></th>";
my $tableFoot = "<tfoot><tr><th><b>$totalSNP SNPs</b></th>";
print DISTANCE "#mega
!TITLE  Genetic distance data from $genotypeNumber genotypes;
!Format DataType=distance DataFormat=LowerLeft NTaxa=$genotypeNumber;
!Description
     Number of polymoprhic loci used = $totalSNP
;
";

foreach my $genotypeId (@genotypeIdAll)
{
	if ($genotypeName->{$genotypeDetails->{$genotypeId}->{'genotype name'}} > 1)
	{
		$tableHead .= "<th>$genotypeDetails->{$genotypeId}->{'genotype name'}<sub>.$genotypeId</sub></th>";
		$tableFoot .= "<th>$genotypeDetails->{$genotypeId}->{'genotype name'}<sub>.$genotypeId</sub></th>";
		print DISTANCE "#$genotypeDetails->{$genotypeId}->{'genotype name'}.$genotypeId\n";
	}
	else
	{
		$tableHead .= "<th>$genotypeDetails->{$genotypeId}->{'genotype name'}</th>";
		$tableFoot .= "<th>$genotypeDetails->{$genotypeId}->{'genotype name'}</th>";
		print DISTANCE "#$genotypeDetails->{$genotypeId}->{'genotype name'}\n";
	}
}
$tableHead .= "</tr></thead>";
$tableFoot .= "</tr></tfoot>";

my $tableBody = "<tbody>";
my $stop = 1;
foreach my $genotypeIdA (@genotypeIdAll)
{
	if ($genotypeName->{$genotypeDetails->{$genotypeIdA}->{'genotype name'}} > 1)
	{
		$tableBody .= "<tr><td>$genotypeDetails->{$genotypeIdA}->{'genotype name'}<sub>.$genotypeIdA</sub></td>";
	}
	else
	{
		$tableBody .= "<tr><td>$genotypeDetails->{$genotypeIdA}->{'genotype name'}</td>";
	}
	$stop = 0;
	foreach my $genotypeIdB (@genotypeIdAll)
	{
		if ($genotypeIdA eq $genotypeIdB)
		{
			$tableBody .= "<td>0</td>";	
			$stop = 1;
		}
		else
		{
			my $distance = (exists $distancePair->{$genotypeIdA}->{$genotypeIdB}) ? $distancePair->{$genotypeIdA}->{$genotypeIdB}/$totalSNP : 1;
			$tableBody .= (!$stop) ? "<td>$distance</td>" : "<td></td>";
			print DISTANCE " $distance" if (!$stop);
		}
	}
	$tableBody .= "</tr>";
	print DISTANCE "\n";
}
$tableBody .= "</tbody>";
my $distances = $tableHead.$tableBody.$tableFoot;
close (DISTANCE);
`gzip -f $commoncfg->{TMPDIR}/genotypeDistance.$$.meg`;

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$\$/$$/g;
$html =~ s/\$distances/$distances/g;
$html =~ s/\$commoncfg->{HTDOCS}/$commoncfg->{HTDOCS}/g;
$html =~ s/\$commoncfg->{TMPURL}/$commoncfg->{TMPURL}/g;

print $html;

__DATA__
<div id="genotypeDistance$$" name="genotypeDistance$$">
	<table width='100%' id='genotypeDistanceMatrix$$' class='display'>
	$distances
	</table>
</div>
<script>
buttonInit();
$('#dialog').dialog("option", "title", "Genotype Distance");
$( "#dialog" ).dialog( "option", "buttons", [ { text: "Download", click: function() { location.href='$commoncfg->{TMPURL}/genotypeDistance.$$.meg.gz'; } }, { text: "Print", click: function() {printDiv('genotypeDistance$$'); } }, { text: "OK", click: function() {closeDialog(); } } ] );
$( "#genotypeDistanceMatrix$$" ).dataTable({
	"scrollY": "400px",
	"scrollCollapse": true,
	"searching": false,
	"paging": false,
	"ordering":  false
});
</script>
