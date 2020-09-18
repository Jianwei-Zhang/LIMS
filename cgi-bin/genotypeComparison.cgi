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
my @headerRows = ("genotype name");
#my @headerRows = ("genotype name","Order number","DArT plate barcode","client plate barcode","well row position","well column position","sample comments");

print header;
my $genotypeIdOneNumber = scalar (@genotypeIdOne);
my $genotypeIdTwoNumber = scalar (@genotypeIdTwo);

if($genotypeIdOneNumber +  $genotypeIdTwoNumber < 2)
{
	print 'Select at least 2 genotypes first!';
	exit;
}
my $genotypeDetails;
my $dartId;
foreach my $genotypeId (@genotypeIdOne,@genotypeIdTwo)
{
	my $getGenotype = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$getGenotype->execute($genotypeId);
	my @getGenotype = $getGenotype->fetchrow_array();
	$genotypeDetails->{$genotypeId} = decode_json $getGenotype[8];
	$dartId = $getGenotype[6];
}
open (TAB,">$commoncfg->{TMPDIR}/genotypeComparison.$$.txt") or die "can't open file: $commoncfg->{TMPDIR}/genotypeComparison.$$.txt";
my $genotypes = "";
my $firstLine = '';
my $secondLine = '';
my $genotypeTabLine = "";
foreach my $headerRow (@headerRows)
{
	$genotypes .= "<tr><td></td><td style='white-space: nowrap;'><b>$headerRow</b></td><td>";
	foreach my $genotypeId (@genotypeIdOne)
	{
		$genotypeDetails->{$genotypeId}->{$headerRow} = '' unless (exists $genotypeDetails->{$genotypeId}->{$headerRow});
		$genotypes .= "<sup>$genotypeDetails->{$genotypeId}->{$headerRow}</sup> ";	
		$firstLine .= "$genotypeDetails->{$genotypeId}->{$headerRow} ";
	}
	$genotypes .= "</td><td>" if ($genotypeIdOneNumber && $genotypeIdTwoNumber);
	$firstLine .= "\t" if ($genotypeIdOneNumber && $genotypeIdTwoNumber);
	foreach my $genotypeId (@genotypeIdTwo)
	{
		$genotypeDetails->{$genotypeId}->{$headerRow} = '' unless (exists $genotypeDetails->{$genotypeId}->{$headerRow});
		$genotypes .= "<sup>$genotypeDetails->{$genotypeId}->{$headerRow}</sup> ";	
		$firstLine .= "$genotypeDetails->{$genotypeId}->{$headerRow} ";
	}
	$genotypes .= "</td></tr>";
	$firstLine .= "\n";
}
$genotypes .= ($genotypeIdOneNumber) ? ($genotypeIdTwoNumber) ? "<tr><td><b>No.</b></td><td><b>SNP</b></td><td><b>G1</b></td><td><b>G2</b></td></tr>" : "<tr><td><b>No.</b></td><td><b>SNP</b></td><td><b>G1</b></td></tr>" : ($genotypeIdTwoNumber) ? "<tr><td><b>No.</b></td><td><b>SNP</b></td><td><b>G2</b></td></tr>" : "<tr><td><b>No.</b></td><td><b>SNP</b></td></tr>";
$secondLine .= ($genotypeIdOneNumber) ? ($genotypeIdTwoNumber) ? "G1\tG2\n" : "G1\n" : ($genotypeIdTwoNumber) ? "G2\n" : "\n";

my $snpNumber = 0;
my $commonGenotype = 0;
my $getSNPs = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'dartSNP' AND z = ?");
$getSNPs->execute($dartId);
while (my @getSNPs = $getSNPs->fetchrow_array())
{
	my $snpDetails = decode_json $getSNPs[8];
	my $snpDetailsLine = '';
	my $snpTabLine = '';
	for (sort {$a cmp $b} keys %$snpDetails)
	{
		$snpDetailsLine .= "$_: $snpDetails->{$_}\n";
		$snpTabLine = ($snpTabLine) ? "$snpDetails->{$_}\t$snpTabLine" : "$snpDetails->{$_}";
		unless($snpNumber > 0)
		{
			$firstLine = "*\t$firstLine";
			$secondLine = "$_\t$secondLine";
		}
	}
	my $lastSnpOne = '';
	my $noDiffOne = 1;
	my $iOne = 0;
	foreach my $genotypeId (@genotypeIdOne)
	{
		$noDiffOne = 0 and last if ($genotypeDetails->{$genotypeId}->{$getSNPs[0]} eq "NN");
		if ($iOne > 0)
		{
			$noDiffOne = 0 and last if ($genotypeDetails->{$genotypeId}->{$getSNPs[0]} ne $lastSnpOne);
		}
		$lastSnpOne = $genotypeDetails->{$genotypeId}->{$getSNPs[0]};
		$iOne++;
	}

	my $lastSnpTwo = '';
	my $noDiffTwo = 1;
	my $iTwo = 0;
	foreach my $genotypeId (@genotypeIdTwo)
	{
		$noDiffTwo = 0 and last if ($genotypeDetails->{$genotypeId}->{$getSNPs[0]} eq "NN");
		$noDiffTwo = 0 and last if ($genotypeDetails->{$genotypeId}->{$getSNPs[0]} eq $lastSnpOne );
		if ($iTwo > 0)
		{
			$noDiffTwo = 0 and last if ($genotypeDetails->{$genotypeId}->{$getSNPs[0]} ne $lastSnpTwo);
		}
		$lastSnpTwo = $genotypeDetails->{$genotypeId}->{$getSNPs[0]};
		$iTwo++;
	}

	if ($genotypeIdOneNumber && $genotypeIdTwoNumber)
	{
		if ($noDiffOne && $noDiffTwo)
		{
			$commonGenotype++;
			$genotypes .= "<tr><td>$commonGenotype</td><td style='white-space: nowrap;' title='$snpDetailsLine'>$getSNPs[2]</td><td>$lastSnpOne</td><td>$lastSnpTwo</td></tr>";
			$genotypeTabLine .= "$snpTabLine\t$lastSnpOne\t$lastSnpTwo\n";
		}
	}
	else
	{
		if ($genotypeIdOneNumber)
		{
			if ($noDiffOne)
			{
				$commonGenotype++;
				$genotypes .= "<tr><td>$commonGenotype</td><td style='white-space: nowrap;' title='$snpDetailsLine'>$getSNPs[2]</td><td>$lastSnpOne</td></tr>";
				$genotypeTabLine .= "$snpTabLine\t$lastSnpOne\n";
			}
		}
		if ( $genotypeIdTwoNumber)
		{
			if ($noDiffTwo)
			{
				$commonGenotype++;
				$genotypes .= "<tr><td>$commonGenotype</td><td style='white-space: nowrap;' title='$snpDetailsLine'>$getSNPs[2]</td><td>$lastSnpTwo</td></tr>";
				$genotypeTabLine .= "$snpTabLine\t$lastSnpTwo\n";
			}
		}
	}
	$snpNumber++;
}
print TAB "$firstLine";
print TAB "$secondLine";
print TAB "$genotypeTabLine";

unless ($commonGenotype)
{
	$genotypes .= ($genotypeIdOneNumber) ? ($genotypeIdTwoNumber) ? "<tr><td colspan='4' class='ui-state-error ui-corner-all' style='padding: .7em;'>No SNP Found.</td></tr>" : "<tr><td colspan='3' class='ui-state-error ui-corner-all' style='padding: .7em;'>No SNP Found.</td></tr>" : ($genotypeIdTwoNumber) ? "<tr><td colspan='3' class='ui-state-error ui-corner-all' style='padding: .7em;'>No SNP Found.</td></tr>" : "<tr><td colspan='2' class='ui-state-error ui-corner-all' style='padding: .7em;'>No SNP Found.</td></tr>";
	print TAB "No SNP Found.\n";
}
close (TAB);
`gzip -f $commoncfg->{TMPDIR}/genotypeComparison.$$.txt`;

open (COMPARISON,">$commoncfg->{TMPDIR}/genotypeComparison.$$.html") or die "can't open file: $commoncfg->{TMPDIR}/genotypeComparison.$$.html";
print COMPARISON "<table width='100%'>\n$genotypes\n</table>";
close (COMPARISON);
`gzip -f $commoncfg->{TMPDIR}/genotypeComparison.$$.html`;

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$\$/$$/g;
$html =~ s/\$genotypes/$genotypes/g;
$html =~ s/\$commoncfg->{HTDOCS}/$commoncfg->{HTDOCS}/g;
$html =~ s/\$commoncfg->{TMPURL}/$commoncfg->{TMPURL}/g;

print $html;

__DATA__
<div id="genotypeComparison$$" name="genotypeComparison$$">
	<table width='100%'>
	$genotypes
	</table>
</div>
<script>
buttonInit();
$('#dialog').dialog("option", "title", "Genotype Comparison");
$( "#dialog" ).dialog( "option", "buttons", [ { text: "Text", click: function() { location.href='$commoncfg->{TMPURL}/genotypeComparison.$$.txt.gz'; } }, { text: "HTML", click: function() { location.href='$commoncfg->{TMPURL}/genotypeComparison.$$.html.gz'; } }, { text: "Print", click: function() {printDiv('genotypeComparison$$'); } }, { text: "OK", click: function() {closeDialog(); } } ] );
</script>
