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

my $nameForMatrix; #build name for matrix due to some duplicated accessions
my $keyForMatrix; #build key for distance matrix
my $nameMaxLength = 10;
my $keyString = "A" x $nameMaxLength;
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
	$keyForMatrix->{$genotypeId} = $keyString;
	$keyString++;
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

`rm -fr $commoncfg->{TMPDIR}/genotypeDistance.$$.in.gz $commoncfg->{TMPDIR}/genotypeDistance.$$.out.gz $commoncfg->{TMPDIR}/genotypeDistance.$$.tre.gz $commoncfg->{TMPDIR}/genotypeDistance.$$.screen`;
open (DISTANCE,">$commoncfg->{TMPDIR}/genotypeDistance.$$.in") or die "can't open file: $commoncfg->{TMPDIR}/genotypeDistance.$$.in";
print DISTANCE (" " x 4) . "$genotypeNumber\n";
foreach my $genotypeIdA (@genotypeIdAll)
{
	$nameForMatrix->{$genotypeIdA} = $genotypeDetails->{$genotypeIdA}->{'genotype name'};
	$nameForMatrix->{$genotypeIdA} =~ s/^\s+|\s+$//; 
	if ($genotypeName->{$genotypeDetails->{$genotypeIdA}->{'genotype name'}} > 1)
	{
		$nameForMatrix->{$genotypeIdA} .= ".$genotypeIdA";
	}
	print DISTANCE $keyForMatrix->{$genotypeIdA};
	foreach my $genotypeIdB (@genotypeIdAll)
	{
		my $distance = ($genotypeIdA eq $genotypeIdB) ? "0.00000" : (exists $distancePair->{$genotypeIdA}->{$genotypeIdB}) ? $distancePair->{$genotypeIdA}->{$genotypeIdB}/$totalSNP : "1.00000";
		$distance = (length ($distance) > 7) ? substr ($distance, 0, 7) : $distance . ("0" x (7- length($distance)));
		print DISTANCE " $distance";
	}
	print DISTANCE "\n";
}
close (DISTANCE);

open (SETTING,">$commoncfg->{TMPDIR}/genotypeDistance.setting") or die "can't open file: $commoncfg->{TMPDIR}/genotypeDistance.setting";
print SETTING "$commoncfg->{TMPDIR}/genotypeDistance.$$.in
F
$commoncfg->{TMPDIR}/genotypeDistance.$$.out
2
Y
F
$commoncfg->{TMPDIR}/genotypeDistance.$$.tre
";
close (SETTING);
system ("/usr/local/bin/neighbor < $commoncfg->{TMPDIR}/genotypeDistance.setting > $commoncfg->{TMPDIR}/genotypeDistance.$$.screen");
foreach my $genotypeIdC (@genotypeIdAll)
{
	`perl -p -i -e 's/$keyForMatrix->{$genotypeIdC}/$nameForMatrix->{$genotypeIdC}/gi' $commoncfg->{TMPDIR}/genotypeDistance.$$.*`;
}
open (OUTFILE,"$commoncfg->{TMPDIR}/genotypeDistance.$$.out") or die "can't open file: $commoncfg->{TMPDIR}/genotypeDistance.$$.out";
my $outfile = <OUTFILE>;
close (OUTFILE);
my $tree = '';
open (TREE,"$commoncfg->{TMPDIR}/genotypeDistance.$$.tre") or die "can't open file: $commoncfg->{TMPDIR}/genotypeDistance.$$.tre";
my $tree = <TREE>;
$tree =~ s/\s+//g;
close (TREE);

`gzip -f $commoncfg->{TMPDIR}/genotypeDistance.$$.in`;
`gzip -f $commoncfg->{TMPDIR}/genotypeDistance.$$.out`;
`gzip -f $commoncfg->{TMPDIR}/genotypeDistance.$$.tre`;

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$\$/$$/g;
$html =~ s/\$outfile/$outfile/g;
$html =~ s/\$tree/$tree/g;
$html =~ s/\$commoncfg->{HTDOCS}/$commoncfg->{HTDOCS}/g;
$html =~ s/\$commoncfg->{TMPURL}/$commoncfg->{TMPURL}/g;

print $html;

#http://bl.ocks.org/kueda/raw/1036776/

__DATA__
<div id="genotypeTree$$" name="genotypeTree$$">
<ul>
	<li><a href='#treePHYLIP$$'>PHYLIP Neighbor-joining</a></li>
	<li><a href='#phylogramTab$$'>Phylogram</a></li>
	<li><a href='#radialtreeTab$$'>Circular Dendrogram</a></li>
</ul>
<div id='treePHYLIP$$'><pre>$outfile</pre></div>
<div id='phylogramTab$$'><div id='phylogram$$'></div></div>
<div id='radialtreeTab$$'><div id='radialtree$$'></div></div>
</div>
<script>
var newick = Newick.parse("$tree");
var newickNodes = [];
function buildNewickNodes(node, callback) {
    newickNodes.push(node);
    if (node.branchset) {
    	for (var i=0; i < node.branchset.length; i++) {
            buildNewickNodes(node.branchset[i]);
        }
    }
}
buildNewickNodes(newick);
d3.phylogram.build('#phylogram$$', newick, {
    width: 500,
    height: 550
});
d3.phylogram.buildRadial('#radialtree$$', newick, {
    width: 550,
    skipLabels: false
});

buttonInit();
$( "#genotypeTree$$" ).tabs();
$('#dialog').dialog("option", "title", "Genotype Tree");
$( "#dialog" ).dialog( "option", "buttons", [ { text: "Distance Matrix", click: function() { location.href='$commoncfg->{TMPURL}/genotypeDistance.$$.in.gz'; } }, { text: "Output", click: function() { location.href='$commoncfg->{TMPURL}/genotypeDistance.$$.out.gz'; } }, { text: "Tree", click: function() { location.href='$commoncfg->{TMPURL}/genotypeDistance.$$.tre.gz'; } }, { text: "Print", click: function() {printDiv('genotypeTree$$'); } }, { text: "OK", click: function() {closeDialog(); } } ] );
</script>
