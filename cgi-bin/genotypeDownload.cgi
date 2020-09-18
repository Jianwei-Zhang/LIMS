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
if($genotypeNumber < 1)
{
	print 'Select at least 1 genotype first!';
	exit;
}

my $readyPhrase = ($genotypeNumber > 1) ? "$genotypeNumber selected genotypes are ready for downloading." : "$genotypeNumber selected genotype is ready for downloading.";

my $genotypeDetails;
my $dartId;
foreach my $genotypeId (@genotypeIdAll)
{
	my $getGenotype = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$getGenotype->execute($genotypeId);
	my @getGenotype = $getGenotype->fetchrow_array();
	$genotypeDetails->{$genotypeId} = decode_json $getGenotype[8];
	$dartId = $getGenotype[6];
}

open (GENOTYPE,">$commoncfg->{TMPDIR}/genotype.$$.txt") or die "can't open file: $commoncfg->{TMPDIR}/genotype.$$.txt";
my $SNPlineNumber = 0;
my $getSNPs = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'dartSNP' AND z = ?");
$getSNPs->execute($dartId);
while (my @getSNPs = $getSNPs->fetchrow_array())
{
	my $snpDetails = decode_json $getSNPs[8];
	unless ($SNPlineNumber > 0)
	{
		my $firstLine = '';
		my $secondLine = '';
 		for (sort {$a cmp $b} keys %$snpDetails)
 		{
			$firstLine .= ($firstLine) ?  "\t*" : "*";
			$secondLine .= ($secondLine) ?  "\t$_" : "$_";
 		}		
		foreach my $genotypeId (@genotypeIdAll)
		{
			$firstLine .= "\t-";
			$secondLine .= "\t$genotypeDetails->{$genotypeId}->{'genotype name'}";
		}
		print GENOTYPE "$firstLine\n$secondLine\n";
	}
	my $genotypeLine = '';
 	for (sort {$a cmp $b} keys %$snpDetails)
 	{
 		$genotypeLine .= ($genotypeLine) ? "\t$snpDetails->{$_}" : "$snpDetails->{$_}";
 	}
	foreach my $genotypeId (@genotypeIdAll)
	{
		$genotypeLine .= "\t$genotypeDetails->{$genotypeId}->{$getSNPs[0]}";
	}
	print GENOTYPE "$genotypeLine\n";
	$SNPlineNumber++;
}

close (GENOTYPE);
`gzip -f $commoncfg->{TMPDIR}/genotype.$$.txt`;

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$\$/$$/g;
$html =~ s/\$readyPhrase/$readyPhrase/g;
$html =~ s/\$commoncfg->{HTDOCS}/$commoncfg->{HTDOCS}/g;
$html =~ s/\$commoncfg->{TMPURL}/$commoncfg->{TMPURL}/g;

print $html;

__DATA__
$readyPhrase
<script>
$('#dialog').dialog("option", "title", "Download Genotype");
$( "#dialog" ).dialog( "option", "buttons", [ { text: "Download", click: function() { location.href='$commoncfg->{TMPURL}/genotype.$$.txt.gz'; } }, { text: "OK", click: function() {closeDialog(); } } ] );
</script>
