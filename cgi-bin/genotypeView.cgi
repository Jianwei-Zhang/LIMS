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

my $genotypeId = param ('genotypeId') || '';
my @headerRows = ("genotype name","Order number","DArT plate barcode","client plate barcode","well row position","well column position","sample comments");

print header;
if ($genotypeId)
{
	my $genotypes = "";
	my $getGenotype = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$getGenotype->execute($genotypeId);
	if($getGenotype->rows < 1)
	{
		print 'No valid Genotype Found!';
		exit;
	}
	else
	{
		my @getGenotype = $getGenotype->fetchrow_array();
		my $genotypeDetails = decode_json $getGenotype[8];
		my $dartToGenebank=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$dartToGenebank->execute($getGenotype[6]);
		my @dartToGenebank = $dartToGenebank->fetchrow_array();

		my $extraItem=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'record' AND name LIKE ?");
		$extraItem->execute($getGenotype[2]);
		while (my @extraItem = $extraItem->fetchrow_array())
		{
			my $extraItemDetails = decode_json $extraItem[8];
			my $extraParent=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
			$extraParent->execute($extraItem[6]);
			my @extraParent = $extraParent->fetchrow_array();
			if ($extraParent[6] == $dartToGenebank[6])
			{
				$genotypes .= "<tr><td style='white-space: nowrap;' colspan='2' class='ui-state-highlight ui-corner-all'><h3>$datasetType{$extraParent[3]} Dataset: $extraParent[2]</h3></td></tr>";
				for (sort {$a <=> $b} keys %$extraItemDetails)
				{
					$extraItemDetails->{$_}->{'field'} = '' unless ($extraItemDetails->{$_}->{'field'});
					$extraItemDetails->{$_}->{'value'} = '' unless ($extraItemDetails->{$_}->{'value'});
					$extraItemDetails->{$_}->{'value'} = escapeHTML($extraItemDetails->{$_}->{'value'});
					$extraItemDetails->{$_}->{'value'} =~ s/\n/<br>/g;
					if($extraItemDetails->{$_}->{'value'} =~ /\.(jpg|jpeg|png|tif|tiff)$/i)
					{
						$genotypes .= "<tr><td style='text-align:right;width:200px;'><b>$extraItemDetails->{$_}->{'field'}</b></td><td>";
						for (split "\;", $extraItemDetails->{$_}->{'value'})
						{
							$_ =~ s/^\s+|\s+$//g;
							$genotypes .= "<img src='$commoncfg->{HTDOCS}/data/images/$_'/>";
						}
						$genotypes .= "</td></tr>";
					}
					else
					{
						$genotypes .= ($extraItemDetails->{$_}->{'value'} =~ /:\/\//) ? "<tr><td style='text-align:right;width:200px;'><b>$extraItemDetails->{$_}->{'field'}</b></td>
										<td><a href='$extraItemDetails->{$_}->{'value'}' target='_blank'>$extraItemDetails->{$_}->{'value'}</a></td>
									</tr>" :
									"<tr><td style='text-align:right;width:200px;'><b>$extraItemDetails->{$_}->{'field'}</b></td>
										<td>$extraItemDetails->{$_}->{'value'}</td>
									</tr>";
					}
				}
			}
		}

		$genotypes .= "<tr><td style='white-space: nowrap;' colspan='2' class='ui-state-highlight ui-corner-all'><h3>Genotype details</h3></td></tr>";

# 		foreach (@headerRows)
# 		{
# 			$genotypeDetails->{$_} = '' unless (exists $genotypeDetails->{$_});
# 			$genotypes .= "<tr><td style='white-space: nowrap;'><b>$_</b></td><td>$genotypeDetails->{$_}</td></tr>";
# 		}


		my $getSNPs = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'dartSNP' AND z = ?");
		$getSNPs->execute($getGenotype[6]);
		while (my @getSNPs = $getSNPs->fetchrow_array())
		{
			my $snpDetails = decode_json $getSNPs[8];
			my $snpDetailsLine = '';
# 			for (keys %$snpDetails)
# 			{
# 				$snpDetailsLine .= "$_: $snpDetails->{$_}\n";
# 			}
			$genotypes .= "<tr><td style='white-space: nowrap;' title='$snpDetailsLine'><b>$getSNPs[2]</b></td><td>$genotypeDetails->{$getSNPs[0]}</td></tr>";
		}

		undef $/;# enable slurp mode
		my $html = <DATA>;
		$html =~ s/\$\$/$$/g;
		$html =~ s/\$genotypeId/$genotypeId/g;		
		$html =~ s/\$genotypeName/$getGenotype[2]/g;
		$html =~ s/\$genotypes/$genotypes/g;
		$html =~ s/\$commoncfg->{HTDOCS}/$commoncfg->{HTDOCS}/g;

		print $html;
	}
}
else
{
	print 'No valid Genotype Found!';
	exit;
}

__DATA__
<div id="viewGenotype$$" name="viewGenotype$$">
	<table width='100%'>
	$genotypes
	</table>
</div>
<script>
buttonInit();
$('#dialog').dialog("option", "title", "View Genotype $genotypeName");
$( "#dialog" ).dialog( "option", "buttons", [ { text: "Print", click: function() {printDiv('viewGenotype$$'); } }, { text: "OK", click: function() {closeDialog(); } } ] );
</script>
