#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
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

my @genomeId = param ('itemId');
my $genomes;
for(@genomeId)
{
	my $genome = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$genome->execute($_);
	my @genome=$genome->fetchrow_array();
	$genomes->{$genome[0]} = $genome[2];
}
my $genomeList = 'Please check genomes first!';
$genomeList = checkbox_group(-name=>'genomes',
-values=>[sort keys %{$genomes}],
-default=>[sort keys %{$genomes}],
-labels=>\%{$genomes},
-columns=>3) if (keys %{$genomes});
my $genomeNameList = '';
$genomeNameList = scrolling_list(-name=>'mergedGenomeId',
-values=>[sort keys %{$genomes}],
-size=>1,
-labels=>\%{$genomes}) if (keys %{$genomes});

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$genomeList/$genomeList/g;
$html =~ s/\$genomeNameList/$genomeNameList/g;

print header;
print $html;

__DATA__
<form id="mergeGenome" name="mergeGenome" action="genomeMerge.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<h3 class='ui-state-error-text'>Are you sure to merge the below genomes?</h3>
	<table>
	<tr><td style='text-align:right'><label for="mergeGenomeList"><b>Genomes</b></label></td><td>$genomeList</td></tr>
	<tr><td style='text-align:right'><label for="mergeGenomeList"><b>Unified name</b></label></td><td>$genomeNameList</td></tr>
	</table>
	<h4 class='ui-state-error-text'>Caution: Merging genomes might cause some unpredictable broken links in the database. It's recommended to upload a new geneme if you are not familiar with this operation.</h4>
</form>
<script>
$('#dialog').dialog("option", "title", "Merge Genome");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Merge", click: function() { submitForm('mergeGenome'); } }, { text: "Cancel", click: function() { closeDialog(); } } ] );
</script>