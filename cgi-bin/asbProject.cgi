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

my $active = 0;
my $activeDetector = 0;
my $cookieAsbProjectContent = cookie('asbProjectContent') || '';

my $asbProjectId = param ('asbProjectId') || '';
my $button;
my $asbProjectContent = '';
if ($asbProjectId)
{
	my $targetInAsbProject=$dbh->prepare("SELECT matrix.* FROM link,matrix WHERE link.type LIKE 'asbProject' AND link.child = matrix.id AND link.parent = $asbProjectId ORDER BY matrix.name");
	$targetInAsbProject->execute();
	while (my @targetInAsbProject = $targetInAsbProject->fetchrow_array())
	{
		$asbProjectContent = "<div id='asbProjectContentInAsbProject$asbProjectId$$'><ul>\n" unless($asbProjectContent);
		$asbProjectContent .= "<li><a href='assembly.cgi?targetId=$targetInAsbProject[0]' title ='$targetInAsbProject[2]'>$targetInAsbProject[2]</a></li>\n";
		$active = $activeDetector if ($cookieAsbProjectContent eq $targetInAsbProject[0]);
		$activeDetector++;
	}
	my $asbProject=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$asbProject->execute($asbProjectId);
	my @asbProject = $asbProject->fetchrow_array();
	$asbProject[8] = escapeHTML($asbProject[8]);
	$asbProject[8] =~ s/\n/<br>/g;
	$button = "<div class='ui-state-highlight ui-corner-all' style='padding: 0 .7em;'>
		<button style='float: right; margin-top: .3em; margin-right: .3em;' onclick='openDialog(\"genomeNew.cgi?asbProjectId=$asbProjectId\")'>New Genome for GPM</button>
		<div style='position: relative;'><h2><a id='asbProject$asbProjectId$$' onmouseover='editIconShow(\"asbProject$asbProjectId$$\")' onmouseout='editIconHide(\"asbProject$asbProjectId$$\")' onclick='openDialog(\"asbProjectEdit.cgi?asbProjectId=$asbProjectId\")' title='Edit/Delete'>$asbProject[2]</a></h2></div>";
	$button .= "<p>$asbProject[8]</p>" if ($asbProject[8]);
	unless($asbProjectContent)
	{
		$button .= "<p class='ui-state-error ui-corner-all' style='padding: .7em;'><span class='ui-icon ui-icon-alert' style='float: left; margin-right: .3em;'></span>
			<strong>No content in this Assembly Project.</strong></p>";
	}
	$button .= "</div>";
	print header(-cookie=>cookie(-name=>'asbProject',-value=>$asbProjectId));
}
else
{
	print header(-type=>'text/html',-status=>'402 Invalid operation');
	exit;
}

$asbProjectContent .= "</ul></div>\n" if($asbProjectContent);

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$asbProjectId/$asbProjectId/g;
$html =~ s/\$\$/$$/g;
$html =~ s/\$button/$button/g;
$html =~ s/\$asbProjectContent/$asbProjectContent/g;
$html =~ s/\$commoncfg->{HTDOCS}/$commoncfg->{HTDOCS}/g;
$html =~ s/\$active/$active/g;

print $html;

__DATA__
$button
$asbProjectContent
<script>
buttonInit();
loadingHide();
$( "#asbProjectContentInAsbProject$asbProjectId$$" ).tabs({
	// loading spinner
	beforeLoad: function(event, ui) {
		ui.panel.html('<img src="$commoncfg->{HTDOCS}/css/images/loading.gif" style="vertical-align:middle;"> Loading...');
	},
	create: function (e, ui) {
		loadingShow();
	},
	activate: function (e, ui) {
		loadingShow();
	},
	active: $active
}).addClass( "ui-tabs-vertical ui-helper-clearfix" );
$( "#asbProjectContentInAsbProject$asbProjectId$$ li" ).removeClass( "ui-corner-top" ).addClass( "ui-corner-left" );
$( "#libraryInfoMenu$asbProjectId$$" ).menu();
</script>