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
my $cookieFreezer = cookie('freezer') || '';

my $roomId = param ('roomId') || '';
my $button;
my $freezers = '';
if ($roomId)
{
	my $freezerInRoom=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'freezer' AND o = $roomId ORDER BY name");
	$freezerInRoom->execute();
	while (my @freezerInRoom = $freezerInRoom->fetchrow_array())
	{
		$freezers = "<div id='freezersInRoom$roomId$$'><ul>\n" unless($freezers);
		$freezers .= "<li><a href='freezer.cgi?freezerId=$freezerInRoom[0]'>Freezer: $freezerInRoom[2]</a></li>\n";
		$active = $activeDetector if ($cookieFreezer eq $freezerInRoom[0]);
		$activeDetector++;
	}
	$freezers .= "</ul></div>\n" if($freezers);
	my $engaged = $freezerInRoom->rows;
	my $room=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$room->execute($roomId);
	my @room = $room->fetchrow_array();
	$room[8] = escapeHTML($room[8]);
	$room[8] =~ s/\n/<br>/g;
	my $maxLoads= $room[4]*$room[5]*$room[6];
	my $disengaged = $maxLoads - $engaged;
	$button = "<div class='ui-state-highlight ui-corner-all' style='padding: 0 .7em;'>";
	$button .= "<div style='float: right; margin-top: .3em; margin-right: .3em;' id='progressbar$roomId$$'><div class='progress-label'>$engaged/$maxLoads loads</div></div>";	
	$button .= "<button style='z-index: 1;float: right; margin-top: .3em; margin-right: .3em;' onclick='openDialog(\"freezerNew.cgi?roomId=$roomId\")'>New Freezer</button>" if ($disengaged > 0);
	$button .= "<div style='position: relative;'><h2><a id='room$roomId$$' onmouseover='editIconShow(\"room$roomId$$\")' onmouseout='editIconHide(\"room$roomId$$\")' onclick='openDialog(\"roomEdit.cgi?roomId=$roomId\")' title='Edit/Delete Room'>$room[2]</a></h2></div>";	
	$button .= "<p>$room[8]</p>" if ($room[8]);	
	$html =~ s/\$maxLoads/$maxLoads/g;
	$html =~ s/\$engaged/$engaged/g;
	unless($freezers)
	{
		$button .= "<p class='ui-state-error ui-corner-all' style='padding: .7em;'><span class='ui-icon ui-icon-alert' style='float: left; margin-right: .3em;'></span>
			<strong>No freezer in this room.</strong></p>";
	}
	$button .= "</div>";
	print header(-cookie=>cookie(-name=>'room',-value=>$roomId));
}
else
{
	print header(-type=>'text/html',-status=>'402 Invalid operation');
	exit;
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$roomId/$roomId/g;
$html =~ s/\$\$/$$/g; #plus extra id
$html =~ s/\$button/$button/g;
$html =~ s/\$freezers/$freezers/g;
$html =~ s/\$active/$active/g;
$html =~ s/\$commoncfg->{HTDOCS}/$commoncfg->{HTDOCS}/g;

print $html;

__DATA__
$button
$freezers
<style>
  .ui-progressbar {
    width: 50%;
    position: relative;
  }
  .progress-label {
    position: absolute;
    left: 45%;
    top: 4px;
    font-weight: bold;
    text-shadow: 1px 1px 0 #fff;
  }
</style>
<script>
buttonInit();
$( "#progressbar$roomId$$" ).progressbar({
	max:$maxLoads,
	value: $engaged
});
loadingHide();
$( "#freezersInRoom$roomId$$" ).tabs({
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
$( "#freezersInRoom$roomId$$ li" ).removeClass( "ui-corner-top" ).addClass( "ui-corner-left" );
</script>