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
if(!$userId)
{
	print header;
	print <<eof;
	<h2>Error</h2>
	Please <a onclick="openDialog('login.cgi');">login</a> first!
	<script type="text/javascript">
		openDialog('login.cgi');
		loadingHide();
	</script>
eof
    exit;
}

my $commoncfg = readConfig("main.conf");
my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});

my $active = 0;
my $activeDetector = 0;
my $cookieRoom = cookie('room') || '';

my $rooms = '';
my $room=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'room'");# ORDER BY name
$room->execute();
if($room->rows > 0)
{
	my $roomResult;
	while (my @room = $room->fetchrow_array())
	{
		@{$roomResult->{$room[2]}} = @room;
	}
	foreach (sort {uc ($a) cmp uc($b)} keys %$roomResult)
	{
		my @room = @{$roomResult->{$_}};
		$room[2] = "Room Name N/A!" unless($room[2]);
		$rooms .= "<li><a href='room.cgi?roomId=$room[0]'>$room[2]</a></li>\n";
		$active = $activeDetector if ($cookieRoom eq $room[0]);
		$activeDetector++;
	}
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$rooms/$rooms/g;
$html =~ s/\$active/$active/g;
$html =~ s/\$commoncfg->{HTDOCS}/$commoncfg->{HTDOCS}/g;

print header(-cookie=>cookie(-name=>'menu',-value=>2));
print $html;

__DATA__
<div id="rooms">
	<ul>
		$rooms
		<li><a href='roomNew.cgi'>New room</a></li>
	</ul>
</div>
<script>
loadingHide();
$( "#rooms" ).tabs({
	// loading spinner
	beforeLoad: function(event, ui) {
		ui.panel.html('<img src="$commoncfg->{HTDOCS}/css/images/loading.gif" style="vertical-align:middle;"> Loading...');
	},
	create: function (e, ui) {
		loadingShow();
	},
	activate: function (e, ui) {
		loadingShow();
		$('html, body').animate({ scrollTop: 0 }, 'slow');
	},
	active: $active
}).addClass( "ui-tabs-vertical ui-helper-clearfix" );
$( "#rooms li" ).removeClass( "ui-corner-top" ).addClass( "ui-corner-left" );
</script>