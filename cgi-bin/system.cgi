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

my $button = "<div class='ui-state-highlight ui-corner-all' style='padding: 0 .7em;'>
	<h2>System tools</h2>
	<button onclick='refresh(\"general\")'>Refresh</button>
	<button class='ui-state-error-text' onclick='deleteItem(0)' title='Delete Orphan Records'>Delete Orphan Records</button>
	<button class='ui-state-error-text' onclick='deleteItem(0,\"delTempFiles\")' title='Delete Cached Files'>Delete Cached Files</button>
	</div>\n";

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$button/$button/g;
$html =~ s/\$\$/$$/g;

print header(-cookie=>cookie(-name=>'general',-value=>7));
print $html;

__DATA__
$button
<script>
buttonInit();
loadingHide();
</script>