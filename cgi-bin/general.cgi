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

my $active = cookie('general') || 0;;

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$active/$active/;
$html =~ s/\$commoncfg->{HTDOCS}/$commoncfg->{HTDOCS}/g;

print header(-cookie=>cookie(-name=>'menu',-value=>4));
print $html;

__DATA__
<div id="general">
	<ul>
		<li><a href="vector.cgi">Vectors</a></li>
		<li><a href="serviceAll.cgi">Services</a></li>
		<li><a href="smrtrun.cgi">SMRT Runs</a></li>
		<li><a href='job.cgi'>PacBio Assembly</a></li>
		<li><a href="genome.cgi">Genomes</a></li>
		<li><a href="dataset.cgi">Datasets</a></li>
		<li><a href="dart.cgi">DArTseq</a></li>
		<li><a href="system.cgi">System</a></li>
	</ul>
</div>
<script>
	loadingHide();
	$( "#general" ).tabs({
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
	$( "#general li" ).removeClass( "ui-corner-top" ).addClass( "ui-corner-left" );
</script>