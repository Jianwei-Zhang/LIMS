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
my $cookieAsbProject = cookie('asbProject') || '';

my $asbProjects = '';
my $asbProject=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'asbProject'");# ORDER BY name
$asbProject->execute();
if ($asbProject->rows > 0)
{
	my $asbProjectResult;
	while (my @asbProject = $asbProject->fetchrow_array())
	{
		@{$asbProjectResult->{$asbProject[2]}} = @asbProject;
	}
	foreach (sort {uc ($a) cmp uc($b)} keys %$asbProjectResult)
	{
		my @asbProject = @{$asbProjectResult->{$_}};
		$asbProject[2] = "AsbProject: Name N/A, please edit!" unless($asbProject[2]);
		$asbProjects .= "<li><a href='asbProject.cgi?asbProjectId=$asbProject[0]'>$asbProject[2]</a></li>\n";
		$active = $activeDetector if ($cookieAsbProject eq $asbProject[0]);
		$activeDetector++;
	}
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$asbProjects/$asbProjects/;
$html =~ s/\$active/$active/;
$html =~ s/\$commoncfg->{HTDOCS}/$commoncfg->{HTDOCS}/g;

print header(-cookie=>cookie(-name=>'menu',-value=>3));
print $html;

__DATA__
<div id="asbProjects">
	<ul>
		$asbProjects
		<li><a href='asbProjectNew.cgi'>New GPM Project</a></li>
	</ul>
</div>
<script>
loadingHide();
$( "#asbProjects" ).tabs({
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
$( "#asbProjects li" ).removeClass( "ui-corner-top" ).addClass( "ui-corner-left" );
</script>