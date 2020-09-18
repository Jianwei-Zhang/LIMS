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

my @vectorId = param ('itemId');
my $vectors;
for(@vectorId)
{
	my $vector = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$vector->execute($_);
	my @vector=$vector->fetchrow_array();
	$vectors->{$vector[0]} = $vector[2];
}
my $vectorList = 'Please check vectors first!';
$vectorList = checkbox_group(-name=>'vectors',
-values=>[sort keys %{$vectors}],
-default=>[sort keys %{$vectors}],
-labels=>\%{$vectors},
-columns=>3) if (keys %{$vectors});
my $vectorNameList = '';
$vectorNameList = scrolling_list(-name=>'vectorName',
-values=>[sort keys %{$vectors}],
-size=>1,
-labels=>\%{$vectors}) if (keys %{$vectors});

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$vectorList/$vectorList/g;
$html =~ s/\$vectorNameList/$vectorNameList/g;

print header;
print $html;

__DATA__
<form id="mergeVector" name="mergeVector" action="vectorMerge.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<h3 class='ui-state-error-text'>Are you sure to merge the below vectors?</h3>
	<table>
	<tr><td style='text-align:right'><label for="mergeVectorList"><b>Vectors</b></label></td><td>$vectorList</td></tr>
	<tr><td style='text-align:right'><label for="mergeVectorList"><b>Unified name</b></label></td><td>$vectorNameList</td></tr>
	</table>
</form>
<script>
$('#dialog').dialog("option", "title", "Merge Vector");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Merge", click: function() { submitForm('mergeVector'); } }, { text: "Cancel", click: function() { closeDialog(); } } ] );
</script>