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

my $fpcCtgId = param ('fpcCtgId') || '';
my $highlight = param ('highlight') || '';
my $fpcCtgDetails = '';
if ($fpcCtgId)
{
	my $fpcCtgList=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$fpcCtgList->execute($fpcCtgId);
	my @fpcCtgList = $fpcCtgList->fetchrow_array();
	$fpcCtgList[8] =~ s/\n/<br>/g;
	my $fpcClone;
	my $fpcCloneSequenced;
	my $fpcCloneMTP;
	my $fpcCloneHighlighted;
	my $fpcCloneLeftEnd;
	my $fpcCloneRightEnd;
	my $fpcCloneMaxEnd = 0;
	my $fpcCloneList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'fpcClone' AND o = ? AND MATCH (note) AGAINST (?)");
	$fpcCloneList->execute($fpcCtgList[3],$fpcCtgList[2]);
	while (my @fpcCloneList = $fpcCloneList->fetchrow_array())
	{
		$fpcClone->{$fpcCloneList[2]} = $fpcCloneList[8];
		$fpcCloneSequenced->{$fpcCloneList[2]} = $fpcCloneList[4];
		$fpcCloneMTP->{$fpcCloneList[2]} = $fpcCloneList[5];
		$fpcCloneHighlighted->{$fpcCloneList[2]} = $fpcCloneList[6];
		$fpcCloneList[8] =~ /Ends Left (\d*)/;
		$fpcCloneLeftEnd->{$fpcCloneList[2]} = $1;
		$fpcCloneList[8] =~ /Ends Right (\d*)/;
		$fpcCloneRightEnd->{$fpcCloneList[2]} = $1;
		$fpcCloneMaxEnd = $1 if ($1 > $fpcCloneMaxEnd);
	}
	my @fpcClone = sort { $fpcCloneLeftEnd->{$a} <=> $fpcCloneLeftEnd->{$b} } keys %$fpcCloneLeftEnd;

	$fpcCtgDetails .= "<table width='100%'>";
	my $col = 3;
	my $colCount=0;
	my $sequenced = 0;
    for (@fpcClone)
    {
		my $sequencedFlag = ($fpcCloneSequenced->{$_}) ? "<span class='ui-icon ui-icon-flag' style='float: left; margin-right: 0;' title='Sequenced Clone'></span>" : "";
		my $mtpCheck = ($fpcCloneMTP->{$_}) ? "<span class='ui-icon ui-icon-check' style='float: left; margin-right: 0;' title='MTP Clone'></span>" : "<span class='ui-icon ui-icon-bullet' style='float: left; margin-right: 0;'></span>";
		my $highlightedClass = ($fpcCloneHighlighted->{$_}) ? " class='ui-state-highlight'" : "";
		$sequenced++ if ($fpcCloneSequenced->{$_});
		if($colCount % $col == 0)
		{
			$fpcCtgDetails .= "<tr><td><div$highlightedClass style='float: left; margin-right: .7em;'>$mtpCheck$sequencedFlag<a onclick='closeDialog();openDialog(\"cloneView.cgi?cloneName=$_\")' title='$fpcClone->{$_}'>$_</a></div></td>";
		}
		elsif($colCount % $col == $col - 1)
		{
			$fpcCtgDetails .= "<td><div$highlightedClass style='float: left; margin-right: .7em;'>$mtpCheck$sequencedFlag<a onclick='closeDialog();openDialog(\"cloneView.cgi?cloneName=$_\")' title='$fpcClone->{$_}'>$_</a></div></td></tr>";
		}
		else
		{
			$fpcCtgDetails .= "<td><div$highlightedClass style='float: left; margin-right: .7em;'>$mtpCheck$sequencedFlag<a onclick='closeDialog();openDialog(\"cloneView.cgi?cloneName=$_\")' title='$fpcClone->{$_}'>$_</a></div></td>";
		}
		$colCount++;
    }
	if($colCount % $col > 0)
	{
		my $colspan = $col - $colCount % $col;
		$fpcCtgDetails .= "<td colspan='$colspan'></td></tr></table>";
	}
	else
	{
		$fpcCtgDetails .= "</table>";
	}
	$fpcCtgDetails ="<div id='fpcCtgList$fpcCtgId$$'>
	<ul class='fpcCtgListMenu' style='width: 150px;'>
		<li><a>About $fpcCtgList[2]</a>
			<ul>
				<li>Total Clones: $colCount ($sequenced sequenced)</li>
				<li>$fpcCtgList[8]</li>
			</ul>
		</li>
	</ul>
	<table width='100%'><tr>
	<td><b>Note:</b></td>
	<td><div style='float: left; margin-right: .7em;'><span class='ui-icon ui-icon-flag' style='float: left; margin-right: 0;' title='Sequenced Clone'></span>Sequenced</div></td>
	<td><div style='float: left; margin-right: .7em;'><span class='ui-icon ui-icon-check' style='float: left; margin-right: 0;' title='MTP Clone'></span>MTP</div></td>
	<td><div class='ui-state-highlight' style='float: left; margin-right: .7em;'>Highlighted</div></td>
	</tr></table>
	$fpcCtgDetails
	</div>";

	undef $/;# enable slurp mode
	my $html = <DATA>;
	$html =~ s/\$fpcCtgDetails/$fpcCtgDetails/g;
	$html =~ s/\$fpcCtgId/$fpcCtgId/g;
	$html =~ s/\$fpcCtg/$fpcCtgList[2]/g;
	$html =~ s/\$\$/$$/g;

	print header(-cookie=>cookie(-name=>'fpcCtg',-value=>$fpcCtgId));
	print $html;
}
else
{
	print header(-type=>'text/html',-status=>'402 Invalid operation');
	exit;
}


__DATA__
$fpcCtgDetails
<script>
buttonInit();
$( ".fpcCtgListMenu" ).menu();
$('#dialog').dialog("option", "title", "FPC $fpcCtg Clone List");
$( "#dialog" ).dialog( "option", "buttons", [ { text: "Print", click: function() {printDiv('fpcCtgList$fpcCtgId$$'); } }, { text: "OK", click: function() {closeDialog(); } } ] );
</script>