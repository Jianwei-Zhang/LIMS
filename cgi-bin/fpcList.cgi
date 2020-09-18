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

my $fpcId = param ('fpcId') || '';
my $fpcDetails = '';
if ($fpcId)
{
	my $fpc=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$fpc->execute($fpcId);
	my @fpc = $fpc->fetchrow_array();
	$fpcDetails .= "<ul class='fpcMenu' style='left: 0px;top: 0px;display:inline-block;width: 120px;'>
				<li><a><span class='ui-icon ui-icon-wrench'></span>FPC Tools</a>
					<ul style='z-index: 1000;white-space: nowrap;'>
						<li><a href='download.cgi?fpcId=$fpc[0]' target='hiddenFrame'><span class='ui-icon ui-icon-disk'></span>Download Clone List</a></li>
						<li><a href='fpcPickMTP.cgi?fpcId=$fpc[0]' target='hiddenFrame'><span class='ui-icon ui-icon-search'></span>Pick Up MTP Clone</a></li>
						<li><a onclick='openDialog(\"fpcMarkCloneForm.cgi?fpcId=$fpc[0]\")'><span class='ui-icon ui-icon-flag'></span>Mark Clones</a></li>
						<li><a onclick='deleteItem($fpc[0])'><span class='ui-icon ui-icon-trash'></span>Delete This FPC</a></li>
					</ul>
				</li>
			</ul><hr>";

	my $fpcCtg;
	my $fpcCtgId;
	my $fpcCtgCloneNumber;
	my $fpcCtgChr;
	my $fpcCtgOrder;

	my $fpcCtgList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'fpcCtg' AND o = ?");
	$fpcCtgList->execute($fpcId);
	while (my @fpcCtgList = $fpcCtgList->fetchrow_array())
	{
		$fpcCtg->{$fpcCtgList[2]} = $fpcCtgList[8];
		$fpcCtgId->{$fpcCtgList[2]} = $fpcCtgList[0];
		$fpcCtgCloneNumber->{$fpcCtgList[2]} = $fpcCtgList[4];
		$fpcCtgChr->{$fpcCtgList[2]} = $fpcCtgList[5];
		$fpcCtgOrder->{$fpcCtgList[2]} = $fpcCtgList[6];
	}
	my @fpcCtg = map  { $_->[0] }
             sort { $a->[1] <=> $b->[1] }
             map  { [$_, $_=~/(\d+)/] }
             keys %$fpcCtg;

	$fpcDetails .= "<ul id='fpcCtgList$fpcId'>";
	my $col = 6;
	my $colCount=0;
    for (@fpcCtg)
    {
    	next if ($_ eq 'Ctg0');
		$fpcDetails .= "<li class='ui-state-default' id='$fpcCtgId->{$_}'><a onclick='closeViewer();openViewer(\"fpcCtgView.cgi?fpcCtgId=$fpcCtgId->{$_}\")'>$_</a><br>$fpcCtgCloneNumber->{$_} clones</li>";
		$colCount++;
    }
	$fpcDetails .= "</ul>";
	
	undef $/;# enable slurp mode
	my $html = <DATA>;
	$html =~ s/\$fpcDetails/$fpcDetails/g;
	$html =~ s/\$fpcId/$fpcId/g;

	print header(-cookie=>cookie(-name=>'fpc',-value=>$fpcId));
	print $html;
}
else
{
	print header(-type=>'text/html',-status=>'402 Invalid operation');
	exit;
}

__DATA__
$fpcDetails

<style>
#fpcCtgList$fpcId { list-style-type: none; display:inline-block;margin: 0; padding: 0; width: 100%; }
#fpcCtgList$fpcId li { margin: 3px 3px 3px 0; padding: 1px; float: left; width: 120px; text-align: center; }
</style>
<script>
buttonInit();
$( ".fpcMenu" ).menu();
loadingHide();
</script>