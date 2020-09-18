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

my $assemblyId = param ('assemblyId') || '';

my $assembly=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$assembly->execute($assemblyId);
my @assembly = $assembly->fetchrow_array();

unless ($assembly[7] == 1 || $assembly[7] == 0) # exit if for frozen or running assembly
{
	print header;
	print <<END;
<script>
	closeDialog();
	errorPop("This assembly is running or frozen.");
</script>	
END
	exit;
}

my @ctgId;
my %ctg;
my %chrName;
my $col = 8;
my $colNumber=0;
my $assemblyCtg= $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblyCtg' AND o = ? ORDER BY name");
$assemblyCtg->execute($assemblyId);
while(my @assemblyCtg = $assemblyCtg->fetchrow_array())
{
	push @ctgId,$assemblyCtg[0];
	$ctg{$assemblyCtg[0]}=$assemblyCtg[2];
	$chrName{$assemblyCtg[0]}=$assemblyCtg[4];

}
@ctgId=sort {$ctg{$a} <=> $ctg{$b}} @ctgId;
my $ctgList = "<table id='ctgList$$' class='display'><thead style='display:none;'><tr>" . "<th></th>" x $col . "</tr></thead><tbody>";
for my $ctgId (@ctgId)
{
	if($colNumber % $col == 0)
	{
		$ctgList .= ( $chrName{$ctgId} > 0 ) ? "<tr><td><input type='checkbox' id='ctgList$ctg{$ctgId}$$' name='assemblyCtgId' value='$ctgId'><label for='ctgList$ctg{$ctgId}$$' title='$chrName{$ctgId}'><u>Ctg$ctg{$ctgId}</u></label></td>" : "<tr><td><input type='checkbox' id='ctgList$ctg{$ctgId}$$' name='assemblyCtgId' value='$ctgId'><label for='ctgList$ctg{$ctgId}$$' title='$chrName{$ctgId}'>Ctg$ctg{$ctgId}</label></td>";
	}
	elsif($colNumber % $col == $col - 1)
	{
		$ctgList .= ( $chrName{$ctgId} > 0 ) ? "<td><input type='checkbox' id='ctgList$ctg{$ctgId}$$' name='assemblyCtgId' value='$ctgId'><label for='ctgList$ctg{$ctgId}$$' title='$chrName{$ctgId}'><u>Ctg$ctg{$ctgId}</u></label></td></tr>" : "<td><input type='checkbox' id='ctgList$ctg{$ctgId}$$' name='assemblyCtgId' value='$ctgId'><label for='ctgList$ctg{$ctgId}$$' title='$chrName{$ctgId}'>Ctg$ctg{$ctgId}</label></td></tr>";
	}
	else
	{
		$ctgList .= ( $chrName{$ctgId} > 0 ) ? "<td><input type='checkbox' id='ctgList$ctg{$ctgId}$$' name='assemblyCtgId' value='$ctgId'><label for='ctgList$ctg{$ctgId}$$' title='$chrName{$ctgId}'><u>Ctg$ctg{$ctgId}</u></label></td>":"<td><input type='checkbox' id='ctgList$ctg{$ctgId}$$' name='assemblyCtgId' value='$ctgId'><label for='ctgList$ctg{$ctgId}$$' title='$chrName{$ctgId}'>Ctg$ctg{$ctgId}</label></td>";
	}
	$colNumber++;
}
my $toBeFilled = $col - ( $colNumber % $col);
$ctgList .= ($toBeFilled < $col ) ? "<td>&nbsp;</td>" x $toBeFilled ."</tr></tbody></table>" : "</tbody></table>";

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$ctgList/$ctgList/g;
$html =~ s/\$assemblyId/$assemblyId/g;
$html =~ s/\$commoncfg->{HTDOCS}/$commoncfg->{HTDOCS}/g;
print header;
print $html;

__DATA__
<button style='float: right; margin-top: .3em; margin-right: .3em;'  onClick='uncheckAll("assemblyCtgId");return false;'>Uncheck All</button>
<button style='float: right; margin-top: .3em; margin-right: .3em;'  onClick='checkAll("assemblyCtgId");return false;'>Check All</button>
<form id="assemblyAssign" name="assemblyAssign" action="assemblyAssignSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
<h3>Please check ctgs for assigning chromosome number:</h3>
$ctgList
<table>
<tr><td style='text-align:right'><label for="chr$$"><b>Chromosome #</b></label></td><td><input class='ui-widget-content ui-corner-all' name="chrNumber" id="editchrNumber" size="4" type="text" maxlength="3" ></td></tr>
</table>
</form>
<script>
buttonInit();
$('#dialog').dialog("option", "title", "Assign Ctg To Chromosome");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Assign", click: function() { submitForm('assemblyAssign'); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
</script>
