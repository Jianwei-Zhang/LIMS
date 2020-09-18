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

my @itemId = param ('itemId');
my $itemList = '';
my $alert = '';
my $col = 5;
my $colCount=0;
for(sort {$a <=> $b} @itemId)
{
	my $item = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$item->execute($_);
	my @item=$item->fetchrow_array();
	my $humanReadable = 'NA';
	if($item[1] eq "plate")
	{
		my $plateToLibrary=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$plateToLibrary->execute($item[6]);
		my @plateToLibrary = $plateToLibrary->fetchrow_array();
		$humanReadable = "$plateToLibrary[2]$item[2].$item[4]<sup class='ui-state-error-text'>*</sup>";
		$alert = '<h5><span class="ui-icon ui-icon-alert" style="float: left; margin-right: .3em;"></span>You will NOT be able to delete any items marked with an asterisk (<a class="ui-state-error-text">*</a>) from this system after printing their labels.</h5>';
	}
	else
	{
		$humanReadable = "$item[1]-$item[2]";
	}
	$itemList = "<table width='100%'>" unless ($itemList);
	if($colCount % $col == 0)
	{
		$itemList .= "<tr><td><input type='checkbox' id='itemList$item[0]$$' name='items' value='$item[0]' checked='checked'><label for='itemList$item[0]$$'>$humanReadable</label></td>";
	}
	elsif($colCount % $col == $col - 1)
	{
		$itemList .= "<td><input type='checkbox' id='itemList$item[0]$$' name='items' value='$item[0]' checked='checked'><label for='itemList$item[0]$$'>$humanReadable</label></td></tr>";
	}
	else
	{
		$itemList .= "<td><input type='checkbox' id='itemList$item[0]$$' name='items' value='$item[0]' checked='checked'><label for='itemList$item[0]$$'>$humanReadable</label></td>";
	}
	$colCount++;
}
if($itemList)
{
	my $toBeFilled = $col - ( $colCount % $col);
	$itemList .= ($toBeFilled < $col ) ? "<td>&nbsp;</td>" x $toBeFilled ."</tr></tbody></table>" : "</tbody></table>";
}
else
{
	$itemList = '<h5 class="ui-state-error-text"><span class="ui-icon ui-icon-alert" style="float: left; margin-right: .3em;"></span>Please check items first!</h5>';
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$alert/$alert/g;
$html =~ s/\$itemList/$itemList/g;

print header;
print $html;

__DATA__
<form id="printItem" name="printItem" action="barcodePrint.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<table>
	<tr><td style='text-align:left'><label for="printItemBlankLabel"><b>Leave</b></label> <input name="blankLabel" id="printItemBlankLabel" size="2" type="text" maxlength="2" VALUE="0" /> <label for="printItemBlankLabel"><b>Blank Label(s) at the top of the Sheet (US Letter)</b></label></td></tr>
	<tr><td style='text-align:left'><label for="printItemList"><b>Then Print Labels for the Selected Items:</b></label></td></tr>
	<tr><td>$itemList</td></tr>
	</table>
	$alert
</form>
<script>
$( "#printItemBlankLabel" ).spinner({ min: 0, max: 84});
$('#dialog').dialog("option", "title", "Print Labels");
$('#dialog').dialog("option", "width", "700");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Print", click: function() { submitForm('printItem'); closeDialog();} }, { text: "Cancel", click: function() { closeDialog(); } } ] );
</script>