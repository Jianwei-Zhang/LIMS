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

my $boxId = param ('boxId') || '';
my $button = "<div>";
my $items = '';
my $engaged = 0;
my $maxLoads = 0;
if ($boxId)
{
	my $itemInBox=$dbh->prepare("SELECT matrix.* FROM link,matrix WHERE link.type LIKE 'box' AND link.child = matrix.id AND link.parent = $boxId ORDER BY matrix.o");
	$itemInBox->execute();
	while (my @itemInBox = $itemInBox->fetchrow_array())
	{
		$itemInBox[8] = escapeHTML($itemInBox[8]);
		$itemInBox[8] =~ s/\n/<br>/g if($itemInBox[1] ne 'plate');;
		$itemInBox[5] = ($itemInBox[4] == $itemInBox[5]) ? "" : "from ".$itemInBox[5];
		$items = "<form id='itemList$boxId$$' name='itemList$boxId$$'>
			<table id='itemsInBox$boxId$$' class='display' cellspacing='0' width='100%'>
				<thead>
					<tr>
						<th>
							<input type='checkbox' id='checkAllBox$$' name='checkAllBox$$' value='Check all' checked='checked' onClick='checkAll(\"itemId\");return false;' title='Check all'>
							<input type='checkbox' id='uncheckAllBox$$' name='uncheckAllBox$$' value='Uncheck all' onClick='uncheckAll(\"itemId\");return false;' title='Uncheck all'>
						</th>
						<th><b>Item</b></th>
						<th><b>Copy Number</b></th>
						<th><b>Creator</b></th>
						<th><b>Creation Date</b></th>
						<th><b>Barcode</b></th>
					</tr>
				</thead>
				<tbody>" unless($items);
		if ($itemInBox[1] eq 'plate')
		{
			my $plateToLibrary=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
			$plateToLibrary->execute($itemInBox[6]);
			my @plateToLibrary = $plateToLibrary->fetchrow_array();
			my $libraryToProject=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
			$libraryToProject->execute($plateToLibrary[6]);
			my @libraryToProject = $libraryToProject->fetchrow_array();
			$items .= "<tr>
				<td style='text-align:center;'><input type='checkbox' id='itemList$itemInBox[0]$$' name='itemId' value='$itemInBox[0]'></td>
				<td title='$itemInBox[8]'><div style='position: relative;'><a id='itemId$itemInBox[0]$$' onmouseover='editIconShow(\"itemId$itemInBox[0]$$\")' onmouseout='editIconHide(\"itemId$itemInBox[0]$$\")' onclick='openDialog(\"plateEdit.cgi?plateId=$itemInBox[0]\")' title='Edit/Delete'>$libraryToProject[2] > $plateToLibrary[2] > Plate $itemInBox[2]</a></div></td>
				<td title='copy number'>Copy $itemInBox[4] $itemInBox[5]</td>
				<td title='Creator'>$itemInBox[9]</td>
				<td title='Creation date'>$itemInBox[10]</td>
				<td style='text-align:center'><img alt='$itemInBox[7]' height='24' src='barcode.cgi?code=$itemInBox[7]'/></td>
				</tr>\n";
		}
		else
		{
			$items .= "<tr>
				<td style='text-align:center;' class='ui-widget-content'><input type='checkbox' id='itemList$itemInBox[0]$$' name='itemId' value='$itemInBox[0]'></td>
				<td class='ui-widget-content' title='$itemInBox[1]'><div style='position: relative;'><a id='itemId$itemInBox[0]$$' onmouseover='editIconShow(\"itemId$itemInBox[0]$$\")' onmouseout='editIconHide(\"itemId$itemInBox[0]$$\")' onclick='openDialog(\"itemEdit.cgi?itemId=$itemInBox[0]\")' title='Edit/Delete'>$itemInBox[8]</a></div></td>
				<td class='ui-widget-content' title='NA'></td>
				<td class='ui-widget-content' title='Creator'>$itemInBox[9]</td>
				<td class='ui-widget-content' title='Creation date'>$itemInBox[10]</td>
				<td class='ui-widget-content' style='text-align:center'><img alt='$itemInBox[7]' height='24' src='barcode.cgi?code=$itemInBox[7]'/></td>
				</tr>\n";
		}
	}
	$engaged = $itemInBox->rows;
	my $box=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$box->execute($boxId);
	my @box = $box->fetchrow_array();
	$box[8] = escapeHTML($box[8]);
	$box[8] =~ s/\n/<br>/g;
	$maxLoads= $box[4]*$box[5]*$box[6];
	my $disengaged = $maxLoads - $engaged;
	$button = "<div>";
	$button .= "<div style='float: right; margin-top: .3em; margin-right: .3em;' id='progressbar$boxId$$'><div class='progress-label'>$engaged/$maxLoads loads</div></div>";	
	$button .= "<button style='z-index: 1;float: right; margin-top: .3em; margin-right: .3em;' onclick='openDialogForm(\"itemDeleteForm.cgi\",\"itemList$boxId$$\")'>Delete item</button>
		<button style='z-index: 1;float: right; margin-top: .3em; margin-right: .3em;' onclick='openDialogForm(\"itemDepotForm.cgi\",\"itemList$boxId$$\")'>Move item</button>" if($items);
	$button .= "<div style='position: relative;'><h4><a id='box$boxId$$' onmouseover='editIconShow(\"box$boxId$$\")' onmouseout='editIconHide(\"box$boxId$$\")' onclick='openDialog(\"boxEdit.cgi?boxId=$boxId\")' title='Edit/Delete'>Box $box[2]</a><img alt='$box[7]' src='barcode.cgi?code=$box[7]'/><sup class='ui-state-disabled'>by $box[9] on $box[10]</sup></h4></div>";
	$button .= "$box[8]" if ($box[8]);	
	unless($items)
	{
		$button .= "<p class='ui-state-error ui-corner-all' style='padding: .7em;'><span class='ui-icon ui-icon-alert' style='float: left; margin-right: .3em;'></span>
			<strong>No item in this box.</strong></p>";
	}
	$items .= "</tbody></table></form>\n" if($items);
	$button .= "</div>";

	undef $/;# enable slurp mode
	my $html = <DATA>;
	$html =~ s/\$boxId/$boxId/g;
	$html =~ s/\$\$/$$/g;
	$html =~ s/\$maxLoads/$maxLoads/g;
	$html =~ s/\$engaged/$engaged/g;
	$html =~ s/\$button/$button/g;
	$html =~ s/\$items/$items/g;
	print header(-cookie=>cookie(-name=>'box',-value=>$boxId));
	print $html;
}
else
{
	print header(-type=>'text/html',-status=>'402 Invalid operation');
	exit;
}

__DATA__
$button
$items
<script>
buttonInit();
$( "#progressbar$boxId$$" ).progressbar({
	max:$maxLoads,
	value: $engaged
});
loadingHide();
$( "#itemsInBox$boxId$$" ).dataTable({
	"scrollY": "450px",
	"scrollCollapse": true,
	"paging": false,
	"order": [ 4, 'desc' ],
	"columnDefs": [
    { "orderable": false, "targets": [0, -1] }
  ]
});
</script>