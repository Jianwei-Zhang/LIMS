#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use JSON::XS;
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

my $libraryId = param ('libraryId') || '';
my $tableHeader;
my $plates = '';
my $nonOrderableTargets = '';
if ($libraryId)
{
	my $library=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$library->execute($libraryId);
	my @library = $library->fetchrow_array();	
	my $plateId;
	my $copyNumber;
	my $allCopyNumber;
	my $barcode;
	my $note;
	my $plateInLibrary=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'plate' AND z = ? ORDER BY o");
	$plateInLibrary->execute($libraryId);
	while (my @plateInLibrary = $plateInLibrary->fetchrow_array())
	{
		$plateInLibrary[2] = "Name N/A" unless($plateInLibrary[2]);
		$plateId->{$plateInLibrary[2]}->{$plateInLibrary[4]} = $plateInLibrary[0];
		$copyNumber->{$plateInLibrary[2]}->{$plateInLibrary[4]} = $plateInLibrary[5];
		$allCopyNumber->{$plateInLibrary[4]} = 1;
		$barcode->{$plateInLibrary[2]}->{$plateInLibrary[4]} = $plateInLibrary[7];
		$note->{$plateInLibrary[2]}->{$plateInLibrary[4]} = $plateInLibrary[8];
	}
	my $engaged = $plateInLibrary->rows;
	my @allCopyNumber = sort {$a <=> $b} keys %{$allCopyNumber};
	my $number = 1;
	for(@allCopyNumber)
	{
		$nonOrderableTargets .= ( $nonOrderableTargets ) ? ",$number" : $number; 
		$number++;
	}
	my $lastCopy = pop @allCopyNumber;
	$tableHeader = "<th style='text-align:left;'>
		<b>Copy $lastCopy</b> <input type='checkbox' id='checkAllPlate$lastCopy$$' name='checkAllPlate$lastCopy$$' value='Check all' checked='checked' onClick='checkClass(\"copy$lastCopy\");return false;' title='Check all Copy $lastCopy'>
		<input type='checkbox' id='uncheckAllPlate$lastCopy$$' name='uncheckAllPlate$lastCopy$$' value='Uncheck all' onClick='uncheckClass(\"copy$lastCopy\");return false;' title='Uncheck all Copy $lastCopy'>
		</th>";
	for (sort {$b <=> $a} @allCopyNumber)
	{
		$tableHeader = "<th style='text-align:left;'>
		<b>Copy $_</b> <input type='checkbox' id='checkAllPlate$_$$' name='checkAllPlate$_$$' value='Check all' checked='checked' onClick='checkClass(\"copy$_\");return false;' title='Check all Copy $_'>
		<input type='checkbox' id='uncheckAllPlate$_$$' name='uncheckAllPlate$_$$' value='Uncheck all' onClick='uncheckClass(\"copy$_\");return false;' title='Uncheck all Copy $_'>
		</th>" . $tableHeader;
	}
	for (sort keys %{$plateId})
	{
		my $plate=$_;
		$plates = "<ul id='plateInfoMenu$libraryId$$' style='margin-top: .3em;width: 250px;'>
				<li><a title='$library[2]'><b>Library '$library[2]' Plates</b></a>
					<ul style='z-index: 1000;white-space: nowrap;'>
						<li><a onclick='openDialog(\"plateNew.cgi?libraryId=$libraryId\")' title='New Plate'>New Plate</a></li>
						<li><a onclick='openDialogForm(\"itemDepotForm.cgi\",\"plateList$libraryId$$\")' title='Depot Plate'>Depot Plate</a></li>
						<li><a onclick='openDialogForm(\"itemDeleteForm.cgi\",\"plateList$libraryId$$\")' title='Delete Plate'>Delete Plate</a></li>
						<li><a onclick='openDialog(\"plateCopyForm.cgi?libraryId=$libraryId\")' title='Copy Plate'>Copy Plate</a></li>
						<li><a title='Print'>Print</a>
							<ul style='z-index: 1000;white-space: nowrap;'>
								<li><a onclick='openDialogForm(\"barcodePrintForm.cgi\",\"plateList$libraryId$$\")' title='Print Labels'>Labels</a></li>
								<li><a onclick='printDiv(\"plateList$libraryId$$\")' title='Print Plate List'>Plate List</a></li>
							</ul>
						</li>
					</ul>
				</li>
			</ul>
			<form id='plateList$libraryId$$' name='plateList$libraryId$$'>
			<table id='plateInLibrary$libraryId$$' name='plateInLibrary$libraryId$$' class='display'>
			<thead><tr><th><b>Plate #</b></th>$tableHeader</tr></thead><tbody>" unless ($plates);
		$plates .= "<tr><td style='text-align:center'>$plate</td>";	
		for(sort {$a <=> $b} keys %{$allCopyNumber})
		{
			if(exists $copyNumber->{$plate}->{$_})
			{
				$plates .= "<td style='text-align:left;'><div style='position: relative;float:left;'><a id='plateId$plateId->{$plate}->{$_}$$' onmouseover='editIconShow(\"plateId$plateId->{$plate}->{$_}$$\")' onmouseout='editIconHide(\"plateId$plateId->{$plate}->{$_}$$\")' onclick='openDialog(\"plateEdit.cgi?plateId=$plateId->{$plate}->{$_}\")' title='View/Edit/Delete'>$library[2]$plate.$_<br><img alt='$barcode->{$plate}->{$_}' src='barcode.cgi?code=$barcode->{$plate}->{$_}&notext=1&height=12'/></a><br><label for='plateList$libraryId$plateId->{$plate}->{$_}$$'>$barcode->{$plate}->{$_}</label><input type='checkbox' class='copy$_' id='plateList$libraryId$plateId->{$plate}->{$_}$$' name='itemId' value='$plateId->{$plate}->{$_}'></div></td>";
			}
			else
			{
				$plates .= "<td class='ui-state-error' style='text-align:left;'>NA</td>";
			}
		}	
		$plates .= "</tr>";	
	}
	if ($plates)
	{
		$plates .= "</tbody></table></form>" 
	}
	else
	{
		$plates = "<p class='ui-state-error ui-corner-all' style='padding: .7em;'><span class='ui-icon ui-icon-alert' style='float: left; margin-right: .3em;'></span>
			<button style='float: right; margin-top: -.4em; margin-right: .3em;' onclick='openDialog(\"plateNew.cgi?libraryId=$libraryId\")'>New plate</button>
			<strong>No plate in this library.</strong></p>";
	}
	
	undef $/;# enable slurp mode
	my $html = <DATA>;
	$html =~ s/\$libraryId/$libraryId/g;
	$html =~ s/\$\$/$$/g;
	$html =~ s/\$plates/$plates/g;
	$html =~ s/\$nonOrderableTargets/$nonOrderableTargets/g;

	print header;
	print $html;
}
else
{
	print header(-type=>'text/html',-status=>'402 Invalid operation');
	exit;
}

__DATA__
$plates
<script>
buttonInit();
$( "#plateInLibrary$libraryId$$" ).dataTable({
	"scrollCollapse": true,
	"paging": false,
	"columnDefs": [
    { "orderable": false, "targets": [$nonOrderableTargets] }
  ]
});
$( "#plateInfoMenu$libraryId$$" ).menu();
loadingHide();
</script>