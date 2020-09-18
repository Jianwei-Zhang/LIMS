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

my $vectors = '';
my $allVector=$dbh->prepare("SELECT * FROM matrix WHERE container Like 'vector'");
$allVector->execute();
while (my @allVector = $allVector->fetchrow_array())
{
	$allVector[8] = escapeHTML($allVector[8]);
	$allVector[8] =~ s/\n/<br>/g;
	$allVector[2] = "Name N/A" unless($allVector[2]);
	my $relatedLibraries = '';
	my $relatedLibrary=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'library' AND x = $allVector[0]");
	$relatedLibrary->execute();
	while (my @relatedLibrary = $relatedLibrary->fetchrow_array())
	{
		$relatedLibraries .= "<a onclick='openDialog(\"libraryView.cgi?libraryId=$relatedLibrary[0]\")'>$relatedLibrary[2]</a> ";
	}
	$vectors = "<form id='vectorList$$' name='vectorList$$'>
		<table id='vectors$$' class='display' style='width: 100%;'>
			<thead>
				<tr>
					<th>
						<input type='checkbox' id='checkAllBox$$' name='checkAllBox$$' value='Check all' checked='checked' onClick='checkAll(\"itemId\");return false;' title='Check all'>
						<input type='checkbox' id='uncheckAllBox$$' name='uncheckAllBox$$' value='Uncheck all' onClick='uncheckAll(\"itemId\");return false;' title='Uncheck all'>
					</th>
					<th style='text-align:left'><b>Vector</b></th>
					<th style='text-align:left'><b>Related libraries</b></th>
					<th style='text-align:left'><b>Creator</b></th>
					<th style='text-align:left'><b>Creation Date</b></th>
				</tr>
			</thead>
			<tbody>" unless($vectors);
	$vectors .= "<tr>
		<td style='text-align:center;'><input type='checkbox' id='vectorList$allVector[0]$$' name='itemId' value='$allVector[0]'></td>
		<td title='Vector'><div style='position: relative;'><a id='vectorId$allVector[0]$$' onmouseover='editIconShow(\"vectorId$allVector[0]$$\")' onmouseout='editIconHide(\"vectorId$allVector[0]$$\")' onclick='openDialog(\"vectorEdit.cgi?vectorId=$allVector[0]\")' title='Edit/Delete'>$allVector[2]</a></div></td>
		<td title='Related libraries'>$relatedLibraries</td>
		<td title='Creator'>$allVector[9]</td>
		<td title='Creation date'>$allVector[10]</td>
		</tr>";
}
$vectors .= "</tbody></table></form>\n" if($vectors);
my $button = "<div class='ui-state-highlight ui-corner-all' style='padding: 0 .7em;'>
	<button style='float: right; margin-top: .3em; margin-right: .3em;' onclick='openDialog(\"vectorNew.cgi\")'>New vector</button>";
$button .= "<button style='float: right; margin-top: .3em; margin-right: .3em;' onclick='openDialogForm(\"itemDeleteForm.cgi\",\"vectorList$$\")'>Delete</button>" if($vectors);
$button .= "<button style='float: right; margin-top: .3em; margin-right: .3em;' onclick='openDialogForm(\"vectorMergeForm.cgi\",\"vectorList$$\")'>Merge</button>" if($allVector->rows > 1);
$button .= "<button style='float: right; margin-top: .3em; margin-right: .3em;' onclick='refresh(\"general\")'>Refresh</button>
			<h2>Vectors</h2>";

unless($vectors)
{
	$button .= "<p class='ui-state-error ui-corner-all' style='padding: .7em;'><span class='ui-icon ui-icon-alert' style='float: left; margin-right: .3em;'></span>
			<strong>No vector, please create one!</strong></p>";
}
$button .= "</div>\n";

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$button/$button/g;
$html =~ s/\$vectors/$vectors/g;
$html =~ s/\$\$/$$/g;


print header(-cookie=>cookie(-name=>'general',-value=>0));
print $html;

__DATA__
$button
$vectors
<script>
buttonInit();
$( "#vectors$$" ).dataTable({
	"scrollY": "600px",
	"scrollCollapse": true,
	"paging": false,
	"order": [ 1, 'asc' ],
	"columnDefs": [
    { "orderable": false, "targets": 0 }
  ]
});
loadingHide();
</script>