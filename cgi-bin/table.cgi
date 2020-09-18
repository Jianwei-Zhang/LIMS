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

my $type = param ('type') || '';
my $parentId = param ('parentId') || '';
my $refresh = param ('refresh') || 'menu';
my $dataset = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$dataset->execute($parentId);
my @dataset = $dataset->fetchrow_array();

my $table = '';
my $allTable=$dbh->prepare("SELECT * FROM matrix WHERE container Like ? AND z = ?");
$allTable->execute($type,$parentId);
while (my @allTable = $allTable->fetchrow_array())
{
 	my $itemDetails = decode_json $allTable[8];
	unless ($table)
	{
		$table .= "<form id='tableList$$' name='tableList$$'>
			<table id='table$$' class='display' style='width: 100%;'>
				<thead>
					<tr>
						<th>
						<input type='checkbox' id='checkAllBox$$' name='checkAllBox$$' value='Check all' checked='checked' onClick='checkAll(\"itemId\");return false;' title='Check all'>
						<input type='checkbox' id='uncheckAllBox$$' name='uncheckAllBox$$' value='Uncheck all' onClick='uncheckAll(\"itemId\");return false;' title='Uncheck all'>
						</th>";

		my $headColumn = 0;
		for (sort {$a <=> $b} keys %$itemDetails)
		{
			if ($headColumn > 5)
			{
				$table .= "<th style='text-align:left'><b>More</b></th>";
				last;
			}
			$itemDetails->{$_}->{'field'} = '' unless ($itemDetails->{$_}->{'field'});
			$itemDetails->{$_}->{'value'} = '' unless ($itemDetails->{$_}->{'value'});
			$itemDetails->{$_}->{'value'} = escapeHTML($itemDetails->{$_}->{'value'});
			$itemDetails->{$_}->{'value'} =~ s/\n/<br>/g;
			$table .= "<th style='text-align:left'><b>$itemDetails->{$_}->{'field'}</b></th>";
			$headColumn++;
		}

		$table .= "<th style='text-align:left'><b>Creation Date</b></th>
					</tr>
				</thead>
				<tbody>";
	}

	$table .= "<tr><td style='text-align:center;'><input type='checkbox' id='record$allTable[0]$$' name='itemId' value='$allTable[0]'></td>";

	my $column = 0;
	for (sort {$a <=> $b} keys %$itemDetails)
	{
		if ($column > 5)
		{
			$table .= "<td><a id='moreId$allTable[0]$$' onclick='openDialog(\"itemView.cgi?itemId=$allTable[0]\")' title='View'>...</a></td>";
			last;
		}
		$itemDetails->{$_}->{'field'} = '' unless ($itemDetails->{$_}->{'field'});
		$itemDetails->{$_}->{'value'} = '' unless ($itemDetails->{$_}->{'value'});
		$itemDetails->{$_}->{'value'} = escapeHTML($itemDetails->{$_}->{'value'});
		$itemDetails->{$_}->{'value'} =~ s/\n/<br>/g;
		if($itemDetails->{$_}->{'value'} =~ /\.(jpg|jpeg|png|tif|tiff)$/i)
		{
			$table .= "<td>";
			for (split "\;", $itemDetails->{$_}->{'value'})
			{
				$_ =~ s/^\s+|\s+$//g;
				$table .= "<img src='$commoncfg->{HTDOCS}/data/images/$_'/>";
			}
			$table .= "</td>";
		}
		else
		{
			$table .= ($itemDetails->{$_}->{'value'} =~ /:\/\//) ? "<td><a href='$itemDetails->{$_}->{'value'}' target='_blank'>$itemDetails->{$_}->{'value'}</a></td>" :
						($column > 1) ? "<td>$itemDetails->{$_}->{'value'}</td>" : "<td><a id='tableId$allTable[0]$$' onclick='openDialog(\"itemView.cgi?itemId=$allTable[0]\")' title='View'>$itemDetails->{$_}->{'value'}</a></td>";
		}
		$column++;
	}

	$table .= "<td title='Creator: $allTable[9]'>$allTable[10]</td>
		</tr>";
}
$table .= "</tbody></table></form>\n" if($table);
my $button = "<div class='ui-state-highlight ui-corner-all' style='padding: 0 .7em;'>
	<button style='float: right; margin-top: .3em; margin-right: .3em;' onclick='openDialogForm(\"itemDeleteForm.cgi\",\"tableList$$\")' title='Delete $type'>Delete $type</button>
	<button style='float: right; margin-top: .3em; margin-right: .3em;' onclick='openDialogForm(\"tableDownload.cgi\",\"tableList$$\")' title='Download $type'>Download $type</button>
	<button style='float: right; margin-top: .3em; margin-right: .3em;' onclick='openDialog(\"tableNew.cgi?type=$type&parentId=$parentId&refresh=$refresh\")'>New $type(not ready yet)</button>
	<h3>List of <a onclick='openDialog(\"datasetView.cgi?datasetId=$dataset[0]\")' title='View'>$dataset[2]</a></h3>
	";

unless($table)
{
	$button .= "<p class='ui-state-error ui-corner-all' style='padding: .7em;'><span class='ui-icon ui-icon-alert' style='float: left; margin-right: .3em;'></span>
			<strong>No data available!</strong></p>";
}
$button .= "</div>\n";

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$button/$button/g;
$html =~ s/\$table/$table/g;
$html =~ s/\$\$/$$/g;

print header();
print $html;

__DATA__
$button
$table
<script>
buttonInit();
$( "#table$$" ).dataTable({
	"scrollY": "600px",
	"scrollCollapse": true,
	"paging": false
});
loadingHide();
</script>