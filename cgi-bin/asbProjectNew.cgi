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

my $col = 3;
my $colCount=0;
my $assemblyTargetIds = "<table id='assemblyTargetIds$$' class='display'><thead style='display:none;'><tr>" . "<th></th>" x $col . "</tr></thead><tbody>";
my $library = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'library'");# ORDER BY name
$library->execute();
if($library->rows > 0)
{
	my $libraryResult;
	while (my @library=$library->fetchrow_array())
	{
		@{$libraryResult->{$library[2]}} = @library;
	}
	foreach (sort {uc ($a) cmp uc($b)} keys %$libraryResult)
	{
		my @library = @{$libraryResult->{$_}};
		if($colCount % $col == 0)
		{
			$assemblyTargetIds .= "<tr><td><input type='checkbox' id='libraryList$library[0]$$' name='targetId' value='$library[0]'><label for='libraryList$library[0]$$' title='library'>$library[2]<sup class='ui-state-disabled'>L</sup></label></td>";
		}
		elsif($colCount % $col == $col - 1)
		{
			$assemblyTargetIds .= "<td><input type='checkbox' id='libraryList$library[0]$$' name='targetId' value='$library[0]'><label for='libraryList$library[0]$$' title='library'>$library[2]<sup class='ui-state-disabled'>L</sup></label></td></tr>";
		}
		else
		{
			$assemblyTargetIds .= "<td><input type='checkbox' id='libraryList$library[0]$$' name='targetId' value='$library[0]'><label for='libraryList$library[0]$$' title='library'>$library[2]<sup class='ui-state-disabled'>L</sup></label></td>";
		}
		$colCount++;
	}
}

my $genome = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'genome'");# ORDER BY name
$genome->execute();
if($genome->rows > 0)
{
	my $genomeResult;
	while (my @genome=$genome->fetchrow_array())
	{
		next if ($genome[4] < 1);
		@{$genomeResult->{$genome[2]}} = @genome;
	}
	foreach (sort {uc ($a) cmp uc($b)} keys %$genomeResult)
	{
		my @genome = @{$genomeResult->{$_}};
		if($colCount % $col == 0)
		{
			$assemblyTargetIds .= "<tr><td><input type='checkbox' id='genomeList$genome[0]$$' name='targetId' value='$genome[0]'><label for='genomeList$genome[0]$$' title='genome'>$genome[2]<sup>G</sup></label></td>";
		}
		elsif($colCount % $col == $col - 1)
		{
			$assemblyTargetIds .= "<td><input type='checkbox' id='genomeList$genome[0]$$' name='targetId' value='$genome[0]'><label for='genomeList$genome[0]$$' title='genome'>$genome[2]<sup class='ui-state-disabled'>G</sup></label></td></tr>";
		}
		else
		{
			$assemblyTargetIds .= "<td><input type='checkbox' id='genomeList$genome[0]$$' name='targetId' value='$genome[0]'><label for='genomeList$genome[0]$$' title='genome'>$genome[2]<sup class='ui-state-disabled'>G</sup></label></td>";
		}
		$colCount++;
	}
}

my $toBeFilled = $col - ( $colCount % $col);
$assemblyTargetIds .= ($toBeFilled < $col ) ? "<td>&nbsp;</td>" x $toBeFilled ."</tr></tbody></table>" : "</tbody></table>";

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$assemblyTargetIds/$assemblyTargetIds/g;
$html =~ s/\$\$/$$/g;

print header;
print $html;

__DATA__
<div id="asbProjectNew" class="ui-widget-content ui-corner-all" style='padding: 0 .7em;'>
	<h3>New GPM Project</h3>
	<form id="newAsbProject" name="newAsbProject" action="asbProjectSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<table>
	<tr><td style='text-align:right'><label for="newAsbProjectName"><b>Project Name</b></label></td><td><input class='ui-widget-content ui-corner-all' name="name" id="newAsbProjectName" size="40" type="text" maxlength="32" /></td></tr>
	<tr><td style='text-align:right'><label for="newAsbProjectDescription"><b>Description</b></label></td><td><textarea class='ui-widget-content ui-corner-all' name="description" id="newAsbProjectDescription" cols="60" rows="10"></textarea></td></tr>
	<tr><td style='text-align:right'><label for="newAsbProjectTarget"><b>Target Libraries/Genomes</b></label></td><td>
	$assemblyTargetIds
	</td></tr>
	<tr><td></td><td><INPUT TYPE="button" VALUE="Save" onclick="submitForm('newAsbProject');"></td></tr>
	</table>
	</form>
</div>
<script>
buttonInit();
loadingHide();
$( "#assemblyTargetIds$$" ).dataTable({
	"dom": 'lfrtip',
	"scrollY": "300px",
	"scrollCollapse": true,
	"paging": false,
	"searching": false,
	"ordering": false,
	"info": false
});
</script>