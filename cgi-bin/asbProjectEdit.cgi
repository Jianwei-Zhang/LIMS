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

my $asbProjectId = param ('asbProjectId') || '';

my $asbProject = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$asbProject->execute($asbProjectId);
my @asbProject=$asbProject->fetchrow_array();

my $checkAsbProject = $dbh->prepare("SELECT child FROM link WHERE parent = ? AND type LIKE 'asbProject'");
$checkAsbProject->execute($asbProjectId);
my $checkedTargetId;
while(my @checkAsbProject=$checkAsbProject->fetchrow_array())
{
	my $assembly = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assembly' AND x = ?");
	$assembly->execute($checkAsbProject[0]);
	$checkedTargetId->{$checkAsbProject[0]} = $assembly->rows;
	if($checkedTargetId->{$checkAsbProject[0]} > 0)
	{
		#check if the same library or genome belongs to different asbProject
		my $checkAsbProjectByChild = $dbh->prepare("SELECT parent FROM link WHERE child = ? AND type LIKE 'asbProject'");
		$checkAsbProjectByChild->execute($checkAsbProject[0]);
		if($checkAsbProjectByChild->rows > 1)
		{
			$checkedTargetId->{$checkAsbProject[0]} = 0; #make target deletable
		}
	}
}

my $col = 2;
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
		my $checked = (exists $checkedTargetId->{$library[0]}) ? ($checkedTargetId->{$library[0]}) ? "onchange='this.checked=true' checked='checked'" : "checked='checked'" : "";
		if($colCount % $col == 0)
		{
			$assemblyTargetIds .= "<tr><td><input type='checkbox' id='libraryList$library[0]$$' name='targetId' value='$library[0]' $checked><label for='libraryList$library[0]$$' title='library'>$library[2]<sup class='ui-state-disabled'>L</sup></label></td>";
		}
		elsif($colCount % $col == $col - 1)
		{
			$assemblyTargetIds .= "<td><input type='checkbox' id='libraryList$library[0]$$' name='targetId' value='$library[0]' $checked><label for='libraryList$library[0]$$' title='library'>$library[2]<sup class='ui-state-disabled'>L</sup></label></td></tr>";
		}
		else
		{
			$assemblyTargetIds .= "<td><input type='checkbox' id='libraryList$library[0]$$' name='targetId' value='$library[0]' $checked><label for='libraryList$library[0]$$' title='library'>$library[2]<sup class='ui-state-disabled'>L</sup></label></td>";
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
		my $checked = (exists $checkedTargetId->{$genome[0]}) ? ($checkedTargetId->{$genome[0]}) ? "onchange='this.checked=true' checked='checked'" : "checked='checked'" : "";
		if($colCount % $col == 0)
		{
			$assemblyTargetIds .= "<tr><td><input type='checkbox' id='genomeList$genome[0]$$' name='targetId' value='$genome[0]' $checked><label for='genomeList$genome[0]$$' title='genome'>$genome[2]<sup class='ui-state-disabled'>G</sup></label></td>";
		}
		elsif($colCount % $col ==  $col - 1)
		{
			$assemblyTargetIds .= "<td><input type='checkbox' id='genomeList$genome[0]$$' name='targetId' value='$genome[0]' $checked><label for='genomeList$genome[0]$$' title='genome'>$genome[2]<sup class='ui-state-disabled'>G</sup></label></td></tr>";
		}
		else
		{
			$assemblyTargetIds .= "<td><input type='checkbox' id='genomeList$genome[0]$$' name='targetId' value='$genome[0]' $checked><label for='genomeList$genome[0]$$' title='genome'>$genome[2]<sup class='ui-state-disabled'>G</sup></label></td>";
		}
		$colCount++;
	}
}

my $toBeFilled = $col - ( $colCount % $col);
$assemblyTargetIds .= ($toBeFilled < $col ) ? "<td>&nbsp;</td>" x $toBeFilled ."</tr></tbody></table>" : "</tbody></table>";

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$asbProjectId/$asbProjectId/g;
$html =~ s/\$asbProjectName/$asbProject[2]/g;
$html =~ s/\$asbProjectDescription/$asbProject[8]/g;
$html =~ s/\$assemblyTargetIds/$assemblyTargetIds/g;
$html =~ s/\$\$/$$/g;
print header;
print $html;

__DATA__
<form id="editAsbProject" name="editAsbProject" action="asbProjectSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="asbProjectId" id="editAsbProjectId" type="hidden" value="$asbProjectId" />
	<table>
	<tr><td style='text-align:right'><label for="editAsbProjectName"><b>Project Name</b></label></td><td><input class='ui-widget-content ui-corner-all' name="name" id="editAsbProjectName" size="40" type="text" maxlength="32" value="$asbProjectName"/></td></tr>
	<tr><td style='text-align:right'><label for="editAsbProjectDescription"><b>Description</b></label></td><td><textarea class='ui-widget-content ui-corner-all' name="description" id="editAsbProjectDescription" cols="50" rows="10">$asbProjectDescription</textarea></td></tr>
	<tr><td style='text-align:right'><label for="editAsbProjectTarget"><b>Target Libraries/Genomes</b></label></td><td>
	$assemblyTargetIds
	</td></tr>
	</table>
</form>
<script>
$( "#assemblyTargetIds$$" ).dataTable({
	"dom": 'lfrtip',
	"scrollY": "300px",
	"scrollCollapse": true,
	"paging": false,
	"searching": false,
	"ordering": false,
	"info": false
});
$('#dialog').dialog("option", "title", "Edit GPM Project");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('editAsbProject'); } }, { text: "Delete", click: function() { deleteItem($asbProjectId); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
</script>