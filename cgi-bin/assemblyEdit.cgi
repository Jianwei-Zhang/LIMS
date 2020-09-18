#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use JSON::XS; #JSON::XS is recommended to be installed for handling JSON string of big size 
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

my $assemblyStatus;
$assemblyStatus->{2} = "Frozen";
$assemblyStatus->{1} = "Assembled";
$assemblyStatus->{0} = "Initialized";
$assemblyStatus->{-1} = "Assembling...";
$assemblyStatus->{-2} = "SeqToSeq Aligning...";
$assemblyStatus->{-3} = "Physical Reference-based Assembling...";
$assemblyStatus->{-4} = "SeqToGenome Aligning...";
$assemblyStatus->{-5} = "Assigning Chr#...";
$assemblyStatus->{-6} = "Filtering Seqs...";
$assemblyStatus->{-7} = "Orienting Seqs...";
$assemblyStatus->{-8} = "Orienting Contigs...";
$assemblyStatus->{-9} = "EndToEnd Merging...";
$assemblyStatus->{-10} = "Filtering Overlaps...";
$assemblyStatus->{-11} = "Estimating Length...";
$assemblyStatus->{-12} = "BesToSeq Aligning...";

my $objectComponent;
$objectComponent->{0} = "Unknown";
$objectComponent->{1} = "Chr-Seq";
$objectComponent->{2} = "Ctg-Seq";

my $assemblyId = param ('assemblyId') || '';
my $assembly=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$assembly->execute($assemblyId);
my @assembly = $assembly->fetchrow_array();
my $target=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$target->execute($assembly[4]);
my @target = $target->fetchrow_array();

my $fpcOrAgpId = '';
my $doNotChangeFpc = ($assembly[7] != 0 && $assembly[6] > 0) ? '<sup class="ui-state-error ui-corner-all">Do NOT Change Unless You Plan To Re-Run Assembly!</sup><br>' : '';
if($target[1] eq "library")
{
	my $fpcList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'fpc' AND barcode = ?");
	$fpcList->execute($target[0]);
	while (my @fpcList = $fpcList->fetchrow_array())
	{
		$fpcOrAgpId .= ($fpcList[0] eq $assembly[6] ) ?
			"<option value='$fpcList[0]' title='$fpcList[8]' selected>FPC: $fpcList[2] v.$fpcList[3]</option>" :
			"<option value='$fpcList[0]' title='$fpcList[8]'>FPC: $fpcList[2] v.$fpcList[3]</option>";
	}
}
elsif($target[1] eq "genome")
{
	my $agpList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'agp' AND x = ?");
	$agpList->execute($target[0]);
	while (my @agpList = $agpList->fetchrow_array())
	{
		$fpcOrAgpId .= ($agpList[0] eq $assembly[6] ) ?
			"<option value='$agpList[0]' selected>AGP: $agpList[2] v.$agpList[3] ($objectComponent->{$agpList[5]})</option>" :
			"<option value='$agpList[0]'>AGP: $agpList[2] v.$agpList[3] ($objectComponent->{$agpList[5]})</option>";
	}
}
if($fpcOrAgpId)
{
	$fpcOrAgpId = "<option class='ui-state-error-text' value='0'>Please select a reference</option>".$fpcOrAgpId;
}
else
{
	$fpcOrAgpId = "<option class='ui-state-error-text' value='0'>No reference available</option>".$fpcOrAgpId;
}


my $checkAsbGenome = $dbh->prepare("SELECT child FROM link WHERE parent = ? AND type LIKE 'asbGenome'");
$checkAsbGenome->execute($assemblyId);
my $checkedExtraId;
while(my @checkAsbGenome=$checkAsbGenome->fetchrow_array())
{
	$checkedExtraId->{$checkAsbGenome[0]} = 1;
}

my $col = 2;
my $colCount=0;
my $assemblyExtraIds = "<table id='assemblyExtraIds$$' class='display'><thead style='display:none;'><tr>" . "<th></th>" x $col . "</tr></thead><tbody>";
# my $library = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'library'");# ORDER BY name
# $library->execute();
# if ($library->rows > 0)
# {
# 	my $libraryResult;
# 	while (my @library=$library->fetchrow_array())
# 	{
# 		next if ($library[0] eq $assembly[4]);
# 		@{$libraryResult->{$library[2]}} = @library;
# 	}
# 	foreach (sort {uc ($a) cmp uc($b)} keys %$libraryResult)
# 	{
# 		my @library = @{$libraryResult->{$_}};
# 		my $checked = (exists $checkedExtraId->{$library[0]}) ? "checked='checked'" : "";
# 		if($colCount % $col == 0)
# 		{
# 			$assemblyExtraIds .= "<tr><td><input type='checkbox' id='libraryList$library[0]$$' name='extraId' value='$library[0]' $checked><label for='libraryList$library[0]$$' title='library'>$library[2]<sup class='ui-state-disabled'>L</sup></label></td>";
# 		}
# 		elsif($colCount % $col == $col - 1)
# 		{
# 			$assemblyExtraIds .= "<td><input type='checkbox' id='libraryList$library[0]$$' name='extraId' value='$library[0]' $checked><label for='libraryList$library[0]$$' title='library'>$library[2]<sup class='ui-state-disabled'>L</sup></label></td></tr>";
# 		}
# 		else
# 		{
# 			$assemblyExtraIds .= "<td><input type='checkbox' id='libraryList$library[0]$$' name='extraId' value='$library[0]' $checked><label for='libraryList$library[0]$$' title='library'>$library[2]<sup class='ui-state-disabled'>L</sup></label></td>";
# 		}
# 		$colCount++;
# 	}
# }
my $doNotChangeGenome = ($assembly[7] != 0 && $assembly[5] > 0) ? '<sup class="ui-state-error ui-corner-all">Do NOT Change Unless You Plan To Re-Run Assembly!</sup><br>' : '';
my $refGenomeId = '';
my $genome = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'genome'");# ORDER BY name
$genome->execute();
if ($genome->rows > 0)
{
	my $genomeResult;
	while (my @genome=$genome->fetchrow_array())
	{
		next if ($genome[0] eq $assembly[4]);
		@{$genomeResult->{$genome[2]}} = @genome;
	}
	foreach (sort {uc ($a) cmp uc($b)} keys %$genomeResult)
	{
		my @genome = @{$genomeResult->{$_}};
		if ($genome[4] > 0) #for assembly
		{
			my $checked = (exists $checkedExtraId->{$genome[0]}) ? "checked='checked'" : "";
			if($colCount % $col == 0)
			{
				$assemblyExtraIds .= "<tr><td><input type='checkbox' id='genomeList$genome[0]$$' name='extraId' value='$genome[0]' $checked><label for='genomeList$genome[0]$$' title='genome'>$genome[2]<sup class='ui-state-disabled'>G</sup></label></td>";
			}
			elsif($colCount % $col ==  $col - 1)
			{
				$assemblyExtraIds .= "<td><input type='checkbox' id='genomeList$genome[0]$$' name='extraId' value='$genome[0]' $checked><label for='genomeList$genome[0]$$' title='genome'>$genome[2]<sup class='ui-state-disabled'>G</sup></label></td></tr>";
			}
			else
			{
				$assemblyExtraIds .= "<td><input type='checkbox' id='genomeList$genome[0]$$' name='extraId' value='$genome[0]' $checked><label for='genomeList$genome[0]$$' title='genome'>$genome[2]<sup class='ui-state-disabled'>G</sup></label></td>";
			}
			$colCount++;
		}
		if ($genome[5] > 0) #as reference
		{
			$refGenomeId .= ($genome[0] eq $assembly[5] ) ?
				"<option value='$genome[0]' selected>$genome[2]</option>" :
				"<option value='$genome[0]'>$genome[2]</option>";
		}
	}
}
if($refGenomeId)
{
	$refGenomeId = "<option class='ui-state-error-text' value='0'>Please select a genome</option>".$refGenomeId;
}
else
{
	$refGenomeId = "<option class='ui-state-error-text' value='0'>No reference genome available</option>".$refGenomeId;
}

my $toBeFilled = $col - ( $colCount % $col);
$assemblyExtraIds .= ($toBeFilled < $col ) ? "<td>&nbsp;</td>" x $toBeFilled ."</tr></tbody></table>" : "</tbody></table>";

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$assemblyId/$assemblyId/g;
$html =~ s/\$assemblyName/$assembly[2]/g;
$html =~ s/\$assemblyVersion/$assembly[3]/g;
$html =~ s/\$doNotChangeFpc/$doNotChangeFpc/g;
$html =~ s/\$fpcOrAgpId/$fpcOrAgpId/g;
$html =~ s/\$doNotChangeGenome/$doNotChangeGenome/g;
$html =~ s/\$refGenomeId/$refGenomeId/g;
$html =~ s/\$assemblyExtraIds/$assemblyExtraIds/g;
my $assemblyDetails = decode_json $assembly[8];
$assemblyDetails->{'description'} = '' if (!exists $assemblyDetails->{'description'});
$html =~ s/\$assemblyDescription/$assemblyDetails->{'description'}/g;
$assemblyDetails->{'autoCheckNewSeq'} = 0 if (!exists $assemblyDetails->{'autoCheckNewSeq'});
if($assemblyDetails->{'autoCheckNewSeq'})
{
	$html =~ s/\$autoCheck/ checked="checked"/g;
}
else
{
	$html =~ s/\$autoCheck//g;
}
$html =~ s/\$assemblyCreator/$assembly[9]/g;
my $assemblyStatusRadios = <<END;
		<div id="editAssemblyStatus">
			<input type="radio" id="editAssemblyStatusRadio1" name="status" value="0" checked="checked"><label for="editAssemblyStatusRadio1">$assemblyStatus->{$assembly[7]} (no change)</label>
			<input type="radio" id="editAssemblyStatusRadio2" name="status" value="1"><label for="editAssemblyStatusRadio2">$assemblyStatus->{1}</label>
			<input type="radio" id="editAssemblyStatusRadio3" name="status" value="2"><label for="editAssemblyStatusRadio3">$assemblyStatus->{2}</label>
		</div>
END
if ($assembly[7] == 1)
{
	$assemblyStatusRadios = <<END;
		<div id="editAssemblyStatus">
			<input type="radio" id="editAssemblyStatusRadio2" name="status" value="1" checked="checked"><label for="editAssemblyStatusRadio2">$assemblyStatus->{1}</label>
			<input type="radio" id="editAssemblyStatusRadio3" name="status" value="2"><label for="editAssemblyStatusRadio3">$assemblyStatus->{2}</label>
		</div>
END
}
elsif ($assembly[7] == 2)
{
	$assemblyStatusRadios = <<END;
		<div id="editAssemblyStatus">
			<input type="radio" id="editAssemblyStatusRadio2" name="status" value="1"><label for="editAssemblyStatusRadio2">$assemblyStatus->{1}</label>
			<input type="radio" id="editAssemblyStatusRadio3" name="status" value="2" checked="checked"><label for="editAssemblyStatusRadio3">$assemblyStatus->{2}</label>
		</div>
END
}
$html =~ s/\$assemblyStatus/$assemblyStatusRadios/g;
$html =~ s/\$assemblyCreationDate/$assembly[10]/g;
print header;
print $html;

__DATA__
	<form id="editAssembly$assemblyId" name="editAssembly$assemblyId" action="assemblySave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="assemblyId" id="editAssemblyId" type="hidden" value="$assemblyId" />
	<table>
	<tr><td style='text-align:right'><label for="editAssemblyName"><b>Assembly Name</b></label></td><td><input class='ui-widget-content ui-corner-all' name="name" id="editAssemblyName" size="30" type="text" maxlength="32" value="$assemblyName"/><br>Version $assemblyVersion <sup class='ui-state-disabled'>by $assemblyCreator on $assemblyCreationDate</sup><br>
	<input type="checkbox" id="editAssemblyAutoCheckNewSeq" name="autoCheckNewSeq" value="1"$autoCheck><label for="editAssemblyAutoCheckNewSeq">AutoCheck New Sequences</label>
	</td></tr>
	<tr><td style='text-align:right'><label for="editAssemblyStatus"><b>Status</b></label></td><td>
$assemblyStatus
	</td></tr>
	<tr><td style='text-align:right'><label for='editAssemblyFpcOrAgp'><b>Physical Reference</b></label></td><td>$doNotChangeFpc<select class='ui-widget-content ui-corner-all' name='fpcOrAgpId' id='editAssemblyFpcOrAgp'>$fpcOrAgpId</select></td></tr>
	<tr><td style='text-align:right'><label for='editAssemblyRefGenome'><b>Reference Genome</b></label></td><td>$doNotChangeGenome<select class='ui-widget-content ui-corner-all' name='refGenomeId' id='editAssemblyRefGenome'>$refGenomeId</select></td></tr>
	<tr><td style='text-align:right'><label for='editAssemblyExtraGenome'><b>Extra Genome</b></label><br><sup class='ui-state-disabled'>(Gap fillers)</sup></td><td>$assemblyExtraIds</td></tr>
	<tr><td style='text-align:right'><label for="editAssemblyDescription"><b>Description</b></label></td><td><textarea class='ui-widget-content ui-corner-all' name="description" id="editAssemblyDescription" cols="50" rows="10" placeholder="Give some information about this assembly. Or you may do it later.">$assemblyDescription</textarea></td></tr>
	</table>
	</form>
<script>
$('#dialog').dialog("option", "title", "Edit Assembly");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('editAssembly$assemblyId'); } }, { text: "Delete", click: function() { deleteItem($assemblyId); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
$( "#editAssemblyStatus" ).buttonset();
</script>