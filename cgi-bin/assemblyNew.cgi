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

my $objectComponent;
$objectComponent->{0} = "Unknown";
$objectComponent->{1} = "Chr-Seq";
$objectComponent->{2} = "Ctg-Seq";
my $targetId = param ('targetId') || '';
my $target=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$target->execute($targetId);
my @target = $target->fetchrow_array();

my $fpcOrAgpId = '';
if($target[1] eq "library")
{
	my $fpcList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'fpc' AND barcode = ? ");
	$fpcList->execute($targetId);
	while (my @fpcList = $fpcList->fetchrow_array())
	{
		$fpcOrAgpId .= "<option value='$fpcList[0]' title='$fpcList[8]'>FPC: $fpcList[2] v.$fpcList[3]</option>";
	}
}
elsif($target[1] eq 'genome')
{
	my $agpList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'agp' AND x = ?");
	$agpList->execute($targetId);
	while (my @agpList = $agpList->fetchrow_array())
	{
		$fpcOrAgpId .= "<option value='$agpList[0]' title='$objectComponent->{$agpList[5]} AGP'>AGP: $agpList[2] v.$agpList[3]</option>";
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

my $col = 3;
my $colCount=0;
my $assemblyExtraIds = "<table id='assemblyExtraIds$$' class='display'><thead style='display:none;'><tr>" . "<th></th>" x $col . "</tr></thead><tbody>";
# my $library = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'library'");# ORDER BY name
# $library->execute();
# if($library->rows > 0)
# {
# 	my $libraryResult;
# 	while (my @library=$library->fetchrow_array())
# 	{
# 		next if ($library[0] eq $targetId);
# 		@{$libraryResult->{$library[2]}} = @library;
# 	}
# 	foreach (sort {uc ($a) cmp uc($b)} keys %$libraryResult)
# 	{
# 		my @library = @{$libraryResult->{$_}};
# 		if($colCount % $col == 0)
# 		{
# 			$assemblyExtraIds .= "<tr><td><input type='checkbox' id='libraryList$library[0]$$' name='extraId' value='$library[0]'><label for='libraryList$library[0]$$' title='library'>$library[2]<sup class='ui-state-disabled'>L</sup></label></td>";
# 		}
# 		elsif($colCount % $col == $col - 1)
# 		{
# 			$assemblyExtraIds .= "<td><input type='checkbox' id='libraryList$library[0]$$' name='extraId' value='$library[0]'><label for='libraryList$library[0]$$' title='library'>$library[2]<sup class='ui-state-disabled'>L</sup></label></td></tr>";
# 		}
# 		else
# 		{
# 			$assemblyExtraIds .= "<td><input type='checkbox' id='libraryList$library[0]$$' name='extraId' value='$library[0]'><label for='libraryList$library[0]$$' title='library'>$library[2]<sup class='ui-state-disabled'>L</sup></label></td>";
# 		}
# 		$colCount++;
# 	}
# }
my $refGenomeId = '';
my $genome = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'genome'");# ORDER BY name
$genome->execute();
if($genome->rows > 0)
{
	my $genomeResult;
	while (my @genome=$genome->fetchrow_array())
	{
		next if ($genome[0] eq $targetId);
		@{$genomeResult->{$genome[2]}} = @genome;
	}
	foreach (sort {uc ($a) cmp uc($b)} keys %$genomeResult)
	{
		my @genome = @{$genomeResult->{$_}};
		if ($genome[4] > 0) # for assembly
		{
			if($colCount % $col == 0)
			{
				$assemblyExtraIds .= "<tr><td><input type='checkbox' id='genomeList$genome[0]$$' name='extraId' value='$genome[0]'><label for='genomeList$genome[0]$$' title='genome'>$genome[2]<sup class='ui-state-disabled'>G</sup></label></td>";
			}
			elsif($colCount % $col ==  $col - 1)
			{
				$assemblyExtraIds .= "<td><input type='checkbox' id='genomeList$genome[0]$$' name='extraId' value='$genome[0]'><label for='genomeList$genome[0]$$' title='genome'>$genome[2]<sup class='ui-state-disabled'>G</sup></label></td></tr>";
			}
			else
			{
				$assemblyExtraIds .= "<td><input type='checkbox' id='genomeList$genome[0]$$' name='extraId' value='$genome[0]'><label for='genomeList$genome[0]$$' title='genome'>$genome[2]<sup class='ui-state-disabled'>G</sup></label></td>";
			}
			$colCount++;
		}
		if ($genome[5] > 0) # as reference
		{
			$refGenomeId .= "<option value='$genome[0]' title='$genome[8]'>$genome[2]</option>";
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
$html =~ s/\$targetId/$targetId/g;
$html =~ s/\$assemblyName/$target[2]/g;
$html =~ s/\$fpcOrAgpId/$fpcOrAgpId/g;
$html =~ s/\$refGenomeId/$refGenomeId/g;
$html =~ s/\$assemblyExtraIds/$assemblyExtraIds/g;

print header;
print $html;

__DATA__
	<form id="newAssembly$targetId" name="newAssembly$targetId" action="assemblySave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="targetId" id="targetId" type="hidden" value="$targetId" />
	<table>
	<tr><td style='text-align:right'><label for="newAssemblyName"><b>Assembly Name</b></label></td><td><input class='ui-widget-content ui-corner-all' name="name" id="newAssemblyName" size="40" type="text" maxlength="32" value="$assemblyName"/> <input type="checkbox" id="newAssemblyAutoCheckNewSeq" name="autoCheckNewSeq" value="1"><label for="newAssemblyAutoCheckNewSeq">AutoCheck New Sequences</label></td></tr>
	<tr><td style='text-align:right'><label for='newAssemblyFpcOrAgp'><b>Physical Reference</b></label></td><td><select class='ui-widget-content ui-corner-all' name='fpcOrAgpId' id='newAssemblyFpcOrAgp'>$fpcOrAgpId</select></td></tr>
	<tr><td style='text-align:right'><label for='newAssemblyRefGenome'><b>Reference Genome</b></label></td><td><select class='ui-widget-content ui-corner-all' name='refGenomeId' id='newAssemblyRefGenome'>$refGenomeId</select></td></tr>
	<tr><td style='text-align:right'><label for='newAssemblyExtraGenome'><b>Extra Genome</b></label><br><sup class='ui-state-disabled'>(Gap fillers)</sup></td><td>$assemblyExtraIds</td></tr>
	<tr><td style='text-align:right'><label for="newAssemblyDescription"><b>Description</b></label></td><td><textarea class='ui-widget-content ui-corner-all' name="description" id="newAssemblyDescription" cols="60" rows="10" placeholder="Give some information about this assembly. Or you may do it later."></textarea></td></tr>
	<tr><td></td><td><INPUT TYPE="button" VALUE="Save" onclick="submitForm('newAssembly$targetId');"></td></tr>
	</table>
	</form>
<script>
buttonInit();
loadingHide();
</script>