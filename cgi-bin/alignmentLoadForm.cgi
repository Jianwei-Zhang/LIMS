#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use DBI;
use lib "lib/";
use lib "lib/pangu";
use pangu;
use userConfig;
use userCookie;

my $userCookie = new userCookie;
my $userId = (cookie('cid')) ? $userCookie->checkCookie(cookie('cid')) : 0;
exit if (!$userId);

my $commoncfg = readConfig("main.conf");
my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});

my $queryId = param ('queryId') || '';
my $subjectId = param ('subjectId') || '';

my $queryGenomeId = '';
my $subjectGenomeId = '';

my $assemblyId = param ('assemblyId') || '';
if ($assemblyId)
{
	my $assembly=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$assembly->execute($assemblyId);
	my @assembly = $assembly->fetchrow_array();
	my $target=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$target->execute($assembly[4]);
	my @target = $target->fetchrow_array();
	my $extraGenomeId = '';
	my $checkAsbGenome = $dbh->prepare("SELECT child FROM link WHERE parent = ? AND type LIKE 'asbGenome'");
	$checkAsbGenome->execute($assemblyId);
	my $checkedExtraId;
	while(my @checkAsbGenome=$checkAsbGenome->fetchrow_array())
	{
		$checkedExtraId->{$checkAsbGenome[0]} = 1;
		my $target=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$target->execute($checkAsbGenome[0]);
		my @target = $target->fetchrow_array();
		$extraGenomeId .= "<option value='$checkAsbGenome[0]' title='Extra data for assembly'>[Extra] $target[2] $target[1]</option>";
	}

	my $allGenomeResult;
	my $genomeList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'genome'");
	$genomeList->execute();
	while (my @genomeList = $genomeList->fetchrow_array())
	{
		next if ($genomeList[0] eq $assembly[4]);
		next if ($genomeList[5] < 1);#remove not for reference
		next if (exists $checkedExtraId->{$genomeList[0]}); #skip extra genomes
		@{$allGenomeResult->{$genomeList[2]}} = @genomeList;
	}
	foreach (sort {uc ($a) cmp uc($b)} keys %$allGenomeResult)
	{
		my @genomeList = @{$allGenomeResult->{$_}};
		if($genomeList[0] eq $assembly[5])
		{
			$subjectGenomeId = "<option value='$genomeList[0]' title='$genomeList[8]' selected>[Assembly Default] $genomeList[2]</option>" . $subjectGenomeId;
		}
		else
		{
			$subjectGenomeId .= "<option value='$genomeList[0]' title='$genomeList[8]'>$genomeList[2]</option>";
		}
	}
	$queryGenomeId = <<END;
	<tr><td style='text-align:right'><label for="newAlignmentQuery"><b>Query Genome</b></label></td>
	<td><input name="queryGenomeId" id="newAlignmentQuery" type="hidden" value="$assembly[4]" />$target[2] ($target[1])</td>
	</tr>	
END
	$subjectGenomeId = <<END;
	<tr><td style='text-align:right'><label for="newAlignmentSubject"><b>Subject Genome</b></label></td>
	<td><select class='ui-widget-content ui-corner-all' name="subjectGenomeId" id="newAlignmentSubject">
	<option value='$assembly[4]' title='SeqToSeq Alignment'>[Self Alignment] $target[2]</option>
	$extraGenomeId
	$subjectGenomeId
	</select></td>
	</tr>
END
}
else
{
	my $allGenomeResult;
	my $genomeList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'genome'");
	$genomeList->execute();
	while (my @genomeList = $genomeList->fetchrow_array())
	{
		@{$allGenomeResult->{$genomeList[2]}} = @genomeList;
	}
	foreach (sort {uc ($a) cmp uc($b)} keys %$allGenomeResult)
	{
		my @genomeList = @{$allGenomeResult->{$_}};
		$queryGenomeId .= ($genomeList[0] eq $queryId ) ?
			"<option value='$genomeList[0]' title='$genomeList[8]' selected>$genomeList[2]</option>" :
			"<option value='$genomeList[0]' title='$genomeList[8]'>$genomeList[2]</option>";
		$subjectGenomeId .= ($genomeList[0] eq $subjectId ) ?
			"<option value='$genomeList[0]' title='$genomeList[8]' selected>$genomeList[2]</option>" :
			"<option value='$genomeList[0]' title='$genomeList[8]'>$genomeList[2]</option>";
	}
	$queryGenomeId = <<END;
	<tr><td style='text-align:right'><label for="newAlignmentQuery"><b>Query Genome</b></label></td>
	<td><select class='ui-widget-content ui-corner-all' name="queryGenomeId" id="newAlignmentQuery">$queryGenomeId</select></td>
	</tr>	
END
	$subjectGenomeId = <<END;
	<tr><td style='text-align:right'><label for="newAlignmentSubject"><b>Subject Genome</b></label></td>
	<td><select class='ui-widget-content ui-corner-all' name="subjectGenomeId" id="newAlignmentSubject">$subjectGenomeId</select></td>
	</tr>
END
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$queryGenomeId/$queryGenomeId/g;
$html =~ s/\$subjectGenomeId/$subjectGenomeId/g;

print header;
print $html;

__DATA__
<form id="newAlignment" name="newAlignment" action="alignmentLoad.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<table>
	$queryGenomeId
	$subjectGenomeId
	<tr><td></td><td><input type="checkbox" id="alignmentRedo" name="redo" value="1" checked="checked"><label for="alignmentRedo">Override Existing Alignments</label></td></tr>
	<tr><td style='text-align:right'><label for="alignEngine"><b>Alignment Engine</b></label></td><td style='text-align:left'><select class='ui-widget-content ui-corner-all' name="alignEngine" id="newAssemblyEngine"><option value="Customized">Customized</option><option value="blastn">blastn</option><option value="BLAT">BLAT</option></select></td></tr>
	<tr><td style='text-align:right' rowspan='2'><label for="newAlignmentFile"><b>Alignment File</b></label></td><td><sup class="ui-state-error-text">Note: Sequence Id should contain no more than 32 characters.</sup><br><input name="alignmentFile" id="newAlignmentFile" type="file" />(in BLAST Tabular format)</td></tr>
	<tr><td>or <input name="alignmentFilePath" id="newAlignmentFilePath" type="text" />(On-server file name with full path)</td></tr>
	</table>
</form>
<script>
$('#dialog').dialog("option", "title", "Load Alignment");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Load Alignment", click: function() { submitForm('newAlignment'); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
</script>