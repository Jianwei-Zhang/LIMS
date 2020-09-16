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
my $userConfig = new userConfig;
my $SEQTOGNMMINOVERLAP = $userConfig->getFieldValueWithUserIdAndFieldName($userId,"SEQTOGNMMINOVERLAP");
my $SEQTOGNMIDENTITY = $userConfig->getFieldValueWithUserIdAndFieldName($userId,"SEQTOGNMIDENTITY");

undef $/;# enable slurp mode
my $html = <DATA>;

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
	<td colspan='2'><input name="queryGenomeId" id="newAlignmentQuery" type="hidden" value="$assembly[4]" />$target[2] ($target[1])</td>
	</tr>	
END
	$subjectGenomeId = <<END;
	<tr><td style='text-align:right'><label for="newAlignmentSubject"><b>Subject Genome</b></label></td>
	<td colspan='2'><select class='ui-widget-content ui-corner-all' name="subjectGenomeId" id="newAlignmentSubject">
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
	<td colspan='2'><select class='ui-widget-content ui-corner-all' name="queryGenomeId" id="newAlignmentQuery">$queryGenomeId</select></td>
	</tr>	
END
	$subjectGenomeId = <<END;
	<tr><td style='text-align:right'><label for="newAlignmentSubject"><b>Subject Genome</b></label></td>
	<td colspan='2'><select class='ui-widget-content ui-corner-all' name="subjectGenomeId" id="newAlignmentSubject">$subjectGenomeId</select></td>
	</tr>
END
}

$html =~ s/\$queryGenomeId/$queryGenomeId/g;
$html =~ s/\$subjectGenomeId/$subjectGenomeId/g;
$html =~ s/\$SEQTOGNMMINOVERLAP/$SEQTOGNMMINOVERLAP/g;
$html =~ s/\$SEQTOGNMIDENTITY/$SEQTOGNMIDENTITY/g;

print header;
unless($subjectGenomeId)
{
	print <<END;
<script>
	closeDialog();
	errorPop("No genome available, please load one first!");
</script>	
END
}
print $html;

__DATA__
<form id="newAlignment" name="newAlignment" action="alignment.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<table>
	$queryGenomeId
	$subjectGenomeId
	<tr><td></td><td><input type="checkbox" id="newAlignmentRedo" name="redo" value="1" checked="checked"><label for="newAlignmentRedo">Override Existing Alignments</label></td></tr>
	<tr><td><b>Alignment Settings:</b></td><td style='text-align:left' colspan='2'><label for="alignEngine"><b>Engine</b></label> <select class='ui-widget-content ui-corner-all' name="alignEngine" id="newEngine"><option value="blastn">blastn</option><option value="BLAT">BLAT</option></select><br></td></tr>
	<tr><td></td><td style='text-align:left' colspan='2'>
		<input type="checkbox" id="newMegablastCheckbox" name="megablast" value="megablast" checked="checked"><label for="newMegablastCheckbox">Megablast</label><br>
		<input type="checkbox" id="newSoftMaskingCheckbox" name="softMasking" value="1"><label for="newSoftMaskingCheckbox">Soft Masking (blastn only)</label><br>
		<input type="checkbox" id="newAlignmentSpeedyMode" name="speedyMode" value="1" checked="checked"><label for="newAlignmentSpeedyMode" title="This will skip 'Check End Match'">Speedy Mode</label><br>
		<input type="checkbox" id="newAlignmentCheckGood" name="checkGood" value="1"><label for="newAlignmentCheckGood">Check End Match for SeqToSeq</label><br>
		<input type="checkbox" id="newMarkRepeatRegionCheckbox" name="markRepeatRegion" value="1" checked="checked"><label for="newMarkRepeatRegionCheckbox">Mark Repeat Region</label></td></tr>
	<tr><td></td><td style='text-align:right'><label for="newAlignmentMinOverlap"><b>Minimum Overlap</b></label></td><td><input name="minOverlapAlignment" id="newAlignmentMinOverlap" size="4" type="text" maxlength="6" VALUE="$SEQTOGNMMINOVERLAP" />(bp)</td></tr>
	<tr><td></td><td style='text-align:right'><label for="newAlignmentIdentity"><b>Minimum Identity</b></label></td><td><input name="identityAlignment" id="newAlignmentIdentity" size="4" type="text" maxlength="4" VALUE="$SEQTOGNMIDENTITY" />(%)</td></tr>
	<tr><td></td><td style='text-align:left' colspan='2'><input type="checkbox" id="newEmailNotificationCheckbox" name="emailNotification" value="1" checked="checked"><label for="newEmailNotificationCheckbox">Email a notification after alignment is done.</label></td></tr>
	</table>
</form>
<script>
$('#dialog').dialog("option", "title", "Genome Alignment");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Run Alignment", click: function() { submitForm('newAlignment'); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
$( "#newAlignmentMinOverlap" ).spinner({
	min: 500,
	max: 99999,
	step: 100,
	start: function( event, ui ) {
		var current = $( "#newAlignmentMinOverlap" ).spinner( "value" );
		if(current >= 2000)
		{
			$( "#newAlignmentMinOverlap" ).spinner( "option", "step", 500 );
		}
		else
		{
			$( "#newAlignmentMinOverlap" ).spinner( "option", "step", 100 );
		}
	}
});
$( "#newAlignmentIdentity" ).spinner({ min: 90, max: 100});
</script>