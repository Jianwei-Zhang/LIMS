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

print header;

my $genomeId = param ('itemId') || '';
unless($genomeId)
{
	print <<END;
<script>
	parent.closeDialog();
	parent.errorPop("Please select a source genome!");
</script>	
END
	exit;
}

my $genome = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$genome->execute($genomeId);
my @genome=$genome->fetchrow_array();

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$genomeId/$genomeId/g;
$html =~ s/\$genomeName/$genome[2]/g;

if($genome[4])
{
	$html =~ s/\$splitGenomeForAssembly/ checked="checked"/g;
}
else
{
	$html =~ s/\$splitGenomeForAssembly//g;
}
if($genome[5])
{
	$html =~ s/\$splitGenomeAsReference/ checked="checked"/g;
}
else
{
	$html =~ s/\$splitGenomeAsReference//g;
}

my $libraryId = "<option class='ui-state-error-text' value='0'>None</option>";
my $libraryList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'library'");# ORDER BY name
$libraryList->execute();
if($libraryList->rows > 0)
{
	my $libraryListResult;
	while (my @libraryList = $libraryList->fetchrow_array())
	{
		@{$libraryListResult->{$libraryList[2]}} = @libraryList;
	}
	foreach (sort {uc ($a) cmp uc($b)} keys %$libraryListResult)
	{
		my @libraryList = @{$libraryListResult->{$_}};
		my $libraryDetails = decode_json $libraryList[8];
		$libraryDetails->{'comments'} = escapeHTML($libraryDetails->{'comments'});
		$libraryId .= ($libraryList[0] eq $genome[6] ) ?
			"<option value='$libraryList[0]' title='$libraryDetails->{'comments'}' selected>Library: $libraryList[2]</option>" :
			"<option value='$libraryList[0]' title='$libraryDetails->{'comments'}'>Library: $libraryList[2]</option>";
	}
}

$html =~ s/\$libraryId/$libraryId/g;
$html =~ s/\$genomeDescription/$genome[8]/g;

my $sequences = "
	<table id='sequences$$' class='display'>
		<thead>
			<tr>
				<th>
					<input type='checkbox' id='checkAllBox$$' name='checkAllBox$$' value='Check all' checked='checked' onClick='checkAll(\"seqId\");return false;' title='Check all'>
					<input type='checkbox' id='uncheckAllBox$$' name='uncheckAllBox$$' value='Uncheck all' onClick='uncheckAll(\"seqId\");return false;' title='Uncheck all'>
				</th>
				<th><b>Seq Name</b></th>
				<th><b>Length (bp)</b></th>
				<th>Description</th>
			</tr>
		</thead>
		<tbody>";
my $getSequences = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 99 AND x = ? ORDER BY y DESC");
$getSequences->execute($genomeId);
while (my @getSequences =  $getSequences->fetchrow_array())
{
	$getSequences[5] = commify($getSequences[5]);
	my $sequenceDetails = decode_json $getSequences[8];
	$sequenceDetails->{'id'} = '' unless (exists $sequenceDetails->{'id'});
	$sequenceDetails->{'description'} = '' unless (exists $sequenceDetails->{'description'});
	my $descriptionLength = length ($sequenceDetails->{'description'});
	if ($descriptionLength > 35)
	{
		$sequenceDetails->{'description'} = "<a title='$sequenceDetails->{'description'}'>". substr($sequenceDetails->{'description'},0,25). "...". substr($sequenceDetails->{'description'},-5). "</a>";
	}
	$sequences .= "<tr onclick='checkOrUncheck(\"seqList$getSequences[0]$$\");'><td style='text-align:center;'><input type='checkbox' id='seqList$getSequences[0]$$' name='seqId' value='$getSequences[0]' onclick='checkOrUncheck(\"seqList$getSequences[0]$$\");'></td><td>$getSequences[2]</td><td>$getSequences[5]</td><td>$sequenceDetails->{'description'}</td></tr>";
}
$sequences .= "</tbody></table>\n";

$html =~ s/\$sequences/$sequences/g;
$html =~ s/\$\$/$$/g;
print $html;

__DATA__
<form id="splitGenome" name="splitGenome" action="genomeSplit.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="genomeId" id="splitGenomeId" type="hidden" value="$genomeId" />
	<b>Select Sequences in $genomeName</b><br>
	$sequences<br>
	<label for="splitGenomeName"><b>Split to</b></label><input class='ui-widget-content ui-corner-all' name="name" id="splitGenomeName" size="40" type="text" maxlength="32" value="$genomeName-Split" placeholder="New Genome Name"/><br>
	<input type="checkbox" id="splitGenomeForAssembly" name="forAssembly" value="1"$splitGenomeForAssembly><label for="splitGenomeForAssembly">Enable this genome to be reassembled with GPM</label><br>
	<input type="checkbox" id="splitGenomeAsReference" name="asReference" value="1"$splitGenomeAsReference><label for="splitGenomeAsReference">Enable this genome to be used as a reference</label><br><hr width="80%">
	<label for="splitGenomeLibraryId"><b>Link to</b></label><select class='ui-widget-content ui-corner-all' name='libraryId' id='splitGenomeLibraryId'>$libraryId</select><br>
	<label for="splitGenomeDescription"><b>Description</b></label><br>
	<textarea class='ui-widget-content ui-corner-all' name="description" id="splitGenomeDescription" cols="80" rows="10">$genomeDescription</textarea>
</form>
<script>
$('#dialog').dialog("option", "title", "Split Genome");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Split", click: function() { submitForm('splitGenome'); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
$( "#sequences$$" ).dataTable({
	"scrollY": "400px",
	"scrollCollapse": true,
	"paging": false,
	"columnDefs": [
    { "orderable": false, "targets": 0 }
  ]
});
</script>
