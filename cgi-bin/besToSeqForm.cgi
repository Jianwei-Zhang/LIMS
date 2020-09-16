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
my $BESTOSEQMINOVERLAP = $userConfig->getFieldValueWithUserIdAndFieldName($userId,"BESTOSEQMINOVERLAP");
my $BESTOSEQIDENTITY = $userConfig->getFieldValueWithUserIdAndFieldName($userId,"BESTOSEQIDENTITY");

undef $/;# enable slurp mode
my $html = <DATA>;

my $libraryId = param ('libraryId') || '';

my $library=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$library->execute($libraryId);
my @library = $library->fetchrow_array();

my $fixedTargetId = param ('targetId') || '';

my $targetId = "<option value='$libraryId'>$library[2] BAC Sequences</option>";

my $genomeList=$dbh->prepare("SELECT * FROM matrix WHERE (z = ? OR id = ?) AND container LIKE 'genome'");
$genomeList->execute($libraryId, $fixedTargetId);
while (my @genomeList = $genomeList->fetchrow_array())
{
	$targetId .= ($genomeList[0] eq $fixedTargetId) ?
		"<option value='$genomeList[0]' selected>$genomeList[2]</option>" :
		"<option value='$genomeList[0]'>$genomeList[2]</option>";
}

$html =~ s/\$libraryId/$libraryId/g;
$html =~ s/\$libraryName/$library[2]/g;
$html =~ s/\$targetId/$targetId/g;
$html =~ s/\$BESTOSEQMINOVERLAP/$BESTOSEQMINOVERLAP/g;
$html =~ s/\$BESTOSEQIDENTITY/$BESTOSEQIDENTITY/g;


print header;
print $html;

__DATA__
<form id="besToSeq" name="besToSeq" action="besToSeq.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<table>
	<tr><td style='text-align:left' colspan='2'><input name="libraryId" id="newLibraryId" type="hidden" value="$libraryId" />Align $libraryName BES to <select class='ui-widget-content ui-corner-all' name='targetId' id='newTargetId'>$targetId</select></td></tr>
	<tr><td></td><td><input type="checkbox" id="newBesToSeqRedo" name="redo" value="1" checked="checked"><label for="newBesToSeqRedo">Override Existing Alignments</label></td></tr>
	<tr><td style='text-align:right'><label for="besToSeqMinOverlap"><b>Minimum Overlap</b></label><br>(length in bp)</td><td><input name="minOverlapBesToSeq" id="besToSeqMinOverlap" size="4" type="text" maxlength="6" VALUE="$BESTOSEQMINOVERLAP" /></td></tr>
	<tr><td style='text-align:right'><label for="besToSeqIdentity"><b>Minimum Identity</b></label><br>(%)</td><td><input name="identityBesToSeq" id="besToSeqIdentity" size="4" type="text" maxlength="4" VALUE="$BESTOSEQIDENTITY" /></td></tr>
	</table>
</form>
<script>
$('#dialog').dialog("option", "title", "Run BES to SEQ Alignment");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Run Alignment", click: function() { submitForm('besToSeq'); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
$( "#besToSeqMinOverlap" ).spinner({
	min: 200,
	max: 1000,
	step: 50,
	start: function( event, ui ) {
		var current = $( "#besToSeqMinOverlap" ).spinner( "value" );
		if(current >= 500)
		{
			$( "#besToSeqMinOverlap" ).spinner( "option", "step", 100 );
		}
		else
		{
			$( "#besToSeqMinOverlap" ).spinner( "option", "step", 50 );
		}
	}
});
$( "#besToSeqIdentity" ).spinner({ min: 90, max: 100});
</script>