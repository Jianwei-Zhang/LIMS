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

my @jobId = param('jobId');
print header;
unless(@jobId)
{
	my $jobList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'job' ORDER BY name");
	$jobList->execute();
	while (my @jobList = $jobList->fetchrow_array())
	{
		next if ($jobList[2] =~ /^-/);
		push @jobId, $jobList[2];
	}
}
my $jobList = 'No Jobs Available or Selected.';
my $jobId;
for(@jobId)
{
	$jobId->{$_} = $_;
}
$jobList = checkbox_group(-name=>'jobId',
	-values=>[@jobId],
	-default=>[@jobId],
	-onclick=>'return false;',
	-labels=>\%{$jobId},
	-columns=>8) if (@jobId);

my $vectorId = "";
my $vectorList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'vector'");
$vectorList->execute();
while (my @vectorList = $vectorList->fetchrow_array())
{
	$vectorId .= "<option value='$vectorList[0]'>$vectorList[2]</option>";
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$jobList/$jobList/g;
$html =~ s/\$vectorId/$vectorId/g;

print $html;

__DATA__
<form id="runPostHGAP" name="runPostHGAP" action="jobPost.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<table>
	<tr><td style='text-align:right;white-space: nowrap;'><b>Circularization Settings:</b></td><td colspan="2"></td></tr>
	<tr><td style='text-align:right' colspan="2"><label for="postHGAPVectorId"><b>Vector(s)</b></label><br><sup>Hold down the control (ctrl for windows) or the <br>command (for Mac) button to select multiple vectors</sup></td><td><select class='ui-widget-content ui-corner-all' name="vectorId" id="postHGAPVectorId" multiple><option value='0' title='Default' selected>(pAGIBAC1_Cut_HindIII)</option>$vectorId</select></td></tr>
	<tr><td></td><td style='text-align:right'><label for="postHGAPMinOverlap"><b>Minimum Overlap</b></label> (length in bp)</td><td><input name="minOverlap" id="postHGAPMinOverlap" size="4" type="text" maxlength="6" VALUE="500" /></td></tr>
	<tr><td></td><td style='text-align:right'><label for="postHGAPIdentity"><b>Minimum Identity</b></label> (%)</td><td><input name="identity" id="postHGAPIdentity" size="4" type="text" maxlength="4" VALUE="95" /></td></tr>
	<tr><td></td><td style='text-align:right'><label for="postHGAPShortCutoff"><b>Filter Sequence Shorter Than </b></label> </td><td><input name="shortCutoff" id="postHGAPShortCutoff" size="4" type="text" maxlength="6" VALUE="10000" /> bp</td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><b>BAC ID Assignment</b></td><td colspan='2'><input type="checkbox" id="runAssignBacGoodOnlyCheckbox" name="goodOnly" value="1"><label for="runAssignBacGoodOnlyCheckbox">Inserted/circularized sequences only, otherwise non-vector and gapped sequences included. </label></td></tr>
	</table>
	<table>
	<tr><td style='text-align:right'><b>Selected Jobs:</b></td><td></td></tr>
	<tr><td style='text-align:right'></td><td>$jobList</td></tr>
	</table>
</form>
<script>
$( "#postHGAPMinOverlap" ).spinner({
	min: 500,
	max: 99999,
	step: 100,
	start: function( event, ui ) {
		var current = $( "#postHGAPMinOverlap" ).spinner( "value" );
		if(current >= 2000)
		{
			$( "#postHGAPMinOverlap" ).spinner( "option", "step", 500 );
		}
		else
		{
			$( "#postHGAPMinOverlap" ).spinner( "option", "step", 100 );
		}
	}
});
$( "#postHGAPIdentity" ).spinner({ min: 90, max: 100});
$('#dialog').dialog("option", "title", "Run postHGAP");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Run", click: function() { submitForm('runPostHGAP'); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
</script>