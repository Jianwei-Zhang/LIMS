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

my $assemblyId = param ('assemblyId') || '';
my $assemblyCtgId = param ('assemblyCtgId') || '';
my $ctgName = '';
if($assemblyCtgId)
{
	my $assemblyCtgList=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$assemblyCtgList->execute($assemblyCtgId);
	my @assemblyCtgList = $assemblyCtgList->fetchrow_array();
	$ctgName = $assemblyCtgList[2];
	$assemblyId = $assemblyCtgList[3];
}
my $assembly=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$assembly->execute($assemblyId);
my @assembly = $assembly->fetchrow_array();

my $autoAssemblyCtgSearchUrl = 'autoAssemblyCtgSearch.cgi?assemblyId=$assemblyId';

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$autoAssemblyCtgSearchUrl/$autoAssemblyCtgSearchUrl/g;
$html =~ s/\$assemblyId/$assemblyId/g;
$html =~ s/\$assemblyName/$assembly[2]/g;
$html =~ s/\$ctgName/$ctgName/g;

print header;
print $html;

__DATA__
<form id="assemblyGapFill" name="assemblyGapFill" action="assemblyGapFill.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="assemblyId" id="assemblyGapFillAssemblyId" type="hidden" value="$assemblyId" />
	<table>
	<tr><td></td><td><label for="assemblyGapFillCtgOne"><b>Contig-1.</b></label><input name="ctgOne" id="assemblyGapFillCtgOne" size="16" type="text" maxlength="16" VALUE="$ctgName" placeholder="First Contig" /></td><td><label for="assemblyGapFillCtgTwo"><b>Contig-2.</b></label><input name="ctgTwo" id="assemblyGapFillCtgTwo" size="16" type="text" maxlength="16" VALUE="" placeholder="Second Contig" /></td></tr>
	</table>
	<sup class='ui-state-highlight'>Leave both Contig-1 and Contig-2 blank to run gap filler for the whole assembly.</sup>
</form>
<script>
$('#dialog').dialog("option", "title", "Gap Filler for Assembly $assemblyName");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Run", click: function() { loadingShow();submitForm('assemblyGapFill'); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
$( "#assemblyGapFillCtgOne" ).autocomplete({
	source: "$autoAssemblyCtgSearchUrl",
	minLength: 1
	});
$( "#assemblyGapFillCtgTwo" ).autocomplete({
	source: "$autoAssemblyCtgSearchUrl",
	minLength: 1
	});
</script>