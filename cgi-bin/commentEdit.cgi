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

my $noun = 'word';
my $itemId = param ('itemId') || '';
my $item=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$item->execute($itemId);
my @item = $item->fetchrow_array();
my $comment = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'comment' AND o = ?");
$comment->execute($itemId);
my @comment=$comment->fetchrow_array();
my $commentDetails = decode_json $comment[8];
$commentDetails->{'description'} = '' unless (exists $commentDetails->{'description'});

$html =~ s/\$commentId/$comment[0]/g;
$html =~ s/\$itemId/$item[0]/g;
$html =~ s/\$itemContainer/$item[1]/g;
$html =~ s/\$itemName/$item[2]/g;
my %commentType = (
	0=>'General',
	1=>'Alert',
	2=>'Favorite'
	);
my $commentType = '';
foreach (sort {$a <=> $b} keys %commentType)
{
	$commentType .= ($comment[4] == $_) ? "<option value='$_' selected>$commentType{$_}</option>" : "<option value='$_'>$commentType{$_}</option>";
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$commentType/$commentType/g;
$html =~ s/\$commentDescription/$commentDetails->{'description'}/g;
$html =~ s/\$noun/$noun/g;

print header;
print $html;

__DATA__
<form id="editComment" name="editComment" action="commentSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="itemId" id="editCommentItemId" type="hidden" value="$itemId" />
	<table>
	<tr><td style='text-align:left;white-space: nowrap;' colspan='2'><label for="editCommentDescription"><b>Comments on $itemContainer - $itemName</b></label></td></tr>
	<tr><td><textarea class='ui-widget-content ui-corner-all word_count' name="description" id="editCommentDescription" cols="50" rows="10">$commentDescription</textarea></td>
	<td><label for='editCommentType'>Comment type</label><select class='ui-widget-content ui-corner-all' name='commentType' id='editCommentType'>$commentType</select><br><br>
		<sub id='editCommentDescription_count' style='display:none'></sub></td></tr>
	</table>
</form>
<script>
$('#dialog').dialog("option", "title", "Edit Comment");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('editComment'); } }, { text: "Delete", click: function() { deleteItem($commentId); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
wordCount('$noun');
</script>