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

my $itemId = param ('itemId') || '';
my $item=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$item->execute($itemId);
my @item = $item->fetchrow_array();

my %commentType = (
	0=>'General',
	1=>'Alert',
	2=>'Favorite'
	);
my $commentType = '';
foreach (sort {$a <=> $b} keys %commentType)
{
	$commentType .= "<option value='$_'>$commentType{$_}</option>";
}

print header;
my $comment = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'comment' AND o = ?");
$comment->execute($itemId);
if ($comment->rows > 0)
{

	my @comment=$comment->fetchrow_array();
	my $commentDetails = decode_json $comment[8];
	$commentDetails->{'description'} = '' unless (exists $commentDetails->{'description'});
	$commentDetails->{'description'} = escapeHTML($commentDetails->{'description'});
	$commentDetails->{'description'} =~ s/\n/<br>/g;

	print <<END;
	<table>
		<tr><td><b>$commentType{$comment[4]} Comments on $item[1] - $item[2]</b> <sup class='ui-state-disabled'>by $comment[9] on $comment[10]</sup></td></tr>
		<tr><td>$commentDetails->{'description'}</td></tr>
	</table>
	<script>
	\$('#dialog').dialog("option", "title", "View Comment");
	\$( "#dialog" ).dialog( "option", "buttons", [{ text: "Edit", click: function() { closeDialog();openDialog("commentEdit.cgi?itemId=$itemId"); } }, { text: "OK", click: function() {closeDialog(); } } ] );
	</script>
END
}
else
{
	my $noun = 'word';
	print <<END;
	<form id="newComment" name="newComment" action="commentSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
		<input name="itemId" id="newCommentItemId" type="hidden" value="$itemId" />
		<table>
		<tr><td style='text-align:left;white-space: nowrap;' colspan='2'><label for="newCommentDescription"><b>Comments on $item[1] - $item[2]</b></label></td></tr>
		<tr><td><textarea class='ui-widget-content ui-corner-all word_count' name="description" id="newCommentDescription" cols="50" rows="10"></textarea></td>
		<td><label for='newCommentType'>Comment type</label><select class='ui-widget-content ui-corner-all' name='commentType' id='newCommentType'>$commentType</select><br><br>
			<sub id='newCommentDescription_count' style='display:none'></sub>
		</td></tr>
		</table>
	</form>
	<script>
	\$('#dialog').dialog("option", "title", "New Comment");
	\$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('newComment'); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
	wordCount('$noun');
	</script>
END
}
