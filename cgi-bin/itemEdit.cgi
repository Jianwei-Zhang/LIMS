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
my $projectId = '';
my $item = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$item->execute($itemId);
my @item=$item->fetchrow_array();
my $itemDetails = decode_json $item[8];

my $parent=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$parent->execute($item[6]);
my @parent = $parent->fetchrow_array();	
my $itemLine = '';
for (sort {$a <=> $b} keys %$itemDetails)
{
	$itemDetails->{$_}->{'field'} = '' unless ($itemDetails->{$_}->{'field'});
	$itemDetails->{$_}->{'value'} = '' unless ($itemDetails->{$_}->{'value'});
	$itemDetails->{$_}->{'value'} = escapeHTML($itemDetails->{$_}->{'value'});
	$itemDetails->{$_}->{'value'} =~ s/\n/<br>/g;
	$itemLine .= "<tr><td style='text-align:right;width:200px;'><b>$itemDetails->{$_}->{'field'}</b></td>
					<td>$itemDetails->{$_}->{'value'}</td>
				</tr>";
}

print header;
print <<END;
<form id="editItem" name="editItem" action="itemSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="itemId" id="editItemId" type="hidden" value="$itemId" />
	<table>
	$itemLine
	</table>
</form>
<script>
\$('#dialog').dialog("option", "title", "Edit $item[1] $item[2] in $parent[1] $parent[2]");
\$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('editItem'); } }, { text: "Delete", click: function() { deleteItem($itemId); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
</script>
END


