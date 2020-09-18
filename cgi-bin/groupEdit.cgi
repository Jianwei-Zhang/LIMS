#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use JSON::XS;
use DBI;
use lib "lib/";
use lib "lib/pangu";
use pangu;
use user;
use config;
use userConfig;
use userCookie;

my $userCookie = new userCookie;
my $userId = (cookie('cid')) ? $userCookie->checkCookie(cookie('cid')) : 0;
exit if (!$userId);

my $user = new user;
my $userDetail = $user->getAllFieldsWithUserId($userId);
my $userName = $userDetail->{"userName"};
my $role = $userDetail->{"role"};
exit if ($role ne "admin");

my $commoncfg = readConfig("main.conf");
my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});

my $config = new config;
my $userConfig = new userConfig;

my $groupId = param ('groupId') || '';
my $group = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
$group->execute($groupId);
my @group=$group->fetchrow_array();
my $groupDetails = decode_json $group[8];
$groupDetails->{'permission'} = '' unless (exists $groupDetails->{'permission'});
$groupDetails->{'description'} = '' unless (exists $groupDetails->{'description'});

my $checkedUserId;
my $checkGroup = $dbh->prepare("SELECT child FROM link WHERE parent = ? AND type LIKE 'group'");
$checkGroup->execute($groupId);
while(my @checkGroup=$checkGroup->fetchrow_array())
{
	$checkedUserId->{$checkGroup[0]} = 1;
}
my $checkedPermissionId;
for (split (",",$groupDetails->{'permission'}))
{
	$checkedPermissionId->{$_} = 1;
}

my $colUser = 3;
my $colCountUser=0;
my $groupUserIds = "<table id='groupUserIds$$' class='display'><thead style='display:none;'><tr>" . "<th></th>" x $colUser . "</tr></thead><tbody>";

my $allUsers = $user->getAllFieldsOfUser();
foreach my $individualUser (@$allUsers)
{
	my $individualUserId = $individualUser->{"id"};
	my $userName = $individualUser->{"userName"};
	my $checked = (exists $checkedUserId->{$individualUserId}) ? "checked='checked'" : "";
	if($colCountUser % $colUser == 0)
	{
		$groupUserIds .= "<tr><td><input type='checkbox' id='userList$individualUserId$$' name='userId' value='$individualUserId' $checked><label for='userList$individualUserId$$' title='User'>$userName</label></td>";
	}
	elsif($colCountUser % $colUser == $colUser - 1)
	{
		$groupUserIds .= "<td><input type='checkbox' id='userList$individualUserId$$' name='userId' value='$individualUserId' $checked><label for='userList$individualUserId$$' title='User'>$userName</label></td></tr>";
	}
	else
	{
		$groupUserIds .= "<td><input type='checkbox' id='userList$individualUserId$$' name='userId' value='$individualUserId' $checked><label for='userList$individualUserId$$' title='User'>$userName</label></td>";
	}
	$colCountUser++;
}
my $toBeFilledUser = $colUser - ( $colCountUser % $colUser);
$groupUserIds .= ($toBeFilledUser < $colUser ) ? "<td>&nbsp;</td>" x $toBeFilledUser ."</tr></tbody></table>" : "</tbody></table>";

my $colPermission = 3;
my $colCountPermission=0;
my $groupPermissionIds = "<table id='groupPermissionIds$$' class='display'><thead style='display:none;'><tr>" . "<th></th>" x $colPermission . "</tr></thead><tbody>";
my $allPermissions = $config->getConfigWithType("Permission");
foreach my $individualPermission (@$allPermissions)
{
	my $individualPermissionId = $individualPermission->{"id"};
	my $permissionName = $individualPermission->{"fieldDescription"};
	my $checked = (exists $checkedPermissionId->{$individualPermissionId}) ? "checked='checked'" : "";
	if($colCountPermission % $colPermission == 0)
	{
		$groupPermissionIds .= "<tr><td><input type='checkbox' id='permissionList$individualPermissionId$$' name='permissionId' value='$individualPermissionId' $checked><label for='permissionList$individualPermissionId$$' title='Permission'>$permissionName</label></td>";
	}
	elsif($colCountPermission % $colPermission == $colPermission - 1)
	{
		$groupPermissionIds .= "<td><input type='checkbox' id='permissionList$individualPermissionId$$' name='permissionId' value='$individualPermissionId' $checked><label for='permissionList$individualPermissionId$$' title='Permission'>$permissionName</label></td></tr>";
	}
	else
	{
		$groupPermissionIds .= "<td><input type='checkbox' id='permissionList$individualPermissionId$$' name='permissionId' value='$individualPermissionId' $checked><label for='permissionList$individualPermissionId$$' title='Permission'>$permissionName</label></td>";
	}
	$colCountPermission++;
}
my $toBeFilledPermission = $colPermission - ( $colCountPermission % $colPermission);
$groupPermissionIds .= ($toBeFilledPermission < $colPermission ) ? "<td>&nbsp;</td>" x $toBeFilledPermission ."</tr></tbody></table>" : "</tbody></table>";

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$groupId/$groupId/g;
$html =~ s/\$groupName/$group[2]/g;
$html =~ s/\$groupUserIds/$groupUserIds/g;
$html =~ s/\$groupPermissionIds/$groupPermissionIds/g;
$html =~ s/\$groupDescription/$groupDetails->{'description'}/g;
$html =~ s/\$\$/$$/g;
print header;
print $html;

__DATA__
<form id="editGroup" name="editGroup" action="groupSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="groupId" id="editGroupId" type="hidden" value="$groupId" />
	<table>
	<tr><td style='text-align:right'><label for="editGroupName"><b>Group Name</b></label></td><td><input class='ui-widget-content ui-corner-all' name="name" id="editGroupName" size="40" type="text" maxlength="32" value="$groupName" /></td></tr>
	<tr><td style='text-align:right'><label for="editGroupUser"><b>Users</b></label></td><td>
	$groupUserIds
	</td></tr>
	<tr><td style='text-align:right'><label for="editGroupPermission"><b>Management Permissions</b></label></td><td>
	$groupPermissionIds
	</td></tr>
	<tr><td style='text-align:right'><label for="editGroupDescription"><b>Description</b></label></td><td><textarea class='ui-widget-content ui-corner-all' name="description" id="editGroupDescription" cols="50" rows="10">$groupDescription</textarea></td></tr>
	</table>
</form>
<script>
$( "#groupUserIds$$" ).dataTable({
	"dom": 'lfrtip',
	"scrollY": "300px",
	"scrollCollapse": true,
	"paging": false,
	"searching": false,
	"ordering": false,
	"info": false
});
$( "#groupPermissionIds$$" ).dataTable({
	"dom": 'lfrtip',
	"scrollY": "300px",
	"scrollCollapse": true,
	"paging": false,
	"searching": false,
	"ordering": false,
	"info": false
});
$('#dialog').dialog("option", "title", "Edit Group");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('editGroup'); } }, { text: "Delete", click: function() { deleteItem($groupId); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
</script>