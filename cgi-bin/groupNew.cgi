#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
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

my $colUser = 3;
my $colCountUser=0;
my $groupUserIds = "<table id='groupUserIds$$' class='display'><thead style='display:none;'><tr>" . "<th></th>" x $colUser . "</tr></thead><tbody>";

my $allUsers = $user->getAllFieldsOfUser();
foreach my $individualUser (@$allUsers)
{
	my $individualUserId = $individualUser->{"id"};
	my $userName = $individualUser->{"userName"};
	if($colCountUser % $colUser == 0)
	{
		$groupUserIds .= "<tr><td><input type='checkbox' id='userList$individualUserId$$' name='userId' value='$individualUserId'><label for='userList$individualUserId$$' title='User'>$userName</label></td>";
	}
	elsif($colCountUser % $colUser == $colUser - 1)
	{
		$groupUserIds .= "<td><input type='checkbox' id='userList$individualUserId$$' name='userId' value='$individualUserId'><label for='userList$individualUserId$$' title='User'>$userName</label></td></tr>";
	}
	else
	{
		$groupUserIds .= "<td><input type='checkbox' id='userList$individualUserId$$' name='userId' value='$individualUserId'><label for='userList$individualUserId$$' title='User'>$userName</label></td>";
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
	if($colCountPermission % $colPermission == 0)
	{
		$groupPermissionIds .= "<tr><td><input type='checkbox' id='permissionList$individualPermissionId$$' name='permissionId' value='$individualPermissionId'><label for='permissionList$individualPermissionId$$' title='Permission'>$permissionName</label></td>";
	}
	elsif($colCountPermission % $colPermission == $colPermission - 1)
	{
		$groupPermissionIds .= "<td><input type='checkbox' id='permissionList$individualPermissionId$$' name='permissionId' value='$individualPermissionId'><label for='permissionList$individualPermissionId$$' title='Permission'>$permissionName</label></td></tr>";
	}
	else
	{
		$groupPermissionIds .= "<td><input type='checkbox' id='permissionList$individualPermissionId$$' name='permissionId' value='$individualPermissionId'><label for='permissionList$individualPermissionId$$' title='Permission'>$permissionName</label></td>";
	}
	$colCountPermission++;
}
my $toBeFilledPermission = $colPermission - ( $colCountPermission % $colPermission);
$groupPermissionIds .= ($toBeFilledPermission < $colPermission ) ? "<td>&nbsp;</td>" x $toBeFilledPermission ."</tr></tbody></table>" : "</tbody></table>";

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$groupUserIds/$groupUserIds/g;
$html =~ s/\$groupPermissionIds/$groupPermissionIds/g;
$html =~ s/\$\$/$$/g;

print header;
print $html;

__DATA__
<form id="newGroup" name="newGroup" action="groupSave.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<table>
	<tr><td style='text-align:right'><label for="newGroupName"><b>Group Name</b></label></td><td><input class='ui-widget-content ui-corner-all' name="name" id="newGroupName" size="40" type="text" maxlength="32" /></td></tr>
	<tr><td style='text-align:right'><label for="newGroupUser"><b>Users</b></label></td><td>
	$groupUserIds
	</td></tr>
	<tr><td style='text-align:right'><label for="newGroupPermission"><b>Management Permissions</b></label></td><td>
	$groupPermissionIds
	</td></tr>
	<tr><td style='text-align:right'><label for="newGroupDescription"><b>Description</b></label></td><td><textarea class='ui-widget-content ui-corner-all' name="description" id="newGroupDescription" cols="60" rows="10"></textarea></td></tr>
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
$('#dialog').dialog("option", "title", "New Group");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Save", click: function() { submitForm('newGroup'); } }, { text: "Cancel", click: function() {closeDialog(); } } ] );
</script>