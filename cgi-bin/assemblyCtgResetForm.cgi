#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use DBI;
use lib "lib/";
use lib "lib/pangu";
use pangu;
use user;
use userCookie;

my $userCookie = new userCookie;
my $userId = (cookie('cid')) ? $userCookie->checkCookie(cookie('cid')) : 0;
exit if (!$userId);

my $user = new user;
my $userDetail = $user->getAllFieldsWithUserId($userId);
my $userName = $userDetail->{"userName"};

my $commoncfg = readConfig("main.conf");
my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});

print header;
my $assemblyId = param ('assemblyId') || '';


if($assemblyId)
{
	my $assembly=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$assembly->execute($assemblyId);
	my @assembly = $assembly->fetchrow_array();

	unless ($assembly[7] == 1 || $assembly[7] == 0) # exit if for frozen or running assembly
	{
		print <<END;
<script>
	closeDialog();
	errorPop("This assembly is running or frozen.");
</script>	
END
		exit;
	}
	undef $/;# enable slurp mode
	my $html = <DATA>;
    $html =~ s/\$\$/$$/g;
	$html =~ s/\$assemblyId/$assemblyId/g;
	$html =~ s/\$commoncfg->{HTDOCS}/$commoncfg->{HTDOCS}/g;

	print $html;
}

__DATA__
<form id="assemblyCtgReset" name="assemblyCtgReset" action="assemblyCtgReset.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
<input name="assemblyId" id="assemblyId$$" type="hidden" value="$assemblyId" />
Click "Reset" to clear all redundancy settings.
</form>
<script>
$('#dialog').dialog("option", "title", "Reset AssemblyCtg Redundancy");
$( "#dialog" ).dialog( "option", "buttons", [{ text: "Reset", click: function() { submitForm('assemblyCtgReset'); } },{ text: "OK", click: function() { closeDialog(); } } ] );
</script>