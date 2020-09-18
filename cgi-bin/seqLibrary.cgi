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

my $commoncfg = readConfig("main.conf");
my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});

my $libraryId = param ('libraryId') || '';
my $active = 0;
my $activeDetector = 0;
my $cookiePlate = cookie('plate') || '';
my $config = new config;
my $userConfig = new userConfig;

my $button = "<ul id='plateInfoMenu$libraryId$$' style='margin-top: .3em; width: 200px;'>
					<li><a><b>Plates in Library</b></a>
						<ul style='z-index: 1000;white-space: nowrap;'>
							<li><a onclick='refresh(\"viewLibraryTabs$libraryId\")'>Refresh</a></li>
						</ul>
					</li>
				</ul>";
my $pools = '';
my $plates = '';
if ($libraryId)
{
#	my $configPoolsPerPage = $config->getFieldDefaultWithFieldName("poolsPerPage");
#	my $poolsPerPage = $userConfig->getFieldValueWithUserIdAndFieldName($userId,"poolsPerPage");
	my $poolsPerPage = 100;
#	my ($inputType, $lengthMenu, $inputDefault) = split(/:/, $configPoolsPerPage);

	my $library=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$library->execute($libraryId);
	my @library = $library->fetchrow_array();
	my $plateId;
	my $plateInLibrary=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'plate' AND z = ? ORDER BY o");
	$plateInLibrary->execute($libraryId);
	while (my @plateInLibrary = $plateInLibrary->fetchrow_array())
	{
		$plateId->{$plateInLibrary[2]} = 1;
	}
	for (sort keys %{$plateId})
	{
		$plates = "
		<div id='platesInLibrary$libraryId$$'><ul>\n" unless($plates);
		$plates .= "<li><a href='plateClone.cgi?plate=$_&libraryId=$libraryId'>$_</a></li>\n";
		$active = $activeDetector if ($cookiePlate eq $libraryId.$_);
		$activeDetector++;
	}
	$plates .= "</ul></div>\n" if ($plates);
	unless($plates)
	{
		$plates .= "<p class='ui-state-error ui-corner-all' style='padding: .7em;'><span class='ui-icon ui-icon-alert' style='float: left; margin-right: .3em;'></span>
			<strong>No plates in this seqLibrary.</strong></p>";
	}
	my $pool=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'pool' AND x = ? ORDER BY name");
	$pool->execute($libraryId);
	while (my @pool = $pool->fetchrow_array())
	{
		$pools = "
		<ul id='poolInfoMenu$libraryId$$' style='margin-top: .3em; width: 200px;'>
			<li><a><b>Sequencing Pools</b></a>
				<ul style='z-index: 1000;white-space: nowrap;'>
					<li><a onclick='printDiv(\"poolList$libraryId$$\")'>Print</a></li>
					<li><a href='download.cgi?seqLibraryId=$libraryId' target='hiddenFrame'>Download</a></li>
					<li><a onclick='refresh(\"viewLibraryTabs$libraryId\")'>Refresh</a></li>
				</ul>
			</li>
		</ul>
		<div id='poolList$libraryId$$' name='poolList$libraryId$$'>
		<table id='poolsInLibrary$libraryId$$' class='display' cellspacing='0' width='100%'>
        <thead>
            <tr>
                <th><b>Pool Name</b></th>
                <th><b>Pool Size</b></th>
                <th><b>BACs</b></th>
                <th><b>Notes</b></th>
                <th><b>Assembly</b></th>
            </tr>
        </thead>
        <tbody>"
		unless($pools);

		my $jobIds;
		my $jobId=$dbh->prepare("SELECT * FROM link WHERE parent = ? AND type LIKE 'poolJob' ORDER BY child DESC");
		$jobId->execute($pool[0]);
		while(my @jobId = $jobId->fetchrow_array())
		{
			$jobIds .= "<a onclick='openDialog(\"jobView.cgi?jobId=$jobId[1]\")'>$jobId[1]</a> ";
		}
		$jobIds = "<a onclick='openDialog(\"poolEdit.cgi?poolId=$pool[0]\")'>link to a new job</a>" unless ($jobIds);

		my $clones;
		my $sequencedNumber = 0;
		my $notSequencedNumber = 0;
		my $clone=$dbh->prepare("SELECT clones.* FROM link,clones WHERE link.type LIKE 'poolClone' AND link.child = clones.name AND link.parent = ? ORDER BY clones.name");
		$clone->execute($pool[0]);
		while(my @clone = $clone->fetchrow_array())
		{
			$clones .= ($clone[6]) ? 
			"<div style='float: left; margin-right: .7em;' title='$clone[5]'><span class='ui-icon ui-icon-check' style='float: left; margin-right: 0;'></span><a onclick='openDialog(\"cloneView.cgi?cloneName=$clone[1]\")'>$clone[1]</a></div>" : 
			"<div class='ui-state-disabled' style='float: left; margin-right: .7em;' title='$clone[5]'><span class='ui-icon ui-icon-cancel' style='float: left; margin-right: 0;'></span><a onclick='openDialog(\"cloneView.cgi?cloneName=$clone[1]\")'>$clone[1]</a></div>";
			($clone[6]) ? $sequencedNumber++ : $notSequencedNumber++;
		}
		$clones = "<a onclick='openDialog(\"poolEdit.cgi?poolId=$pool[0]\")'>Add clones</a>" unless ($clones);
		my $totalCloneNumber = $sequencedNumber + $notSequencedNumber;
		$pools .= ($notSequencedNumber > 0) ? "<tr><td><div style='position: relative;'><a id='pool$pool[0]$$' onmouseover='editIconShow(\"pool$pool[0]$$\")' onmouseout='editIconHide(\"pool$pool[0]$$\")' onclick='openDialog(\"poolEdit.cgi?poolId=$pool[0]\")' title='$pool[8]'>$pool[2]</a></div></td>
		<td>$totalCloneNumber</td>
		<td>$clones</td>
		<td><sup>Sequenced: $sequencedNumber; Not sequenced: $notSequencedNumber.</sup><br>$pool[8]</td>
		<td><div id='jobId$pool[0]'>$jobIds</div></td></tr>" : "<tr><td><div style='position: relative;'><a id='pool$pool[0]$$' onmouseover='editIconShow(\"pool$pool[0]$$\")' onmouseout='editIconHide(\"pool$pool[0]$$\")' onclick='openDialog(\"poolEdit.cgi?poolId=$pool[0]\")' title='$pool[8]'>$pool[2]</a></div></td>
		<td>$totalCloneNumber</td>
		<td>$clones</td>
		<td><sup>Sequenced: $sequencedNumber.</sup><br>$pool[8]</td>
		<td><div id='jobId$pool[0]'>$jobIds</div></td></tr>";
	}
	$pools .= "</tbody></table></div>\n" if ($pools);
	
	undef $/;# enable slurp mode
	my $html = <DATA>;
	$html =~ s/\$button/$button/g;
	$html =~ s/\$plates/$plates/g;
	$html =~ s/\$poolsPerPage/$poolsPerPage/g;
	$html =~ s/\$pools/$pools/g;
	$html =~ s/\$libraryId/$libraryId/g;
	$html =~ s/\$\$/$$/g;
	$html =~ s/\$active/$active/g;
	$html =~ s/\$commoncfg->{HTDOCS}/$commoncfg->{HTDOCS}/g;

	print header;
	print $html;
}
else
{
	print header(-type=>'text/html',-status=>'402 Invalid operation');
	exit;
}


__DATA__
$button
$plates
$pools
<script>
buttonInit();
loadingHide();
$( "#platesInLibrary$libraryId$$" ).tabs({
	// loading spinner
	beforeLoad: function(event, ui) {
		ui.panel.html('<img src="$commoncfg->{HTDOCS}/css/images/loading.gif" style="vertical-align:middle;"> Loading...');
	},
	create: function (e, ui) {
		loadingShow();
	},
	activate: function (e, ui) {
		loadingShow();
	},
	active: $active
});
$( "#plateInfoMenu$libraryId$$" ).menu();
$( "#poolInfoMenu$libraryId$$" ).menu();
$( "#poolsInLibrary$libraryId$$" ).dataTable({
	"lengthMenu": [ 5, 10, 25, 50, 100 ],
	"pageLength": $poolsPerPage,
});
</script>