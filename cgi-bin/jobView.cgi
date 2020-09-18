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

my $jobId = param ('jobId');

my $oneJob;
my $jobMeta = '';
my $jobDetail = 'No details available.';
my $circularizedTab = '';
my $circularizedDetails = '';
my $nonVectorTab = '';
my $nonVectorDetails = '';
my $gappedTab = '';
my $gappedDetails = '';
my $bacAssignedTab = '';
my $bacAssignedDetails = '';
my $refreshButton = <<END;
{ text: "Refresh", click: function() {closeDialog();openDialog('jobView.cgi?jobId=$jobId'); } },
END
my $printButton = <<END;
{ text: "Print", click: function() {printDiv('viewJobTabs$$'); } },
END
my $autoRefresh = <<END;
myTimeout = setTimeout(refreshDialog, 5000);
function refreshDialog() {
	closeDialog();
	openDialog('jobView.cgi?jobId=$jobId');
}
END

print header;
if($jobId)
{
	my $job=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'job' AND (name LIKE '$jobId' OR name LIKE '-$jobId')");
	$job->execute();
	if($job->rows > 0)
	{
		my @job = $job->fetchrow_array();
		$oneJob = decode_json $job[8];
		if ($job[7] != -1 )
		{
			$refreshButton = '';
			$autoRefresh = '';
		}
		else
		{
			$printButton = '';
		}
		$jobDetail = "<table>
			<tr><td style='text-align:right'><b>jobId</b></td><td>$jobId</td></tr>
			<tr><td style='text-align:right'><b>name</b></td><td>$oneJob->{'name'}</td></tr>";
		$jobDetail .= ($job[3]) ? "<tr><td style='text-align:right'><b>Polished Contigs</b></td><td><a href='download.cgi?jobId=$jobId&seqType=0' title='Download Polished Sequences' target='hiddenFrame'>$job[3]</a></td></tr>"
			: "<tr><td style='text-align:right'><b>Polished Contigs</b></td><td>$job[3]</td></tr>";
		$jobDetail .= ($job[4]) ? "<tr><td style='text-align:right'><b>Circularized</b></td><td><a href='download.cgi?jobId=$jobId&seqType=1' title='Download Circularized Sequences' target='hiddenFrame'>$job[4]</a></td></tr>"
			: "<tr><td style='text-align:right'><b>Circularized</b></td><td>$job[4]</td></tr>";
		$jobDetail .= ($job[5]) ? "<tr><td style='text-align:right'><b>Non-vector</b></td><td><a href='download.cgi?jobId=$jobId&seqType=3' title='Download Non-vector Sequences' target='hiddenFrame'>$job[5]</a></td></tr>"
			: "<tr><td style='text-align:right'><b>Non-vector</b></td><td>$job[5]</td></tr>";
		$jobDetail .= ($job[6]) ? "<tr><td style='text-align:right'><b>Gapped Contigs</b></td><td><a href='download.cgi?jobId=$jobId&seqType=4' title='Download Gapped Sequences' target='hiddenFrame'>$job[6]</a></td></tr>"
			: "<tr><td style='text-align:right'><b>Gapped Contigs</b></td><td>$job[6]</td></tr>";
		$jobDetail .= ($job[7] > -1) ? "<tr><td style='text-align:right'><b>BAC Id Assigned</b></td><td><a href='download.cgi?jobId=$jobId&seqType=bacId' title='Download BAC Id Assigned Sequences' target='hiddenFrame'>$job[7]</a></td></tr></table>"
			: "<tr><td style='text-align:right'><b>BAC Id Assigned</b></td><td>NA</td></tr></table>";
		my %seqType = (
			0=>'Assembled',
			1=>'Insert',
			2=>'Circularized',
			3=>'NonVector',
			4=>'Gapped',
			5=>'Partial',
			6=>'Vector/Mixer',
			7=>'Mixer',
			8=>'SHORT',
			97=>'Piece',
			98=>'BES',
			99=>'Genome'
			);

		my %bacAssignType = (
			0 => "",
			1 => "TagValid",
			2 => "BesValid",
			3 => "TagValid+BesValid",
			4 => "TagForced"
			);
		
		my $sequencesOfJob=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o > 0 AND o < 5 AND x = ?");
		$sequencesOfJob->execute($jobId);
		while(my @sequencesOfJob = $sequencesOfJob->fetchrow_array())
		{
			$sequencesOfJob[5] = commify ($sequencesOfJob[5]);
			my $sequenceDetails = decode_json $sequencesOfJob[8];
			$sequenceDetails->{'id'} = '' unless (exists $sequenceDetails->{'id'});
			$sequenceDetails->{'description'} = '' unless (exists $sequenceDetails->{'description'});
			my $descriptionLength = length ($sequenceDetails->{'description'});
			if ($descriptionLength > 15)
			{
				$sequenceDetails->{'description'} = "<sup title='$sequenceDetails->{'description'}'>". substr($sequenceDetails->{'description'},0,10). "...". substr($sequenceDetails->{'description'},-3). "</sup>";
			}
			else
			{
				$sequenceDetails->{'description'} = "<sup title='$sequenceDetails->{'description'}'>$sequenceDetails->{'description'}</sup>";
			}
			$sequenceDetails->{'sequence'} = '' unless (exists $sequenceDetails->{'sequence'});
			$sequenceDetails->{'sequence'} =~ tr/a-zA-Z/N/c; #replace nonword characters.
			$sequenceDetails->{'gapList'} = '' unless (exists $sequenceDetails->{'gapList'});
			my $seqTitle = $sequenceDetails->{'id'}.$sequenceDetails->{'description'};
			if($sequencesOfJob[3] == 1 || $sequencesOfJob[3] == 2)
			{
				if($sequencesOfJob[2])
				{
					if($sequencesOfJob[2] =~ /^-/)
					{
						$bacAssignedDetails .= "<tr><td></td><td>Revoked$sequencesOfJob[2] ($bacAssignType{$sequencesOfJob[7]})</td><td>$seqType{$sequencesOfJob[3]}</td><td><a onclick='closeDialog();openDialog(\"seqView.cgi?seqId=$sequencesOfJob[0]\")' title='View this sequence'>$sequencesOfJob[5]</a></td><td><a href='download.cgi?seqId=$sequencesOfJob[0]' title='Download this sequence' target='hiddenFrame'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-disk'></span>$seqTitle</a></td></tr>";
						$circularizedDetails .= "<tr><td></td><td>Revoked$sequencesOfJob[2] ($bacAssignType{$sequencesOfJob[7]})</td><td>$seqType{$sequencesOfJob[3]}</td><td><a onclick='closeDialog();openDialog(\"seqView.cgi?seqId=$sequencesOfJob[0]\")' title='View this sequence'>$sequencesOfJob[5]</a></td><td><a href='download.cgi?seqId=$sequencesOfJob[0]' title='Download this sequence' target='hiddenFrame'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-disk'></span>$seqTitle</a></td></tr>";
					}
					else
					{
						$bacAssignedDetails .= "<tr><td><a onclick='revokeAssignment($sequencesOfJob[0]);' title='Revoke Assigned BAC ID'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-cancel' title='Revoke Assigned BAC ID'></span></a></td><td><a onclick='closeDialog();openDialog(\"cloneView.cgi?cloneName=$sequencesOfJob[2]\")'>$sequencesOfJob[2]</a> ($bacAssignType{$sequencesOfJob[7]})</td><td>$seqType{$sequencesOfJob[3]}</td><td><a onclick='closeDialog();openDialog(\"seqView.cgi?seqId=$sequencesOfJob[0]\")' title='View this sequence'>$sequencesOfJob[5]</a></td><td><a href='download.cgi?seqId=$sequencesOfJob[0]' title='Download this sequence' target='hiddenFrame'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-disk'></span>$seqTitle</a></td></tr>";
						$circularizedDetails .= "<tr><td><a onclick='revokeAssignment($sequencesOfJob[0]);' title='Revoke Assigned BAC ID'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-cancel' title='Revoke Assigned BAC ID'></span></a></td><td><a onclick='closeDialog();openDialog(\"cloneView.cgi?cloneName=$sequencesOfJob[2]\")'>$sequencesOfJob[2]</a> ($bacAssignType{$sequencesOfJob[7]})</td><td>$seqType{$sequencesOfJob[3]}</td><td><a onclick='closeDialog();openDialog(\"seqView.cgi?seqId=$sequencesOfJob[0]\")' title='View this sequence'>$sequencesOfJob[5]</a></td><td><a href='download.cgi?seqId=$sequencesOfJob[0]' title='Download this sequence' target='hiddenFrame'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-disk'></span>$seqTitle</a></td></tr>";
					}
				}
				else
				{
					$circularizedDetails .= "<tr><td></td><td>NA</td><td>$seqType{$sequencesOfJob[3]}</td><td><a onclick='closeDialog();openDialog(\"seqView.cgi?seqId=$sequencesOfJob[0]\")' title='View this sequence'>$sequencesOfJob[5]</a></td><td><a href='download.cgi?seqId=$sequencesOfJob[0]' title='Download this sequence' target='hiddenFrame'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-disk'></span>$seqTitle</a></td></tr>";
				}
			}
			elsif($sequencesOfJob[3] == 3)
			{
				if($sequencesOfJob[2])
				{
					if($sequencesOfJob[2] =~ /^-/)
					{
						$bacAssignedDetails .= "<tr><td></td><td>Revoked$sequencesOfJob[2] ($bacAssignType{$sequencesOfJob[7]})</td><td>$seqType{$sequencesOfJob[3]}</td><td><a onclick='closeDialog();openDialog(\"seqView.cgi?seqId=$sequencesOfJob[0]\")' title='View this sequence'>$sequencesOfJob[5]</a></td><td><a href='download.cgi?seqId=$sequencesOfJob[0]' title='Download this sequence' target='hiddenFrame'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-disk'></span>$seqTitle</a></td></tr>";
						$nonVectorDetails .= "<tr><td></td><td>Revoked$sequencesOfJob[2] ($bacAssignType{$sequencesOfJob[7]})</td><td>$seqType{$sequencesOfJob[3]}</td><td><a onclick='closeDialog();openDialog(\"seqView.cgi?seqId=$sequencesOfJob[0]\")' title='View this sequence'>$sequencesOfJob[5]</a></td><td><a href='download.cgi?seqId=$sequencesOfJob[0]' title='Download this sequence' target='hiddenFrame'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-disk'></span>$seqTitle</a></td></tr>";
					}
					else
					{
						$bacAssignedDetails .= "<tr><td><a onclick='revokeAssignment($sequencesOfJob[0]);' title='Revoke Assigned BAC ID'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-cancel' title='Revoke Assigned BAC ID'></span></a></td><td><a onclick='closeDialog();openDialog(\"cloneView.cgi?cloneName=$sequencesOfJob[2]\")'>$sequencesOfJob[2]</a> ($bacAssignType{$sequencesOfJob[7]})</td><td>$seqType{$sequencesOfJob[3]}</td><td><a onclick='closeDialog();openDialog(\"seqView.cgi?seqId=$sequencesOfJob[0]\")' title='View this sequence'>$sequencesOfJob[5]</a></td><td><a href='download.cgi?seqId=$sequencesOfJob[0]' title='Download this sequence' target='hiddenFrame'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-disk'></span>$seqTitle</a></td></tr>";
						$nonVectorDetails .= "<tr><td><a onclick='revokeAssignment($sequencesOfJob[0]);' title='Revoke Assigned BAC ID'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-cancel' title='Revoke Assigned BAC ID'></span></a></td><td><a onclick='closeDialog();openDialog(\"cloneView.cgi?cloneName=$sequencesOfJob[2]\")'>$sequencesOfJob[2]</a> ($bacAssignType{$sequencesOfJob[7]})</td><td>$seqType{$sequencesOfJob[3]}</td><td><a onclick='closeDialog();openDialog(\"seqView.cgi?seqId=$sequencesOfJob[0]\")' title='View this sequence'>$sequencesOfJob[5]</a></td><td><a href='download.cgi?seqId=$sequencesOfJob[0]' title='Download this sequence' target='hiddenFrame'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-disk'></span>$seqTitle</a></td></tr>";
					}
				}
				else
				{
					$nonVectorDetails .= "<tr><td></td><td>NA</td><td>$seqType{$sequencesOfJob[3]}</td><td><a onclick='closeDialog();openDialog(\"seqView.cgi?seqId=$sequencesOfJob[0]\")' title='View this sequence'>$sequencesOfJob[5]</a></td><td><a href='download.cgi?seqId=$sequencesOfJob[0]' title='Download this sequence' target='hiddenFrame'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-disk'></span>$seqTitle</a></td></tr>";
				}
			}
			elsif($sequencesOfJob[3] == 4)
			{
				if($sequencesOfJob[2])
				{
					if($sequencesOfJob[2] =~ /^-/)
					{
						$bacAssignedDetails .= "<tr><td></td><td>Revoked$sequencesOfJob[2] ($bacAssignType{$sequencesOfJob[7]})</td><td>$seqType{$sequencesOfJob[3]}</td><td><a onclick='closeDialog();openDialog(\"seqView.cgi?seqId=$sequencesOfJob[0]\")' title='View this sequence'>$sequencesOfJob[5]</a></td><td><a href='download.cgi?seqId=$sequencesOfJob[0]' title='Download this sequence' target='hiddenFrame'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-disk'></span>$seqTitle</a></td></tr>";
						$gappedDetails .= "<tr><td></td><td>Revoked$sequencesOfJob[2] ($bacAssignType{$sequencesOfJob[7]})</td><td>$seqType{$sequencesOfJob[3]}</td><td><a onclick='closeDialog();openDialog(\"seqView.cgi?seqId=$sequencesOfJob[0]\")' title='View this sequence'>$sequencesOfJob[5]</a></td><td><a href='download.cgi?seqId=$sequencesOfJob[0]' title='Download this sequence' target='hiddenFrame'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-disk'></span>$seqTitle</a></td></tr>";
					}
					else
					{
						$bacAssignedDetails .= "<tr><td><a onclick='revokeAssignment($sequencesOfJob[0]);' title='Revoke Assigned BAC ID'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-cancel' title='Revoke Assigned BAC ID'></span></a></td><td><a onclick='closeDialog();openDialog(\"cloneView.cgi?cloneName=$sequencesOfJob[2]\")'>$sequencesOfJob[2]</a> ($bacAssignType{$sequencesOfJob[7]})</td><td>$seqType{$sequencesOfJob[3]}</td><td><a onclick='closeDialog();openDialog(\"seqView.cgi?seqId=$sequencesOfJob[0]\")' title='View this sequence'>$sequencesOfJob[5]</a></td><td><a href='download.cgi?seqId=$sequencesOfJob[0]' title='Download this sequence' target='hiddenFrame'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-disk'></span>$seqTitle</a></td></tr>";
						$gappedDetails .= "<tr><td><a onclick='revokeAssignment($sequencesOfJob[0]);' title='Revoke Assigned BAC ID'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-cancel' title='Revoke Assigned BAC ID'></span></a></td><td><a onclick='closeDialog();openDialog(\"cloneView.cgi?cloneName=$sequencesOfJob[2]\")'>$sequencesOfJob[2]</a> ($bacAssignType{$sequencesOfJob[7]})</td><td>$seqType{$sequencesOfJob[3]}</td><td><a onclick='closeDialog();openDialog(\"seqView.cgi?seqId=$sequencesOfJob[0]\")' title='View this sequence'>$sequencesOfJob[5]</a></td><td><a href='download.cgi?seqId=$sequencesOfJob[0]' title='Download this sequence' target='hiddenFrame'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-disk'></span>$seqTitle</a></td></tr>";
					}
				}
				else
				{
					$gappedDetails .= "<tr><td></td><td>NA</td><td>$seqType{$sequencesOfJob[3]}</td><td><a onclick='closeDialog();openDialog(\"seqView.cgi?seqId=$sequencesOfJob[0]\")'><a onclick='closeDialog();openDialog(\"seqView.cgi?seqId=$sequencesOfJob[0]\")' title='View this sequence'>$sequencesOfJob[5]</a></a></td><td><a href='download.cgi?seqId=$sequencesOfJob[0]' title='Download this sequence' target='hiddenFrame'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-disk'></span>$seqTitle</a></td></tr>";
				}
			}
		}

		if($job[4])
		{
			$circularizedTab = "<li><a href='#viewJobTabs-Circularized'>Circularized ($job[4])</a></li>";
			$circularizedDetails = "<div id='viewJobTabs-Circularized'>
				<a style='float: right;' href='download.cgi?jobId=$jobId&seqType=1' title='Download Circularized Sequences' target='hiddenFrame'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-disk' title='Download'></span>Download Circularized Sequences</a>
				<table class='jobDetails display'><thead><tr><th></th><th>BAC ID</th><th>seqType</th><th>Length (bp)</th><th>Seq Details</th></tr></thead><tbody>$circularizedDetails</tbody></table>
				</div>";
		}
		if($job[5])
		{
			$nonVectorTab = "<li><a href='#viewJobTabs-nonVector'>Non-vector ($job[5])</a></li>";
			$nonVectorDetails = "<div id='viewJobTabs-nonVector'>
				<a style='float: right;' href='download.cgi?jobId=$jobId&seqType=3' title='Download Non-vector Sequences' target='hiddenFrame'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-disk' title='Download'></span>Download Non-vector Sequences</a>
				<table class='jobDetails display'><thead><tr><th></th><th>BAC ID</th><th>seqType</th><th>Length (bp)</th><th>Seq Details</th></tr></thead><tbody>$nonVectorDetails</tbody></table>
				</div>";
		}
		if($job[6])
		{
			$gappedTab = "<li><a href='#viewJobTabs-Gapped'>Gapped ($job[6])</a></li>";
			$gappedDetails = "<div id='viewJobTabs-Gapped'>
				<a style='float: right;' href='download.cgi?jobId=$jobId&seqType=4' title='Download Gapped Sequences' target='hiddenFrame'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-disk' title='Download'></span>Download Gapped Sequences</a>
				<table class='jobDetails display'><thead><tr><th></th><th>BAC ID</th><th>seqType</th><th>Length (bp)</th><th>Seq Details</th></tr></thead><tbody>$gappedDetails</tbody></table>
				</div>";
		}
		if($job[7] > -1)
		{
			$bacAssignedTab = "<li><a href='#viewJobTabs-bacAssigned'>BAC ID Assigned ($job[7])</a></li>";
			$bacAssignedDetails = "<div id='viewJobTabs-bacAssigned'>
				<a style='float: right;' href='download.cgi?jobId=$jobId&seqType=bacId' title='Download BAC Id Assigned Sequences' target='hiddenFrame'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-disk' title='Download'></span>Download BAC ID Assigned Sequences</a>
				<table class='jobDetails display'><thead><tr><th></th><th>BAC ID</th><th>seqType</th><th>Length (bp)</th><th>Seq Details</th></tr></thead><tbody>$bacAssignedDetails</tbody></table>
				</div>";
		}
	}
	else
	{
		print <<END;
		<script>
			errorPop("No jobId found!");
			closeDialog();
		</script>	
END
		exit;
	}

	$jobMeta = <<END;
	<table>
		<tr><td style='text-align:right'><b>jobId</b></td><td>$jobId</td></tr>
		<tr><td style='text-align:right'><b>name</b></td><td>$oneJob->{'name'}</td></tr>
		<tr><td style='text-align:right'><b>sampleName</b></td><td>$oneJob->{'sampleName'}</td></tr>
		<tr><td style='text-align:right'><b>protocolName</b></td><td>$oneJob->{'protocolName'}</td></tr>
		<tr><td style='text-align:right'><b>version</b></td><td>$oneJob->{'version'}</td></tr>
		<tr><td style='text-align:right'><b>collectionProtocol</b></td><td>$oneJob->{'collectionProtocol'}</td></tr>
		<tr><td style='text-align:right'><b>plateId</b></td><td>$oneJob->{'plateId'}</td></tr>
		<tr><td style='text-align:right'><b>primaryProtocol</b></td><td>$oneJob->{'primaryProtocol'}</td></tr>
		<tr><td style='text-align:right'><b>comments</b></td><td>$oneJob->{'comments'}</td></tr>
		<tr><td style='text-align:right'><b>createdBy</b></td><td>$oneJob->{'createdBy'}</td></tr>
		<tr><td style='text-align:right'><b>modifiedBy</b></td><td>$oneJob->{'modifiedBy'}</td></tr>
		<tr><td style='text-align:right'><b>lastHeartbeat</b></td><td>$oneJob->{'lastHeartbeat'}</td></tr>
		<tr><td style='text-align:right'><b>whenCreated</b></td><td>$oneJob->{'whenCreated'}</td></tr>
		<tr><td style='text-align:right'><b>whenStarted</b></td><td>$oneJob->{'whenStarted'}</td></tr>
		<tr><td style='text-align:right'><b>whenEnded</b></td><td>$oneJob->{'whenEnded'}</td></tr>
		<tr><td style='text-align:right'><b>whenModified</b></td><td>$oneJob->{'whenModified'}</td></tr>
	</table>
END

	undef $/;# enable slurp mode
	my $html = <DATA>;
	$html =~ s/\$\$/$$/g;
	$html =~ s/\$jobId/$jobId/g;
	$html =~ s/\$jobDetail/$jobDetail/g;
	$html =~ s/\$circularizedTab/$circularizedTab/g;
	$html =~ s/\$circularizedDetails/$circularizedDetails/g;
	$html =~ s/\$nonVectorTab/$nonVectorTab/g;
	$html =~ s/\$nonVectorDetails/$nonVectorDetails/g;
	$html =~ s/\$gappedTab/$gappedTab/g;
	$html =~ s/\$gappedDetails/$gappedDetails/g;
	$html =~ s/\$bacAssignedTab/$bacAssignedTab/g;
	$html =~ s/\$bacAssignedDetails/$bacAssignedDetails/g;
	$html =~ s/\$jobMeta/$jobMeta/g;
	$html =~ s/\$refreshButton/$refreshButton/g;
	$html =~ s/\$printButton/$printButton/g;
	$html =~ s/\$autoRefresh/$autoRefresh/g;
	print $html;
}
else
{
	print <<END;
	<script>
		errorPop("No jobId provided!");
		closeDialog();
	</script>	
END
	exit;
}


__DATA__
<div id="viewJobTabs$$" name="viewJobTabs$$">
	<ul>
	<li><a href="#viewJobTabs-Job">Summary</a></li>
	$circularizedTab
	$nonVectorTab
	$gappedTab
	$bacAssignedTab
	<li><a href="#viewJobTabs-Meta">Meta</a></li>
	</ul>
	<div id="viewJobTabs-Job">
	$jobDetail
	</div>
	$circularizedDetails
	$nonVectorDetails
	$gappedDetails
	$bacAssignedDetails
	<div id="viewJobTabs-Meta">
	$jobMeta
	</div>
</div>
<script>
$autoRefresh
buttonInit();
$( "#viewJobTabs$$" ).tabs();
$( ".jobDetails" ).dataTable({
	"paging": false,
	"searching": false,
	"info": false,
	"columnDefs": [
    { "orderable": false, "targets": 0 }
  ]
});
$('#dialog').dialog("option", "title", "View Assembly Job $jobId");
$( "#dialog" ).dialog( "option", "buttons", [ $refreshButton{ text: "reRun Post-HGAP", click: function() {closeDialog();openDialog('jobPost.cgi?jobId=$jobId'); } }, $printButton { text: "OK", click: function() {closeDialog(); } } ] );
</script>

