#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use JSON::XS;
use DBI;
use lib "lib/";
use lib "lib/pangu";
use pangu;
use config;
use userCookie;

my $userCookie = new userCookie;
my $userId = (cookie('cid')) ? $userCookie->checkCookie(cookie('cid')) : 0;
exit if (!$userId);

my $commoncfg = readConfig("main.conf");
my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});
my $config = new config;
my $JOBURL = $config->getFieldValueWithFieldName("JOBURL");
my $JOBREALTIME = $config->getFieldValueWithFieldName("JOBREALTIME");

my $job;
my $hiddenJob;
my $assignedSequence;
my $jobList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'job'");
$jobList->execute();
while (my @jobList = $jobList->fetchrow_array())
{
	if($jobList[2] =~ /^-/)
	{
		$jobList[2] =~ s/-//g;
		$hiddenJob->{$jobList[2]} = 1;
	}
	@{$job->{$jobList[2]}} = @jobList;
}

my $poolClones;
my $poolCloneNumber;
my $poolClone=$dbh->prepare("SELECT parent,child FROM link WHERE type LIKE 'poolClone'");
$poolClone->execute();
while(my @poolClone = $poolClone->fetchrow_array())
{
	$poolClones->{$poolClone[0]} .= "$poolClone[1] ";
	$poolCloneNumber->{$poolClone[0]}++;
}

my $jobToPool;
my $jobToPoolId;
my $jobToPoolClones;

my $poolJob=$dbh->prepare("SELECT link.child,matrix.* FROM link,matrix WHERE link.type LIKE 'poolJob' AND link.parent = matrix.id ORDER BY matrix.o");
$poolJob->execute();
while (my @poolJob = $poolJob->fetchrow_array())
{
	$jobToPoolClones->{$poolJob[0]} = (exists $poolClones->{$poolJob[1]}) ? "<a onclick='openDialog(\"poolView.cgi?poolId=$poolJob[1]\")' title='$poolClones->{$poolJob[1]}'>$poolCloneNumber->{$poolJob[1]}</a>" : "";
	$jobToPoolId->{$poolJob[0]} = $poolJob[1];
	my $poolToLibrary = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$poolToLibrary->execute($poolJob[5]);
	my @poolToLibrary = $poolToLibrary->fetchrow_array();
	$jobToPool->{$poolJob[0]} = "<a onclick='openDialog(\"poolView.cgi?poolId=$poolJob[1]\")'>$poolToLibrary[2]$poolJob[3]</a>";
}

my $jobs = "<div id='jobListDiv$$' name='jobListDiv$$'>
	<form id='jobList$$' name='jobList$$'>
	<table id='jobs$$' class='display'>
		<thead>
			<tr>
				<th>
					<input type='checkbox' id='checkAllBox$$' name='checkAllBox$$' value='Check all' checked='checked' onClick='checkAll(\"jobId\");return false;' title='Check all'>
					<input type='checkbox' id='uncheckAllBox$$' name='uncheckAllBox$$' value='Uncheck all' onClick='uncheckAll(\"jobId\");return false;' title='Uncheck all'>
				</th>
				<th><b>Job Id</b></th>
				<th><b>Pool Name</b></th>
				<th><b>Pooled Clones</b></th>
				<th><b>Polished Contigs</b></th>
				<th><b>Circularized</b>(%)</th>
				<th><b>Non-vector</b></th>
				<th><b>Gapped Contigs</b></th>
				<th title='Including Non-vector and Gapped Unitigs'><b>BAC Id Assigned</b>(%)</th>
				<th><b>Job Name</b></th>
			</tr>
		</thead>
		<tbody>";
if($JOBREALTIME)
{
	my $allJobsOnPac = `curl -d 'options={"page":1,"rows":0,"sortOrder":"desc","sortBy":"jobId"}' $JOBURL`;
	my $allJobs  = decode_json $allJobsOnPac;
	for my $eachJob (@{$allJobs->{'rows'}})
	{
		next if ($eachJob->{'jobStatus'}  ne 'Completed');
		my $jobDetails = encode_json $eachJob;
		my $jobId = $eachJob->{'jobId'};
		next if (exists $hiddenJob->{$jobId});
		$jobToPool->{$jobId} = "NA" unless (exists $jobToPool->{$jobId});
		$jobToPoolClones->{$jobId} = "NA" unless (exists $jobToPoolClones->{$jobId});
		my $circularizedPercent = "NA";
		my $bacAssignedPercent = "NA";
		if (exists $jobToPoolId->{$jobId})
		{
			$circularizedPercent = int ($job->{$jobId}[4] * 100 / $poolCloneNumber->{$jobToPoolId->{$jobId}}) if (exists $poolCloneNumber->{$jobToPoolId->{$jobId}});
			$bacAssignedPercent = int ($job->{$jobId}[7] * 100 / $poolCloneNumber->{$jobToPoolId->{$jobId}}) if (exists $poolCloneNumber->{$jobToPoolId->{$jobId}});
		}
		$jobs .= (exists $job->{$jobId}) ? "<tr>
			<td style='text-align:center;'><input type='checkbox' id='jobList$jobId$$' name='jobId' value='$jobId'></td>
			<td title=''><a onclick='openDialog(\"jobView.cgi?jobId=$jobId\")'>$jobId</a></td>
			<td title='Pool'>$jobToPool->{$jobId}</a></td>
			<td title='Clones'>$jobToPoolClones->{$jobId}</td>
			<td title='Polished Contigs'>$job->{$jobId}[3]</td>
			<td title='Circularized Contigs'>$job->{$jobId}[4] ($circularizedPercent)</td>
			<td title='Non-vector'>$job->{$jobId}[5]</td>
			<td title='Gapped Contigs'>$job->{$jobId}[6]</td>
			<td title='Assigned Contigs Including Non-vector and Gapped Unitigs'>$job->{$jobId}[7] ($bacAssignedPercent)</td>
			<td title='Sample: $eachJob->{'sampleName'}\nProtocol: $eachJob->{'protocolName'}\n$eachJob->{'whenCreated'}'><a href='$JOBURL/$jobId/contents/data/polished_assembly.fasta.gz' target='hiddenFrame'>$eachJob->{'name'}</a></td>
			</tr>\n" : "<tr>
			<td style='text-align:center;'><input type='checkbox' id='jobList$jobId$$' name='jobId' value='$jobId'></td>
			<td title=''><a onclick='openDialog(\"jobView.cgi?jobId=$jobId\")'>$jobId</a></td>
			<td title='Pool'>$jobToPool->{$jobId}</a></td>
			<td title='Clones'>$jobToPoolClones->{$jobId}</td>
			<td title='Polished Contigs'></td>
			<td title='Circularized Contigs'></td>
			<td title='Non-vector'></td>
			<td title='Gapped Contigs'></td>
			<td title='Assigned Contigs Including Non-vector and Gapped Unitigs'></td>
			<td title='Sample: $eachJob->{'sampleName'}\nProtocol: $eachJob->{'protocolName'}\n$eachJob->{'whenCreated'}'><a href='$JOBURL/$jobId/contents/data/polished_assembly.fasta.gz' target='hiddenFrame'>$eachJob->{'name'}</a></td>
			</tr>\n";
		unless (exists $job->{$jobId})
		{
			my $insertJob=$dbh->prepare("INSERT INTO matrix VALUES ('', 'job', ?, 0, 0, 0, 0, -2, ?, 'AutoUser', NOW())");
			$insertJob->execute($jobId,$jobDetails);
		}
	}
}
else
{
	foreach my $jobId (sort {$b <=> $a} keys %$job)
	{
		next if (exists $hiddenJob->{$jobId});
		$jobToPool->{$jobId} = "NA" unless (exists $jobToPool->{$jobId});
		$jobToPoolClones->{$jobId} = "NA" unless (exists $jobToPoolClones->{$jobId});
		my $circularizedPercent = "NA";
		if (exists $poolCloneNumber->{$jobToPoolId->{$jobId}})
		{
			$circularizedPercent = int ($job->{$jobId}[4] * 100 / $poolCloneNumber->{$jobToPoolId->{$jobId}});
		}
		my $bacAssignedPercent = "NA";
		if (exists $poolCloneNumber->{$jobToPoolId->{$jobId}})
		{
			$bacAssignedPercent = int ($job->{$jobId}[7] * 100 / $poolCloneNumber->{$jobToPoolId->{$jobId}});
		}
		my $jobDetails = decode_json $job->{$jobId}[8];
		$jobs .= (-e "$commoncfg->{POLISHED}/$jobId.fasta.gz") ? "<tr>
			<td style='text-align:center;'><input type='checkbox' id='jobList$jobId$$' name='jobId' value='$jobId'></td>
			<td title=''><a onclick='openDialog(\"jobView.cgi?jobId=$jobId\")'>$jobId</a></td>
			<td title='Pool'>$jobToPool->{$jobId}</a></td>
			<td title='Clones'>$jobToPoolClones->{$jobId}</td>
			<td title='Polished Contigs'>$job->{$jobId}[3]</td>
			<td title='Circularized Contigs'>$job->{$jobId}[4] ($circularizedPercent)</td>
			<td title='Non-vector'>$job->{$jobId}[5]</td>
			<td title='Gapped Contigs'>$job->{$jobId}[6]</td>
			<td title='Assigned Contigs'>$job->{$jobId}[7] ($bacAssignedPercent)</td>
			<td title='Sample: $jobDetails->{'sampleName'}\nProtocol: $jobDetails->{'protocolName'}\n$jobDetails->{'whenCreated'}'><a href='$commoncfg->{POLISHEDURL}/$jobId.fasta.gz' target='hiddenFrame'>$jobDetails->{'name'}</a></td>
			</tr>\n" : "<tr>
			<td style='text-align:center;'><input type='checkbox' id='jobList$jobId$$' name='jobId' value='$jobId'></td>
			<td title=''><a onclick='openDialog(\"jobView.cgi?jobId=$jobId\")'>$jobId</a></td>
			<td title='Pool'>$jobToPool->{$jobId}</a></td>
			<td title='Clones'>$jobToPoolClones->{$jobId}</td>
			<td title='Polished Contigs'>$job->{$jobId}[3]</td>
			<td title='Circularized Contigs'>$job->{$jobId}[4] ($circularizedPercent)</td>
			<td title='Non-vector'>$job->{$jobId}[5]</td>
			<td title='Gapped Contigs'>$job->{$jobId}[6]</td>
			<td title='Assigned Contigs Including Non-vector and Gapped Unitigs'>$job->{$jobId}[7] ($bacAssignedPercent)</td>
			<td title='Sample: $jobDetails->{'sampleName'}\nProtocol: $jobDetails->{'protocolName'}\n$jobDetails->{'whenCreated'}'>$jobDetails->{'name'}</td>
			</tr>\n";
	}
}
$jobs .= "</tbody></table></form></div>\n" if($jobs);
my $button = "<ul id='jobInfoMenu$$' style='margin-top: .3em; width: 250px;'>
		<li><a><span class='ui-icon ui-icon-triangle-1-e'></span><b>PacBio HGAP Jobs</b></a>
			<ul style='z-index: 1000;white-space: nowrap;'>
				<li><a><span class='ui-icon ui-icon-bullet'></span>Job List</a>
					<ul style='z-index: 1000;white-space: nowrap;'>
						<li><a onclick='printDiv(\"jobListDiv$$\")'><span class='ui-icon ui-icon-print'></span>Print</a></li>
						<li><a href='jobExport.cgi' target='hiddenFrame'><span class='ui-icon ui-icon-disk'></span>Export</a></li>
						<li><a onclick='openDialog(\"jobImportForm.cgi\")'><span class='ui-icon ui-icon-arrowthickstop-1-n'></span>Import From File</a></li>
					</ul>
				</li>
				<li><a onclick='openDialogForm(\"jobHideForm.cgi\",\"jobList$$\")'><span class='ui-icon ui-icon-trash'></span>Hide Jobs</a></li>
				<li><a onclick='openDialog(\"jobShowForm.cgi\")'><span class='ui-icon ui-icon-lightbulb'></span>Manage Hidden Jobs</a></li>
				<li><a onclick='openDialogForm(\"jobPostForm.cgi\",\"jobList$$\")'><span class='ui-icon ui-icon-play'></span>Run postHGAP</a></li>
				<li><a><span class='ui-icon ui-icon-pin-w'></span>Assign BAC Id</a>
					<ul style='z-index: 1000;white-space: nowrap;'>
						<li><a onclick='openDialogForm(\"jobAssignBacForm.cgi\",\"jobList$$\")'><span class='ui-icon ui-icon-bullet'></span>Speedy Mode</a></li>
						<li><a onclick='openDialogForm(\"jobForceBacForm.cgi\",\"jobList$$\")'><span class='ui-icon ui-icon-bullet'></span>Forcing Mode</a></li>
					</ul>
				</li>
				<li><a onclick='refresh(\"general\")'><span class='ui-icon ui-icon-refresh'></span>Refresh</a></li>
			</ul>
		</li>
	</ul>";

unless($jobs)
{
	$button .= "<p class='ui-state-error ui-corner-all' style='padding: .7em;'><span class='ui-icon ui-icon-alert' style='float: left; margin-right: .3em;'></span>
			<strong>No job, please create one!</strong></p>";
}

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$button/$button/g;
$html =~ s/\$jobs/$jobs/g;
$html =~ s/\$\$/$$/g;

print header(-cookie=>cookie(-name=>'general',-value=>3));
print $html;

__DATA__
$button
$jobs
<script>
buttonInit();
loadingHide();
$( "#jobs$$" ).dataTable({
	"scrollY": "600px",
	"scrollCollapse": true,
	"paging": false,
	"order": [ 1, 'desc' ],
	"columnDefs": [
    { "orderable": false, "targets": 0 }
  ]
});
$( "#jobInfoMenu$$" ).menu();
</script>