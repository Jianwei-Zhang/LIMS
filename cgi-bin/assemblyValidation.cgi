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
my $validateDetail = '';
my $noSequenceTab = '';
my $nullAssemblySeqTab = '';
my $errorAssemblySeqTab = '';
my $noAssemblySeqTab = '';
my $repeatedAssemblySeqTab = '';
my $assemblySeqWithoutCtgTab = '';
my $nullAssemblySeqDetails= '';
my $errorAssemblySeqDetails= '';
my $noAssemblySeqDetails = '';
my $repeatedAssemblySeqDetails = '';
my $assemblySeqWithoutCtgDetails = '';
my $noSequenceDetails = '';
my %assemblySeqToCtg;
my $newSeqIdList;
my %seqCount;
my $seq;
my %seqInList;


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

	my $noSequenceNum = 0;
	my $nullAssemblySeqNum = 0;
	my $errorAssemblySeqNum = 0;
	my $noAssemblySeqNum = 0;
	my $repeatedAssemblySeqNum = 0;
	my $assemblySeqWithoutCtgNum = 0;

	my $emptyAssemblyCtg = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblyCtg' AND note LIKE '' AND o = ?");
	$emptyAssemblyCtg->execute($assemblyId);
#	$noAssemblySeqNum = $emptyAssemblyCtg->rows;
	while (my @emptyAssemblyCtg = $emptyAssemblyCtg->fetchrow_array())
	{
		$noAssemblySeqDetails .= "<tr><td>$emptyAssemblyCtg[0]</td><td>Ctg$emptyAssemblyCtg[2]</td></tr>";
		$noAssemblySeqNum++;
	}

	my $validateSeqToCtg=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblyCtg' AND note NOT LIKE '' AND o = ? ORDER BY name");
	$validateSeqToCtg->execute($assemblyId);
	while(my @validateSeqToCtg = $validateSeqToCtg->fetchrow_array())
	{
		foreach (split ",", $validateSeqToCtg[8])
		{
			$nullAssemblySeqNum++ unless ($_=~/\S/);
			$nullAssemblySeqDetails .= "<tr><td>Ctg$validateSeqToCtg[2]</td><td>$validateSeqToCtg[8]</td></tr>" unless ($_=~/\S/);
			next unless ($_=~/\S/);
			$errorAssemblySeqDetails .= "<tr><td>Ctg$validateSeqToCtg[2]</td><td>$_</td><td>$validateSeqToCtg[8]</td></tr>" if($_=~/^(\-{2,})/);
			$errorAssemblySeqNum++ if($_=~/^(\-{2,})/);
			$seq=$_;
			$_ =~ s/[^a-zA-Z0-9]//g;
			$seqInList{$_}{$validateSeqToCtg[2]}=$seq;
			$seqCount{$_}+=1;
			$assemblySeqToCtg{$_}{$validateSeqToCtg[2]}=1;
			my $assemblySeq=$dbh->prepare("SELECT * FROM matrix WHERE id = ? ");
			$assemblySeq->execute($_);			
			my @assemblySeq = $assemblySeq->fetchrow_array();
			my $assemblySequence=$dbh->prepare("SELECT * FROM matrix WHERE id = ? ");
			$assemblySequence->execute($assemblySeq[5]);
			$noSequenceDetails .= "<tr><td>Ctg$validateSeqToCtg[2]</td><td>$assemblySeq[2]</td><td>sequence$assemblySeq[5]</td></tr>" unless ($assemblySequence->rows > 0);
			$noSequenceNum++ unless ($assemblySequence->rows > 0);
		}
	}
	#whether assemblySeq appears only once
	my @seqInCtg=keys %seqCount;
	for my $seqInCtg(@seqInCtg)
	{
		if ($seqCount{$seqInCtg}>1)
		{
			$repeatedAssemblySeqDetails .= "<tr><td>$seqInCtg </td><td>$seqCount{$seqInCtg} times</td><td>";
			$repeatedAssemblySeqNum++;
			my @ctgName=keys %{$assemblySeqToCtg{$seqInCtg}};
			for my $ctgName(@ctgName)
			{
				$repeatedAssemblySeqDetails .= "[Ctg$ctgName $seqInList{$seqInCtg}{$ctgName}]<br />";
			}
			$repeatedAssemblySeqDetails .="</td></tr>";
		}
	}
	
	#whether assemblySeq appears in the assemblyCtg
	my $validateAssemblySeq=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblySeq' AND o = ?");
	$validateAssemblySeq->execute($assemblyId);
	while (my @validateAssemblySeq = $validateAssemblySeq->fetchrow_array())
	{
		$assemblySeqWithoutCtgDetails .= "<tr><td> $validateAssemblySeq[0]</td></tr>" unless (exists $seqCount{$validateAssemblySeq[0]});
		$assemblySeqWithoutCtgNum++ unless (exists $seqCount{$validateAssemblySeq[0]});
	}

	$validateDetail = "<table>
	<tr><td style='text-align:right'><b>No sequence</b></td><td>$noSequenceNum</td></tr>
	<tr><td style='text-align:right'><b>Null assemblySeq</b></td><td>$nullAssemblySeqNum</td></tr>
	<tr><td style='text-align:right'><b>Error assemblySeq</b></td><td>$errorAssemblySeqNum</td></tr>
	<tr><td style='text-align:right'><b>Empty assemblyCtg</b></td><td>$noAssemblySeqNum</td></tr>
	<tr><td style='text-align:right'><b>Repeated assemblySeq</b></td><td>$repeatedAssemblySeqNum</td></tr>
	<tr><td style='text-align:right'><b>AssemblySeq without Ctg</b></td><td>$assemblySeqWithoutCtgNum</td></tr>
	</table>";
	if ($noSequenceNum > 0)
	{
		$noSequenceTab = "<li><a href='#viewValidateTabs-noSequence'>No sequence ($noSequenceNum)</a></li>";
		$noSequenceDetails = "<div id='viewValidateTabs-noSequence'>
		<table class='validateDetails display'><thead><tr><th>Ctg</th><th>Assembly Seq</th><th>No Sequence</th></tr></thead><tbody>$noSequenceDetails</tbody></table>
		</div>";
	}

	if ($nullAssemblySeqNum > 0)
	{
		$nullAssemblySeqTab = "<li><a href='#viewValidateTabs-nullAssemblySeq'>Null assemblySeq ($nullAssemblySeqNum)</a></li>";
		$nullAssemblySeqDetails = "<div id='viewValidateTabs-nullAssemblySeq'>
		<table class='validateDetails display'><thead><tr><th>Ctg</th><th>Assembly Seq List</th></tr></thead><tbody>$nullAssemblySeqDetails</tbody></table>
		</div>";
	}
	if ($errorAssemblySeqNum > 0)
	{
		$errorAssemblySeqTab = "<li><a href='#viewValidateTabs-errorAssemblySeq'>Error assemblySeq ($errorAssemblySeqNum)</a></li>";
		$errorAssemblySeqDetails = "<div id='viewValidateTabs-errorAssemblySeq'>
		<table class='validateDetails display'><thead><tr><th>Ctg</th><th>Error Assembly Seq</th><th>Assembly Seq List</th></tr></thead><tbody>$errorAssemblySeqDetails</tbody></table>
		</div>";
	}
	if ($noAssemblySeqNum > 0)
	{
		$noAssemblySeqTab = "<li><a href='#viewValidateTabs-noAssemblySeq'>Empty assemblyCtg ($noAssemblySeqNum)</a></li>";
		$noAssemblySeqDetails = "<div id='viewValidateTabs-noAssemblySeq'>
		<table class='validateDetails display'><thead><tr><th>Empty assemblyCtg (ids)</th><th>Ctg</th></tr></thead><tbody>$noAssemblySeqDetails</tbody></table>
		</div>";
	}
	if ($repeatedAssemblySeqNum > 0)
	{
		$repeatedAssemblySeqTab = "<li><a href='#viewValidateTabs-repeatedAssemblySeq'>Repeated assemblySeq ($repeatedAssemblySeqNum)</a></li>";
		$repeatedAssemblySeqDetails = "<div id='viewValidateTabs-repeatedAssemblySeq'>
		<table class='validateDetails display'><thead><tr><th>Assembly Seq</th><th>Numbers</th><th>Ctgs</th></tr></thead><tbody>$repeatedAssemblySeqDetails</tbody></table>
		</div>";
	}
	if ($assemblySeqWithoutCtgNum > 0)
	{
		$assemblySeqWithoutCtgTab = "<li><a href='#viewValidateTabs-assemblySeqWithoutCtg'>AssemblySeq without Ctg ($assemblySeqWithoutCtgNum)</a></li>";
		$assemblySeqWithoutCtgDetails = "<div id='viewValidateTabs-assemblySeqWithoutCtg'>
		<table class='validateDetails display'><thead><tr><th>Assembly seq does not appeare in the assemblyCtg </th></tr></thead><tbody>$assemblySeqWithoutCtgDetails</tbody></table>
		</div>";
	}

    #print "$assemblySeqWithoutCtg\n";
	undef $/;# enable slurp mode
	my $html = <DATA>;
    $html =~ s/\$\$/$$/g;
	$html =~ s/\$noSequenceNum/$noSequenceNum/g;
	$html =~ s/\$noSequenceTab/$noSequenceTab/g;
	$html =~ s/\$noSequenceDetails/$noSequenceDetails/g;
    $html =~ s/\$validateDetail/$validateDetail/g;
	$html =~ s/\$nullAssemblySeqTab/$nullAssemblySeqTab/g;
	$html =~ s/\$errorAssemblySeqTab/$errorAssemblySeqTab/g;	
	$html =~ s/\$noAssemblySeqTab/$noAssemblySeqTab/g;
	$html =~ s/\$repeatedAssemblySeqTab/$repeatedAssemblySeqTab/g;
	$html =~ s/\$assemblySeqWithoutCtgTab/$assemblySeqWithoutCtgTab/g;
	$html =~ s/\$nullAssemblySeqNum/$nullAssemblySeqNum/g;
	$html =~ s/\$errorAssemblySeqNum/$errorAssemblySeqNum/g;	
	$html =~ s/\$noAssemblySeqNum/$noAssemblySeqNum/g;
	$html =~ s/\$repeatedAssemblySeqNum/$repeatedAssemblySeqNum/g;
	$html =~ s/\$assemblySeqWithoutCtgNum/$assemblySeqWithoutCtgNum/g;
	$html =~ s/\$nullAssemblySeqDetails/$nullAssemblySeqDetails/g;
	$html =~ s/\$errorAssemblySeqDetails/$errorAssemblySeqDetails/g;	
	$html =~ s/\$noAssemblySeqDetails/$noAssemblySeqDetails/g;
	$html =~ s/\$repeatedAssemblySeqDetails/$repeatedAssemblySeqDetails/g;
	$html =~ s/\$assemblySeqWithoutCtgDetails/$assemblySeqWithoutCtgDetails/g;
	$html =~ s/\$assemblyId/$assemblyId/g;
	$html =~ s/\$commoncfg->{HTDOCS}/$commoncfg->{HTDOCS}/g;
	
	print $html;
}

__DATA__
<form id="assemblyCorrect" name="assemblyCorrect" action="assemblyCorrect.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
<input name="assemblyId" id="assemblyId$$" type="hidden" value="$assemblyId" />
</form>
<div id="viewValidateTabs$$" name="viewValidateTabs$$">
	<ul>
	<li><a href="#viewValidateTabs-Num">Summary</a></li>
	$noSequenceTab
	$nullAssemblySeqTab
	$errorAssemblySeqTab
	$noAssemblySeqTab
	$repeatedAssemblySeqTab
	$assemblySeqWithoutCtgTab
	</ul>
	<div id="viewValidateTabs-Num">
	$validateDetail
	</div>
	$noSequenceDetails
	$nullAssemblySeqDetails
	$errorAssemblySeqDetails
	$noAssemblySeqDetails
	$repeatedAssemblySeqDetails
	$assemblySeqWithoutCtgDetails
</div>
<script>
buttonInit();
$( "#viewValidateTabs$$" ).tabs();
$( ".validateDetails" ).dataTable({
	"paging": false,
	"searching": false,
	"info": false,
	"columnDefs": [
    { "orderable": false, "targets": 0 }
  ]
});

$('#dialog').dialog("option", "title", "Validate Assembly");
if ($noSequenceNum > 0 || $nullAssemblySeqNum > 0 || $errorAssemblySeqNum > 0 || $noAssemblySeqNum > 0 || $repeatedAssemblySeqNum > 0 || $assemblySeqWithoutCtgNum > 0)
{
	$( "#dialog" ).dialog( "option", "buttons", [{ text: "Validate", click: function() { submitForm('assemblyCorrect'); } },{ text: "OK", click: function() { closeDialog(); } } ] );
}
else
{
	$( "#dialog" ).dialog( "option", "buttons", [{ text: "OK", click: function() { closeDialog(); } } ] );
}
</script>
