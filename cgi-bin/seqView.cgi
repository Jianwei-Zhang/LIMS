#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use JSON::XS;
use Bio::SeqIO;
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

my $seqId = param ('seqId') || '';

my %seqType = (
	0=>'Assembled',
	1=>'BAC-Insert',
	2=>'BAC-Circularized',
	3=>'BAC-NonVector',
	4=>'BAC-Gapped',
	5=>'Partial',
	6=>'Vector/Mixer',
	7=>'Mixer',
	8=>'SHORT',
	97=>'Piece',
	98=>'BES',
	99=>'Genome'
	);

my %bacAssignType = (
	0=>'',
	1=>'TagValid',
	2=>'BesValid',
	3=>'TagValid+BesValid',
	4=>'TagForced'
	);

my %seqDir = (
	0=>'NA',
	1=>'f',
	2=>'r'
	);

print header;
if ($seqId)
{
	my $getSeq = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$getSeq->execute($seqId);
	if($getSeq->rows < 1)
	{
		print 'No valid seq Found!';
		exit;
	}
	else
	{
		my @getSeq = $getSeq->fetchrow_array();
		my $cloneLink = '';
		my $targetLink = '';
		my $assemblyView = '';
		my $sequenceDetails = decode_json $getSeq[8];
		$sequenceDetails->{'id'} = '' unless (exists $sequenceDetails->{'id'});
		$sequenceDetails->{'description'} = '' unless (exists $sequenceDetails->{'description'});
		my $descriptionLength = length ($sequenceDetails->{'description'});
		if ($descriptionLength > 30)
		{
			$sequenceDetails->{'description'} = "<a title='$sequenceDetails->{'description'}'>". substr($sequenceDetails->{'description'},0,20). "...". substr($sequenceDetails->{'description'},-5). "</a>";
		}
		$sequenceDetails->{'gapList'} = '' unless (exists $sequenceDetails->{'gapList'});
		$sequenceDetails->{'filter'} = '' unless (exists $sequenceDetails->{'filter'});
		my @gaps = split (",",$sequenceDetails->{'gapList'});
		my $longSeq = 1200;
		my $sequence = 'ERROR: NO SEQUENCE FOUND! PLEASE CONTACT YOUR ADMINISTRATOR.';
		my $in = Bio::SeqIO->new(-file => "$commoncfg->{DATADIR}/$sequenceDetails->{'sequence'}",
								-format => 'Fasta');
		while ( my $seq = $in->next_seq() )
		{
			$sequence = $seq->seq;
		}
		my $sequenceLeft = ($getSeq[5] > $longSeq) ? substr($sequence,0,$longSeq/2) : $sequence;
		my $sequenceRight =  ($getSeq[5] > $longSeq) ? substr($sequence,-$longSeq/2) : '';
		$sequence = multiLineSeq($sequenceLeft,60);
		if ($sequenceRight)
		{
			$sequence .= "...\n";
			$sequence .= multiLineSeq($sequenceRight,60);
		}
		$sequence .= "(<a href='download.cgi?seqId=$getSeq[0]' title='Download this sequence' target='hiddenFrame'>".commify($getSeq[5])." bp<span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-disk'></span></a>)" if ($getSeq[5] > $longSeq);

		if($getSeq[3] < 50) #all pacBio Sequences
		{
			$targetLink = "<a style='float: right;' onclick='closeDialog();openDialog(\"jobView.cgi?jobId=$getSeq[4]\")'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-arrow-1-ne' title='View'></span>Related Job</a>";
			my $getClone = $dbh->prepare("SELECT * FROM clones WHERE name LIKE ?");
			$getClone->execute($getSeq[2]);
			$cloneLink = "<a style='float: right;' onclick='closeDialog();openDialog(\"cloneView.cgi?cloneName=$getSeq[2]\")'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-arrow-1-ne' title='View'></span>Related Clone</a>" if($getClone->rows > 0);
		}
		else
		{
			my $target=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
			$target->execute($getSeq[4]);
			my @target = $target->fetchrow_array();
			if($target[1] eq 'library')
			{
				my $getClone = $dbh->prepare("SELECT * FROM clones WHERE name LIKE ?");
				$getClone->execute($getSeq[2]);
				$cloneLink = "<a style='float: right;' onclick='closeDialog();openDialog(\"cloneView.cgi?cloneName=$getSeq[2]\")'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-arrow-1-ne' title='View'></span>Related Clone</a>";
			}
			if($target[1] eq 'genome')
			{
				$targetLink = "<a style='float: right;' onclick='closeDialog();openDialog(\"genomeView.cgi?genomeId=$getSeq[4]\")'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-arrow-1-ne' title='View'></span>Related Genome</a>";
			}
		}

		my $assemblyViewTabs;
		my $assemblySeq = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblySeq' AND y = ?");
		$assemblySeq->execute($seqId);
		while (my @assemblySeq = $assemblySeq->fetchrow_array())
		{
			$assemblyView = "<tr><td style='text-align:right'><b>Assembly</b></td><td><div id='assemblySeqTabs$$'><ul>" unless ($assemblyView);
			my $assemblyCtg=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblyCtg' AND MATCH (note) AGAINST (?)");
			$assemblyCtg->execute($assemblySeq[0]);
			my @assemblyCtg = $assemblyCtg->fetchrow_array();

			my $assembly=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
			$assembly->execute($assemblyCtg[3]);
			my @assembly = $assembly->fetchrow_array();
			
			my $assemblySeqList = '';
			my $i = 1;
			for (split ",", $assemblyCtg[8])
			{
				my $hidden = ($_ =~ /^-/) ? 1 : 0;
				$_ =~ s/[^a-zA-Z0-9]//g;
				my $assemblySeqMember = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
				$assemblySeqMember->execute($_);
				my @assemblySeqMember = $assemblySeqMember->fetchrow_array();
				my $assemblySeqOrientation = ($assemblySeqMember[7] < 0) ? "<sup title='Orientation: -'>-</sup>" : "<sup title='Orientation: +'>+</sup>";
				$assemblySeqList .= ($_ eq $assemblySeq[0]) ? ($hidden) ? "$i.<a class='ui-state-error-text' onclick='closeDialog();openDialog(\"seqView.cgi?seqId=$assemblySeqMember[5]\")' title='SeqId:$assemblySeqMember[5]'>$assemblySeqMember[2]</a> " : "$i.<a class='ui-state-error-text' onclick='closeDialog();openDialog(\"seqView.cgi?seqId=$assemblySeqMember[5]\")' title='SeqId:$assemblySeqMember[5]'>$assemblySeqMember[2]</a>$assemblySeqOrientation "
					: ($hidden) ? "$i.<a onclick='closeDialog();openDialog(\"seqView.cgi?seqId=$assemblySeqMember[5]\")' title='SeqId:$assemblySeqMember[5]' class='ui-state-disabled'>$assemblySeqMember[2]</a> " : "$i.<a onclick='closeDialog();openDialog(\"seqView.cgi?seqId=$assemblySeqMember[5]\")' title='SeqId:$assemblySeqMember[5]'>$assemblySeqMember[2]</a>$assemblySeqOrientation ";
				$i++;
			}
			
			$assemblyView .= "<li><a href='#assemblyTabs-$assembly[0]-$assemblyCtg[0]'>$assembly[2] v.$assembly[3] Ctg$assemblyCtg[2]</a></li>";
			$assemblyViewTabs .= "<div id='assemblyTabs-$assembly[0]-$assemblyCtg[0]'>
				<a style='float: right;' onclick='closeDialog();closeViewer();openViewer(\"assemblyCtgView.cgi?assemblyCtgId=$assemblyCtg[0]&highlight=$assemblySeq[0]\")'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-image' title='View'></span>View Ctg$assemblyCtg[2]</a><br>
				Sequence List:<br>$assemblySeqList</div>";
		}
		$assemblyView .= "</ul>$assemblyViewTabs<div></td></tr>" if ($assemblyView);
		$getSeq[2] = 'no name' unless ($getSeq[2]); 
		$getSeq[5] = commify($getSeq[5]);
		
		undef $/;# enable slurp mode
		my $html = <DATA>;
		$html =~ s/\$seqId/$getSeq[0]/g;
		$html =~ s/\$seqName/$getSeq[2]/g;
		$html =~ s/\$seqType/$seqType{$getSeq[3]}/g;
		if(@gaps)
		{
			my $gaps = @gaps;
			if($gaps > 1)
			{
				$html =~ s/\$gaps/$gaps gaps,/g;
			}
			else
			{
				$html =~ s/\$gaps/$gaps gap,/g;
			}
		}
		else
		{
			$html =~ s/\$gaps//g;
		}
		my $assignChr = '';
		if($getSeq[3] == 99)
		{
			$assignChr = "<tr><td style='text-align:right;white-space: nowrap;'><b>Chromosome</b></td><td><div id='seqChr$seqId' style='position: relative;'><a id='seqChr$seqId$$' onmouseover='editIconShow(\"seqChr$seqId$$\")' onmouseout='editIconHide(\"seqChr$seqId$$\")' onclick='loaddiv(\"seqChr$seqId\",\"seqChrEdit.cgi?seqId=$seqId\")' title='Edit this chromosome number'>$getSeq[6]</a></div></td></tr>";
		}
		
		$html =~ s/\$assignChr/$assignChr/g;
		$html =~ s/\$targetLink/$targetLink/g;
		$html =~ s/\$cloneLink/$cloneLink/g;
		$html =~ s/\$seqLength/$getSeq[5]/g;
		$html =~ s/\$assemblyView/$assemblyView/g;
		my $filter;
		if ($sequenceDetails->{'filter'})
		{
			$filter=$sequenceDetails->{'filter'};
		}
		else
		{
			$filter='none';
		}
		$html =~ s/\$filter/$filter/g;
		if($getSeq[2] eq $sequenceDetails->{'id'})
		{
			$html =~ s/\$sequence/>$getSeq[2] $sequenceDetails->{'description'}\n$sequence/g;
		}
		else
		{
			$html =~ s/\$sequence/>$getSeq[2]<span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-circlesmall-plus' title='Original ID: $sequenceDetails->{'id'}'><\/span> $sequenceDetails->{'description'}\n$sequence/g;
		}
		$html =~ s/\$\$/$$/g;
		$html =~ s/\$seqId/$seqId/g;
		$html =~ s/\$commoncfg->{HTDOCS}/$commoncfg->{HTDOCS}/g;

		print $html;
	}
}
else
{
	print 'No valid seq Found!';
	exit;
}

__DATA__
<div id="viewSeq$$" name="viewSeq$$">
	<table width='100%'>
	<tr><td style='text-align:right;white-space: nowrap;'><b>Seq Name</b></td><td><div id='seqName$seqId' style='position: relative;'><a id='seqName$seqId$$' onmouseover='editIconShow("seqName$seqId$$")' onmouseout='editIconHide("seqName$seqId$$")' onclick="loaddiv('seqName$seqId','seqNameEdit.cgi?seqId=$seqId')" title="Edit this name">$seqName</a></div>($seqType) $gaps $seqLength bp $cloneLink $targetLink</td></tr>
	$assignChr
	$assemblyView
	<tr><td style='text-align:right;'><b>Sequence</b></td><td><PRE>$sequence</PRE></td></tr>
	<tr><td style='text-align:right;white-space: nowrap;'><b>Filter region(s)</b><br><sub>(e.g. 120-130,134-201)</sub></td><td><div id='seqFilter$seqId' style='position: relative;'><a id='seqFilter$seqId$$' onmouseover='editIconShow("seqFilter$seqId$$")' onmouseout='editIconHide("seqFilter$seqId$$")' onclick="loaddiv('seqFilter$seqId','seqFilterEdit.cgi?seqId=$seqId')" title="Edit this filter">$filter</a></div></td></tr>
	</table>
</div>
<script>
buttonInit();
$( "#assemblySeqTabs$$" ).tabs();
$('#dialog').dialog("option", "title", "View Sequence $seqName");
$( "#dialog" ).dialog( "option", "buttons", [ { text: "Browse", click: function() {closeDialog(); closeViewer();openViewer("seqBrowse.cgi?seqId=$seqId")} }, { text: "Delete", click: function() { deleteItem($seqId); } }, { text: "Print", click: function() {printDiv('viewSeq$$'); } }, { text: "OK", click: function() {closeDialog(); } } ] );
</script>
