#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use JSON::XS; #JSON::XS is recommended to be installed for handling JSON string of big size 
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

undef $/;# enable slurp mode
my $html = <DATA>;

my %seqDir = (
	0=>'NA',
	1=>'f',
	2=>'r'
	);

if(param ('libraryId'))
{
	my $libraryId = param ('libraryId');
	my $targetGenome;
	my @targetGenome;
	my $besToGenome;
	my $besIdToGenome;
	my $genomeList=$dbh->prepare("SELECT * FROM matrix WHERE z = ? AND container LIKE 'genome'");
	$genomeList->execute($libraryId);
	while (my @genomeList = $genomeList->fetchrow_array())
	{
		$targetGenome->{$genomeList[0]} = $genomeList[2];
		push @targetGenome,$genomeList[0];

		my $refSequence=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 99 AND x = ? ORDER BY y DESC");
		$refSequence->execute($genomeList[0]);
		while (my @refSequence = $refSequence->fetchrow_array())
		{
			my $refSequenceDetails = decode_json $refSequence[8];
			$refSequenceDetails->{'id'} = '' unless (exists $refSequenceDetails->{'id'});
			$refSequenceDetails->{'description'} = '' unless (exists $refSequenceDetails->{'description'});
			$refSequenceDetails->{'gapList'} = '' unless (exists $refSequenceDetails->{'gapList'});
			$refSequenceDetails->{'filter'} = '' unless (exists $refSequenceDetails->{'filter'});
			my $besLeftPosition;
			my $besRightPosition;
			my $besLeftDirection;
			my $besRightDirection;
			my $besLeftAlignment;
			my $besRightAlignment;
			my $fpcView;
			my $fpcCloneLeftEnd;
			my $fpcCloneRightEnd;
			my $besList=$dbh->prepare("SELECT * FROM alignment WHERE subject = ? AND program LIKE 'BES%' AND perc_indentity >= 95 ORDER BY s_start");
			$besList->execute($refSequence[0]);
			while (my @besList = $besList->fetchrow_array())
			{
				$besIdToGenome->{$genomeList[0]}->{$besList[2]} = 1;
				my $besSequence=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
				$besSequence->execute($besList[2]);
				my @besSequence = $besSequence->fetchrow_array();

				$fpcView->{$besSequence[2]} = "None" unless (exists $fpcView->{$besSequence[2]});
				$fpcCloneLeftEnd->{$besSequence[2]} = -1 unless (exists $fpcCloneLeftEnd->{$besSequence[2]});
				$fpcCloneRightEnd->{$besSequence[2]} = -1 unless (exists $fpcCloneRightEnd->{$besSequence[2]});
				my $getFpcClone = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'fpcClone' AND name LIKE ?");
				$getFpcClone->execute($besSequence[2]);
				while (my @getFpcClone = $getFpcClone->fetchrow_array())
				{
					$fpcView->{$besSequence[2]} = 'Ctg0';
					$fpcCloneLeftEnd->{$besSequence[2]} = 0;
					$fpcCloneRightEnd->{$besSequence[2]} = 0;
					if ($getFpcClone[8] =~ /Map "(.*)" Ends Left (\d*)/)
					{
						$fpcView->{$besSequence[2]} = ucfirst ($1);
						$fpcCloneLeftEnd->{$besSequence[2]} = $2;
					}
					if ($getFpcClone[8] =~ /Ends Right (\d*)/)
					{
						$fpcCloneRightEnd->{$besSequence[2]} = $1;
					}
				}

				if (exists $besLeftPosition->{$refSequence[0]}->{$besSequence[2]})
				{
					my $besDistance = $besList[10] - $besLeftPosition->{$refSequence[0]}->{$besSequence[2]};
					next if($besDistance > 300000 || $besDistance < 25000);
					next if($besLeftDirection->{$refSequence[0]}->{$besSequence[2]} == $besSequence[6]);
					$besRightDirection->{$refSequence[0]}->{$besSequence[2]} = $besSequence[6];
					$besRightPosition->{$refSequence[0]}->{$besSequence[2]} = ($besList[11] > $besList[10]) ? $besList[11] : $besList[10];
					$besRightAlignment->{$refSequence[0]}->{$besSequence[2]} = ($besList[11] > $besList[10]) ? "+" : "-";
					$besToGenome->{$genomeList[0]}->{$besSequence[2]} = ($besLeftAlignment->{$refSequence[0]}->{$besSequence[2]} eq $besRightAlignment->{$refSequence[0]}->{$besSequence[2]}) ? "$refSequence[2]\t$besLeftPosition->{$refSequence[0]}->{$besSequence[2]}\t$besDistance\t$seqDir{$besLeftDirection->{$refSequence[0]}->{$besSequence[2]}}\t$seqDir{$besSequence[6]}\t$besLeftAlignment->{$refSequence[0]}->{$besSequence[2]}\t=\t$fpcView->{$besSequence[2]}\t$fpcCloneLeftEnd->{$besSequence[2]}\t$fpcCloneRightEnd->{$besSequence[2]}" :
						"$refSequence[2]\t$besLeftPosition->{$refSequence[0]}->{$besSequence[2]}\t$besDistance\t$seqDir{$besLeftDirection->{$refSequence[0]}->{$besSequence[2]}}\t$seqDir{$besSequence[6]}\t$besLeftAlignment->{$refSequence[0]}->{$besSequence[2]}\t$besRightAlignment->{$refSequence[0]}->{$besSequence[2]}\t$fpcView->{$besSequence[2]}\t$fpcCloneLeftEnd->{$besSequence[2]}\t$fpcCloneRightEnd->{$besSequence[2]}";
				}
				else
				{
					$besLeftPosition->{$refSequence[0]}->{$besSequence[2]} = ($besList[11] > $besList[10]) ? $besList[10] : $besList[11];
					$besLeftDirection->{$refSequence[0]}->{$besSequence[2]} = $besSequence[6];
					$besLeftAlignment->{$refSequence[0]}->{$besSequence[2]} = ($besList[11] > $besList[10]) ? "+" : "-";
				}
			}
		}
	}

	my $besClone;
	my $paredBes;
	my $singleEndCount;
	my $getSequences = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 98 AND x = ? ORDER BY name,id");
	$getSequences->execute($libraryId);
	while(my @getSequences =  $getSequences->fetchrow_array())
	{
		my $sequenceDetails = decode_json $getSequences[8];
		$sequenceDetails->{'id'} = '' unless (exists $sequenceDetails->{'id'});
		$sequenceDetails->{'description'} = '' unless (exists $sequenceDetails->{'description'});
		$sequenceDetails->{'gapList'} = '' unless (exists $sequenceDetails->{'gapList'});
		$besClone->{$getSequences[2]} .= "$seqDir{$getSequences[6]}:$getSequences[5] ";
		$paredBes->{$getSequences[2]} .= $seqDir{$getSequences[6]};
		$singleEndCount->{$seqDir{$getSequences[6]}}++;
	}
	my $singleEnd;
	for (sort keys %$singleEndCount)
	{
		$singleEnd .= ($singleEnd) ? " $_:$singleEndCount->{$_}" : "$_:$singleEndCount->{$_}";
	}

	my $totalClone = keys %$besClone;
	my $totalBes = $getSequences->rows; 
	$totalBes = "$totalBes ($singleEnd) BESs";
	my $besDetails = "<table id='bes$$' class='display' style='width: 100%;'>
			<thead>
				<tr style='text-align:left'>
					<th><b>Clone Name</b></th>
					<th><b>Direction:Length</b> (bp)</th>
					<th><b>BAC Hits</b></th>";
	for(@targetGenome)
	{
			$besDetails .= "<th><b>$targetGenome->{$_}</b></th>";
	}
	$besDetails .= "</tr>
			</thead>
			<tbody>";
	my $paredBesNumber = 0;
	my $paredBesMapped;
	for (sort keys %$besClone)
	{
		my $besCloneName = $_;
		my $targetGenomeMatch = '';
		for(@targetGenome)
		{
				$targetGenomeMatch .= (exists $besToGenome->{$_}->{$besCloneName})? "<td>$besToGenome->{$_}->{$besCloneName}</td>" : "<td></td>";
				$paredBesMapped->{$_}++ if (exists $besToGenome->{$_}->{$besCloneName});
		}
		$besDetails .= "<tr><td>$besCloneName</td><td>$besClone->{$besCloneName}</td><td></td>$targetGenomeMatch</tr>";
		$paredBesNumber++ if ($paredBes->{$besCloneName} =~ /fr/ || $paredBes->{$besCloneName} =~ /rf/);
	}	

	my $mappedPairs = '';
	$besDetails .= "</tbody><tfoot style='text-align:left'><tr><th>$totalClone clones</th><th>$totalBes, $paredBesNumber pairs</th><th></th>"; 
	for(@targetGenome)
	{
		my $targetGenomeId = $_;
		my $besIdToGenomeTotal = scalar (keys %{$besIdToGenome->{$targetGenomeId}});

		my $besToBeFlipped;
		open (BES,">$commoncfg->{TMPDIR}/BES-$targetGenome->{$targetGenomeId}.txt") or die "can't open file: $commoncfg->{TMPDIR}/BES-$targetGenome->{$targetGenomeId}.txt";
		for (sort keys %{$besToGenome->{$targetGenomeId}})
		{
			print BES "$_\t$besToGenome->{$targetGenomeId}->{$_}\n";
			$besToBeFlipped->{$targetGenomeId}++ if ($besToGenome->{$targetGenomeId}->{$_} =~ /=/);
		}
		close (BES);
		`gzip -f $commoncfg->{TMPDIR}/BES-$targetGenome->{$targetGenomeId}.txt`;
		$besDetails .= "<th>Mapped: $besIdToGenomeTotal ($paredBesMapped->{$targetGenomeId} paired, including $besToBeFlipped->{$targetGenomeId} conflict ones)</th>";
		$totalBes .= ($mappedPairs) ? " and $besIdToGenomeTotal<sup>$targetGenome->{$targetGenomeId}</sup>" : ": $besIdToGenomeTotal<sup>$targetGenome->{$targetGenomeId}</sup>";
		$mappedPairs .= ($mappedPairs) ? " and $paredBesMapped->{$targetGenomeId} ($besToBeFlipped->{$targetGenomeId} conflict)<a href='$commoncfg->{TMPURL}/BES-$targetGenome->{$targetGenomeId}.txt.gz' target='hiddenFrame'><sup>$targetGenome->{$targetGenomeId}</sup></a>" : "$paredBesMapped->{$targetGenomeId} ($besToBeFlipped->{$targetGenomeId} conflict)<a href='$commoncfg->{TMPURL}/BES-$targetGenome->{$targetGenomeId}.txt.gz' target='hiddenFrame'><sup>$targetGenome->{$targetGenomeId}</sup></a>";
	}
	$besDetails .= "</tr>"; 
	$besDetails .= "</tfoot></table>"; 

	open (BESTEXT,">$commoncfg->{TMPDIR}/BES-report.$libraryId.html") or die "can't open file: $commoncfg->{TMPDIR}/BES-report.$libraryId.html";
	print BESTEXT $besDetails;
	close (BESTEXT);
	`gzip -f $commoncfg->{TMPDIR}/BES-report.$libraryId.html`;
	my $besReport = "<a href='$commoncfg->{TMPURL}/BES-report.$libraryId.html.gz' target='hiddenFrame'>Download Details</a>" if (-e "$commoncfg->{TMPDIR}/BES-report.$libraryId.html.gz");
	$html =~ s/\$besReport/$besReport/g;
	$html =~ s/\$totalClone/$totalClone/g;
	$html =~ s/\$totalBes/$totalBes/g;
	$html =~ s/\$paredBesNumber/$paredBesNumber/g;
	$html =~ s/\$mappedPairs/$mappedPairs/g;
	$html =~ s/\$\$/$$/g;
	print header;
	print $html;
}
else
{
	print header(-type=>'text/html',-status=>'402 Invalid operation');
}


__DATA__
<div id="besReport$$"  name="besReport$$">
$totalClone clones <br><br>
$totalBes mapped. <br><br>
$paredBesNumber pairs: $mappedPairs mapped. <br><br>

($besReport)
</div>
<script>
buttonInit();
$('#dialog').dialog("option", "title", "BES Report");
$( "#dialog" ).dialog( "option", "buttons", [ { text: "Print", click: function() {printDiv('besReport$$'); } }, { text: "OK", click: function() {closeDialog(); } } ] );
</script>