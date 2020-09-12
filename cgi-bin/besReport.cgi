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

my %seqDir = (
	0=>'NA',
	1=>'f',
	2=>'r'
	);

if(param ('libraryId'))
{
	my $libraryId = param ('libraryId');
	my $besClone;
	my $paredBes;
	my $singleEndCount;
	my $getSequences = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 98 AND x = ?");
	$getSequences->execute($libraryId);
	while(my @getSequences =  $getSequences->fetchrow_array())
	{
		my $sequenceDetails = decode_json $getSequences[8];
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
		my $queryDirLibrary;
		for (my $position = 0; $position < length($libraryId); $position += 2)
		{
			$queryDirLibrary .= "/q". substr($libraryId,$position,2);
		}
		my $subjectDirGenome;
		for (my $position = 0; $position < length($genomeList[0]); $position += 2)
		{
			$subjectDirGenome .= "/s". substr($genomeList[0],$position,2);
		}
		if (-e "$commoncfg->{DATADIR}/alignments/setToSet$queryDirLibrary$subjectDirGenome/$libraryId-$genomeList[0].list")
		{
			my $besLeftPosition;
			my $besRightPosition;
			my $besLeftDirection;
			my $besRightDirection;
			my $besLeftAlignment;
			my $besRightAlignment;
			my $fpcView;
			my $fpcCloneLeftEnd;
			my $fpcCloneRightEnd;
			open (LIST,"$commoncfg->{DATADIR}/alignments/setToSet$queryDirLibrary$subjectDirGenome/$libraryId-$genomeList[0].list") or die "can't open file: $commoncfg->{DATADIR}/alignments/setToSet$queryDirLibrary$subjectDirGenome/$libraryId-$genomeList[0].list";
			while (<LIST>)
			{
				chop;
				if(-e "$commoncfg->{DATADIR}/$_")
				{
					open (TBL, "$commoncfg->{DATADIR}/$_") or die "can't open file: $commoncfg->{DATADIR}/$_";
					while(<TBL>)
					{
						chop;
						/^#/ and next;
						my @besList = split("\t",$_);
						next if ($besList[12] > 0);
						next if ($besList[2] < 95);

						$besIdToGenome->{$genomeList[0]}->{$besList[0]} = 1;
						my $besSequence=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
						$besSequence->execute($besList[0]);
						my @besSequence = $besSequence->fetchrow_array();
						my $refSequence=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
						$refSequence->execute($besList[1]);
						my @refSequence = $refSequence->fetchrow_array();

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

						if (exists $besLeftPosition->{$besList[1]}->{$besSequence[2]})
						{
							my $besDistance = $besList[8] - $besLeftPosition->{$besList[1]}->{$besSequence[2]};
							next if($besDistance > 300000 || $besDistance < 25000);
							next if($besLeftDirection->{$besList[1]}->{$besSequence[2]} == $besSequence[6]);
							$besRightDirection->{$besList[1]}->{$besSequence[2]} = $besSequence[6];
							$besRightPosition->{$besList[1]}->{$besSequence[2]} = ($besList[9] > $besList[8]) ? $besList[9] : $besList[8];
							$besRightAlignment->{$besList[1]}->{$besSequence[2]} = ($besList[9] > $besList[8]) ? "+" : "-";
							$besToGenome->{$genomeList[0]}->{$besSequence[2]} = ($besLeftAlignment->{$besList[1]}->{$besSequence[2]} eq $besRightAlignment->{$besList[1]}->{$besSequence[2]}) ? "$refSequence[2]\t$besLeftPosition->{$besList[1]}->{$besSequence[2]}\t$besDistance\t$seqDir{$besLeftDirection->{$besList[1]}->{$besSequence[2]}}\t$seqDir{$besSequence[6]}\t$besLeftAlignment->{$besList[1]}->{$besSequence[2]}\t=\t$fpcView->{$besSequence[2]}\t$fpcCloneLeftEnd->{$besSequence[2]}\t$fpcCloneRightEnd->{$besSequence[2]}" :
								"$refSequence[2]\t$besLeftPosition->{$besList[1]}->{$besSequence[2]}\t$besDistance\t$seqDir{$besLeftDirection->{$besList[1]}->{$besSequence[2]}}\t$seqDir{$besSequence[6]}\t$besLeftAlignment->{$besList[1]}->{$besSequence[2]}\t$besRightAlignment->{$besList[1]}->{$besSequence[2]}\t$fpcView->{$besSequence[2]}\t$fpcCloneLeftEnd->{$besSequence[2]}\t$fpcCloneRightEnd->{$besSequence[2]}";
						}
						else
						{
							$besLeftPosition->{$besList[1]}->{$besSequence[2]} = ($besList[9] > $besList[8]) ? $besList[8] : $besList[9];
							$besLeftDirection->{$besList[1]}->{$besSequence[2]} = $besSequence[6];
							$besLeftAlignment->{$besList[1]}->{$besSequence[2]} = ($besList[9] > $besList[8]) ? "+" : "-";
						}
					}
					close(TBL);
				}
			}
			close(LIST);
		}
	}

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

	undef $/;# enable slurp mode
	my $html = <DATA>;
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