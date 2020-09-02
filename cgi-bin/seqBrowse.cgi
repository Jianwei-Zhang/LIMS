#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser );
use JSON::XS; #JSON::XS is recommended to be installed for handling JSON string of big size
use Bio::SeqIO;
use DBI;
use SVG;
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

my $svg;

my $seqId = param ('seqId') || '';
my $scrollLeft = param ('scrollLeft') || '0';
my $highlight = param ('highlight') || '';
my $seqDetails = '';
my $dialogWidth = 600;
my $textFontSize = 10;
my $textFontWidth = 7;
my $pixelUnit = 500;
my $barY = 25;
my $rulerY = 20;
my $margin = 4;
my $barHeight = 12;
my $hiddenSeqPosition = 50;
my $barSpacing = 400; #space between reference and ctgs
my $unitLength = 10000;
my $maxCol = 0;
my $totalLength;
if ($seqId)
{
	my $refSequence=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$refSequence->execute($seqId);
	my @refSequence = $refSequence->fetchrow_array();
	my $refSequenceDetails = decode_json $refSequence[8];
	$refSequenceDetails->{'id'} = '' unless (exists $refSequenceDetails->{'id'});
	$refSequenceDetails->{'description'} = '' unless (exists $refSequenceDetails->{'description'});
	$refSequenceDetails->{'gapList'} = '' unless (exists $refSequenceDetails->{'gapList'});
	$refSequenceDetails->{'filter'} = '' unless (exists $refSequenceDetails->{'filter'});
	my $totalSeqs = 0;


	# create an SVG object
	$svg= SVG->new(width=>'$svgWidth',height=>'$svgHeight');
	$svg->rectangle(
		id    => $$,
		style => { stroke => 'white',
					fill => 'white'
					},
		x     => 0,
		y     => 0,
		width => '$svgWidth',
		height=> '$svgHeight'
	);


	# use explicit element constructor to generate a group element
    my $ruler=$svg->group(
        id    => 'ruler',
        style => { stroke=>'black'}
    );

	# use explicit element constructor to generate a group element
    my $refSeq=$svg->group(
        id    => 'refSeq'
    );
    my $refSeqGap=$svg->group(
        id    => 'refSeqGap'
    );
    my $refSeqFilter=$svg->group(
        id    => 'refSeqFilter'
    );

	if ($refSequence[5] < $pixelUnit*500)
	{
		$pixelUnit = $refSequence[5]/500;
	}

	$refSeq->rectangle(
			x     => $margin,
			y     => $barY,
			width => $refSequence[5] / $pixelUnit,
			height=> $barHeight,
			style => { stroke => 'black',
						fill => 'lightgreen'
					},
			onclick => "closeDialog();openDialog('seqView.cgi?seqId=$refSequence[0]')",
			id    => "refSeq$refSequence[0]"
		);
	if($refSequenceDetails->{'gapList'})
	{
		foreach (split ",", $refSequenceDetails->{'gapList'})
		{
			my ($gapStart,$gapEnd) = split "-", $_;
			$refSeqGap->rectangle(
				x     => $margin + $gapStart / $pixelUnit,
				y     => $barY,
				width => ($gapEnd - $gapStart + 1) / $pixelUnit,
				height=> $barHeight,
				style => { stroke => 'red',
							fill => 'red'
					},
				id    => "refSeqGap$gapStart$gapEnd"
			);
		}
	}
	if($refSequenceDetails->{'filter'})
	{
		foreach (split ",", $refSequenceDetails->{'filter'})
		{
			my ($filterStart,$filterEnd) = split "-", $_;
			$refSeqFilter->rectangle(
				x     => $margin + $filterStart / $pixelUnit,
				y     => $barY,
				width => ($filterEnd - $filterStart + 1) / $pixelUnit,
				height=> $barHeight,
				style => { stroke => 'black',
							fill => 'lightgrey'
					},
				id    => "refSeqFilter$filterStart$filterEnd"
			);
		}
	}

	#list BES
    my $besSeq=$svg->group(
        id    => 'besSeq'
    );
	my $besLeftPosition;
	my $besRightPosition;
	my $besLeftDirection;
	my $besRightDirection;
	my $besLeftAlignment;
	my $besRightAlignment;
	my $fpcContig;
	my $fpcCloneLeftEnd;
	my $fpcCloneRightEnd;
	my @lengthList;
	my $totalLength = 0;
	my $besId = 0;
	open (BES,">$commoncfg->{TMPDIR}/BES-$refSequence[2].txt") or die "can't open file: $commoncfg->{TMPDIR}/BES-$refSequence[2].txt";
	my $besList=$dbh->prepare("SELECT * FROM alignment WHERE subject = ? AND program LIKE 'BES%' AND perc_indentity >= 95 ORDER BY s_start");
	$besList->execute($seqId);
	while (my @besList = $besList->fetchrow_array())
	{
		$besId++;
		my $besBarX = ($besList[11] > $besList[10]) ? $margin + $besList[10] / $pixelUnit : $margin + $besList[11] / $pixelUnit;
		my $besBarY = $barY + $barHeight*1.5;
		$besSeq->rectangle(
			id    => $besList[2].$besId,
			x     => $besBarX,
			y     => $besBarY,
			width => ($besList[11] > $besList[10]) ? ($besList[11] - $besList[10]) / $pixelUnit : ($besList[10] - $besList[11]) / $pixelUnit,
			height=> $barHeight,
			style => { stroke => ($besList[11] > $besList[10]) ? 'black' : 'red',
						fill => 'blue',
						'fill-opacity' => 0.3
					}
		);
		my $besSequence=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$besSequence->execute($besList[2]);
		my @besSequence = $besSequence->fetchrow_array();

		$fpcContig->{$besSequence[2]} = "None" unless (exists $fpcContig->{$besSequence[2]});
		$fpcCloneLeftEnd->{$besSequence[2]} = -1 unless (exists $fpcCloneLeftEnd->{$besSequence[2]});
		$fpcCloneRightEnd->{$besSequence[2]} = -1 unless (exists $fpcCloneRightEnd->{$besSequence[2]});

		my $getFpcClone = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'fpcClone' AND name LIKE ?");
		$getFpcClone->execute($besSequence[2]);
		while (my @getFpcClone = $getFpcClone->fetchrow_array())
		{
			$fpcContig->{$besSequence[2]} = 'Ctg0';
			$fpcCloneLeftEnd->{$besSequence[2]} = 0;
			$fpcCloneRightEnd->{$besSequence[2]} = 0;
			if ($getFpcClone[8] =~ /Map "(.*)" Ends Left (\d*)/)
			{
				$fpcContig->{$besSequence[2]} = ucfirst ($1);
				$fpcCloneLeftEnd->{$besSequence[2]} = $2;
			}
			if ($getFpcClone[8] =~ /Ends Right (\d*)/)
			{
				$fpcCloneRightEnd->{$besSequence[2]} = $1;
			}
		}

		if (exists $besLeftPosition->{$besSequence[2]})
		{
			my $besDistance = $besList[10] - $besLeftPosition->{$besSequence[2]};
 			next if($besDistance > 300000 || $besDistance < 25000);
 			next if($besLeftDirection->{$besSequence[2]} == $besSequence[6]);
			push @lengthList,$besDistance;
			$besRightDirection->{$besSequence[2]} = $besSequence[6];
			$besRightPosition->{$besSequence[2]} = ($besList[11] > $besList[10]) ? $besList[11] : $besList[10];
			$besRightAlignment->{$besSequence[2]} = ($besList[11] > $besList[10]) ? "+" : "-";
			$totalLength += $besDistance;
			if($besLeftAlignment->{$besSequence[2]} eq $besRightAlignment->{$besSequence[2]})
			{
				print BES "$besSequence[2]\t$refSequence[2]\t$besLeftPosition->{$besSequence[2]}\t$besDistance\t$seqDir{$besLeftDirection->{$besSequence[2]}}\t$seqDir{$besSequence[6]}\t$besLeftAlignment->{$besSequence[2]}\t=\t$fpcContig->{$besSequence[2]}\t$fpcCloneLeftEnd->{$besSequence[2]}\t$fpcCloneRightEnd->{$besSequence[2]}\n";
			}
			else
			{
				print BES "$besSequence[2]\t$refSequence[2]\t$besLeftPosition->{$besSequence[2]}\t$besDistance\t$seqDir{$besLeftDirection->{$besSequence[2]}}\t$seqDir{$besSequence[6]}\t$besLeftAlignment->{$besSequence[2]}\t$besRightAlignment->{$besSequence[2]}\t$fpcContig->{$besSequence[2]}\t$fpcCloneLeftEnd->{$besSequence[2]}\t$fpcCloneRightEnd->{$besSequence[2]}\n";
			}
		}
		else
		{
			$besLeftPosition->{$besSequence[2]} = ($besList[11] > $besList[10]) ? $besList[10] : $besList[11];
			$besLeftDirection->{$besSequence[2]} = $besSequence[6];
			$besLeftAlignment->{$besSequence[2]} = ($besList[11] > $besList[10]) ? "+" : "-";
		}
	}
	close (BES);
	`gzip -f '$commoncfg->{TMPDIR}/BES-$refSequence[2].txt'`;
	my $besEvaluation = "<a href='$commoncfg->{TMPURL}/BES-$refSequence[2].txt.gz' target='hiddenFrame'><span class='ui-icon ui-icon-bullet'></span>BES Evaluation</a>" if (-e "$commoncfg->{TMPDIR}/BES-$refSequence[2].txt.gz");


	my @besCloneList = sort { $besLeftPosition->{$a} <=> $besLeftPosition->{$b} } keys %$besRightPosition;

    my $besClone=$svg->group(
        id    => 'besClone'
    );
	my $colCount=0;
	my $toBeFlipped = 0;
	my @colPostion;
	push @colPostion, -1;
    for my $currentClone (@besCloneList)
    {
		next if ($fpcCloneRightEnd->{$currentClone} ne '-1');
    	my $col = 0;
    	my $goodCol = @colPostion;
    	for(@colPostion)
    	{
    		if($_ < $besLeftPosition->{$currentClone})
    		{
    			$goodCol = $col;
    			last;
    		}
    		$col++;
    	}
		$besClone->anchor(
			id      => 'clone'.$currentClone.$colCount,
			onclick => "closeDialog();openDialog('cloneView.cgi?cloneName=$currentClone')",
			title   => $currentClone
			)->rectangle(
				x     => $margin + $besLeftPosition->{$currentClone} / $pixelUnit,
				y     => $barY + $barHeight * (1.3 * ($goodCol + 2) + 1),
				width => ($besRightPosition->{$currentClone} - $besLeftPosition->{$currentClone}) / $pixelUnit,
				height=> $barHeight,
				style => { stroke => ($besLeftAlignment->{$currentClone} eq $besRightAlignment->{$currentClone}) ? 'red' :'grey',
							'stroke-dasharray' => '3,1',
							fill => ($besLeftAlignment->{$currentClone} eq $besRightAlignment->{$currentClone}) ? 'yellow': 'white'
							},
				id    => "besClone$currentClone$colCount$$"
			);
		$toBeFlipped++ if ($besLeftAlignment->{$currentClone} eq $besRightAlignment->{$currentClone});
		my $besDistanceText =  commify($besRightPosition->{$currentClone} - $besLeftPosition->{$currentClone});
		my $lable = ($fpcCloneRightEnd->{$currentClone} eq '-1') ?  $currentClone." (~$besDistanceText bp)" : $currentClone." (~$besDistanceText bp) FPC: $fpcContig->{$currentClone}($fpcCloneLeftEnd->{$currentClone}-$fpcCloneRightEnd->{$currentClone})";
		my $lableLength = length $lable;

    	my $textX = $margin + $barHeight + $besLeftPosition->{$currentClone} / $pixelUnit;
    	my $textY = $barY + $barHeight * (1.3 * ($goodCol + 2 ) + 2) - 2;
    	my $textEnd = $besLeftPosition->{$currentClone}  + $lableLength * $textFontWidth * $pixelUnit;
		$besClone->text(
			id      => 'cloneName'.$currentClone.$colCount,
			onclick => "closeDialog();openDialog('cloneView.cgi?cloneName=$currentClone')",
			x       => $textX,
			y       => $textY,
			style   => {
				'font-size'   => $textFontSize,
				'stroke'        => 'black'
			}
		)->cdata($lable);
    	$colPostion[$goodCol] = ($textEnd > $besRightPosition->{$currentClone}) ? $textEnd : $besRightPosition->{$currentClone};
		$colCount++;
    }

	my $nonFpcCloneNumber = $colCount;
	my $nonFpcCloneToBeFlippedNumber = $toBeFlipped;

	my $colNumber = 0;
    for(@colPostion)
	{
		$colPostion[$colNumber] = $refSequence[5];
    	$colNumber++;
    }

	#zone line
	$ruler->line(
		id    => "fpcCloneOrNot",
		style => {
			stroke => 'blue',
			'stroke-dasharray' => '3,3'
			},
		x1    => $margin,
		y1    => $barY + $barHeight * (1.3 * ($colNumber + 2) + 1) - 2,
		x2    => $margin +  $refSequence[5] / $pixelUnit,
		y2    => $barY + $barHeight * (1.3 * ($colNumber + 2) + 1) - 2
	);
	#nonFpcClone
	$ruler->text(
		id      => "nonFpcCloneNumber",
		x       => $margin,
		y       => $barY + $barHeight * (1.3 * ($colNumber + 2) + 1) - 2,
		style   => {
			'font-size'   =>  11,
			'stroke'        => 'blue'
		}
	)->cdata("$nonFpcCloneNumber ($nonFpcCloneToBeFlippedNumber conflict) non-FPC clones");

    for my $currentClone (@besCloneList)
    {
		next if ($fpcCloneRightEnd->{$currentClone} eq '-1');
    	my $col = 0;
    	my $goodCol = @colPostion;
    	for(@colPostion)
    	{
    		if($_ < $besLeftPosition->{$currentClone})
    		{
    			$goodCol = $col;
    			last;
    		}
    		$col++;
    	}
		$besClone->anchor(
			id      => 'clone'.$currentClone.$colCount,
			onclick => "closeDialog();openDialog('cloneView.cgi?cloneName=$currentClone')",
			title   => $currentClone
			)->rectangle(
				x     => $margin + $besLeftPosition->{$currentClone} / $pixelUnit,
				y     => $barY + $barHeight * (1.3 * ($goodCol + 2) + 1),
				width => ($besRightPosition->{$currentClone} - $besLeftPosition->{$currentClone}) / $pixelUnit,
				height=> $barHeight,
				style => { stroke => ($besLeftAlignment->{$currentClone} eq $besRightAlignment->{$currentClone}) ? 'red' :'grey',
							'stroke-dasharray' => '3,1',
							fill => ($besLeftAlignment->{$currentClone} eq $besRightAlignment->{$currentClone}) ? 'yellow': 'white'
							},
				id    => "besClone$currentClone$colCount$$"
			);
		$toBeFlipped++ if ($besLeftAlignment->{$currentClone} eq $besRightAlignment->{$currentClone});
		my $besDistanceText =  commify($besRightPosition->{$currentClone} - $besLeftPosition->{$currentClone});
		my $lable = ($fpcCloneRightEnd->{$currentClone} eq '0') ?  $currentClone." (~$besDistanceText bp)" : $currentClone." (~$besDistanceText bp) FPC: $fpcContig->{$currentClone}($fpcCloneLeftEnd->{$currentClone}-$fpcCloneRightEnd->{$currentClone})";
		my $lableLength = length $lable;

    	my $textX = $margin + $barHeight + $besLeftPosition->{$currentClone} / $pixelUnit;
    	my $textY = $barY + $barHeight * (1.3 * ($goodCol + 2 ) + 2) - 2;
    	my $textEnd = $besLeftPosition->{$currentClone}  + $lableLength * $textFontWidth * $pixelUnit;
		$besClone->text(
			id      => 'cloneName'.$currentClone.$colCount,
			onclick => "closeDialog();openDialog('cloneView.cgi?cloneName=$currentClone')",
			x       => $textX,
			y       => $textY,
			style   => {
				'font-size'   => $textFontSize,
				'stroke'        => ($fpcCloneRightEnd->{$currentClone} eq '0') ? 'red' : 'black'
			}
		)->cdata($lable);
    	$colPostion[$goodCol] = ($textEnd > $besRightPosition->{$currentClone}) ? $textEnd : $besRightPosition->{$currentClone};
		$colCount++;
    }
	my $fpcCloneNumber = $colCount - $nonFpcCloneNumber;
	my $fpcCloneToBeFlippedNumber = $toBeFlipped - $nonFpcCloneToBeFlippedNumber;

	#fpcClone
	$ruler->text(
		id      => "fpcCloneNumber",
		x       => $margin,
		y       => $barY + $barHeight * (1.3 * (@colPostion + 2) + 1) - 2,
		style   => {
			'font-size'   =>  11,
			'stroke'        => 'green'
		}
	)->cdata("$fpcCloneNumber ($fpcCloneToBeFlippedNumber conflict) FPC clones");


	@lengthList = sort {$b <=> $a} @lengthList;

	my $medianLength = 0;
	if($#lengthList % 2 == 1)
	{
		my $median = int ($#lengthList/2);
		$medianLength = ($lengthList[$median]+$lengthList[$median+1])/2;
	}
	else
	{
		my $median = $#lengthList/2;
		$medianLength = $lengthList[$median];
	}

	my $n50Length = 0;
	my $subtotal = 0;
	foreach (@lengthList)
	{
		$subtotal += $_;
		if($subtotal > $totalLength/2 && $n50Length == 0)
		{
			$n50Length = $_;
		}
	}

	for (my $rulerNumber = 0;$rulerNumber <= $refSequence[5];$rulerNumber += $unitLength)
	{
		#dash lines
		$ruler->line(
			id    => "rulerDash$rulerNumber",
			style => {
				stroke=> ($rulerNumber % (5*$unitLength) == 0) ? 'grey' : 'lightgrey',
				'stroke-dasharray' => '3,3'
				},
			x1    => $margin + $rulerNumber / $pixelUnit,
			y1    => $rulerY,
			x2    => $margin + $rulerNumber / $pixelUnit,
			y2    => '$svgHeight'
		);
		#region
		$ruler->line(
			id    => "rulerRegion$rulerNumber",
			x1    => $margin + $rulerNumber / $pixelUnit,
			y1    => $rulerY,
			x2    => $margin +  ($rulerNumber + $unitLength ) / $pixelUnit,
			y2    => $rulerY
		);
		#ticks
		$ruler->line(
			id    => "rulerTick$rulerNumber",
			x1    => $margin + $rulerNumber / $pixelUnit,
			y1    => ($rulerNumber % (5*$unitLength) == 0) ? $rulerY - 5 : $rulerY - 3,
			x2    => $margin + $rulerNumber / $pixelUnit,
			y2    => $rulerY

		);
		if($rulerNumber % (5*$unitLength) == 0)
		{
			my $commifiedRulerNumber = commify($rulerNumber);
			$commifiedRulerNumber =~ s/,000$/k/g;
			$commifiedRulerNumber =~ s/,000k$/M/g;
			my $textX =$margin + $rulerNumber / $pixelUnit;
			$ruler->text(
				id      => "rulerNumber$rulerNumber",
				x       => $textX,
				y       => $rulerY - 2,
				style   => {
					'font-size'   =>  11,
					'stroke'        => ($commifiedRulerNumber =~ /M$/ ) ? 'red' : 'black'
				}
			)->cdata($commifiedRulerNumber);
		}
	}

	# now render the SVG object, implicitly use svg namespace
	$seqDetails = $svg->xmlify;
	my $svgWidth = $refSequence[5] / $pixelUnit + $margin * 2;
	my $svgHeight = $barY + $barHeight * (1.3 * (@colPostion +2) + 2) + $margin;
	$dialogWidth = ($svgWidth > 1000 ) ? 1050 : ($svgWidth < 550) ? 600 : $svgWidth + 50;
	$seqDetails =~ s/\$svgWidth/$svgWidth/g;
	$seqDetails =~ s/\$svgHeight/$svgHeight/g;
	$seqDetails =~ s/&/&amp;/g;

	open (SVGFILE,">$commoncfg->{TMPDIR}/BES-$refSequence[2].svg") or die "can't open file: $commoncfg->{TMPDIR}/BES-$refSequence[2].svg";
	print SVGFILE $seqDetails;
	close (SVGFILE);
	`gzip -f '$commoncfg->{TMPDIR}/BES-$refSequence[2].svg'`;
	my $besSvg = "<a href='$commoncfg->{TMPURL}/BES-$refSequence[2].svg.gz' target='hiddenFrame'><span class='ui-icon ui-icon-bullet'></span>SVG</a>" if (-e "$commoncfg->{TMPDIR}/BES-$refSequence[2].svg.gz");

	$html =~ s/\$besEvaluation/$besEvaluation/g;
	$html =~ s/\$besSvg/$besSvg/g;
	$html =~ s/\$seqDetails/$seqDetails/g;
	$html =~ s/\$dialogWidth/$dialogWidth/g;
	$html =~ s/\$seqName/$refSequence[2]/g;
	$html =~ s/\$scrollLeft/$scrollLeft/g;
	$html =~ s/\$\$/$$/g;

	print header();
	print $html;
}
else
{
	print header(-type=>'text/html',-status=>'402 Invalid operation');
	exit;
}

__DATA__

<ul class='seqBrowserMenu' style='width: 100px;float: left; margin-right: .3em; white-space: nowrap;'>
	<li><a><span class='ui-icon ui-icon-disk'></span>Download</a>
		<ul style='z-index: 1000;white-space: nowrap;'>
			<li>$besEvaluation</li>
			<li>$besSvg</li>
		</ul>
	</li>
</ul>
<br><br>
$seqDetails
<script>
$( ".seqBrowserMenu" ).menu();
$('#viewer').dialog("option", "title", "Seq Browser: $seqName");
$('#viewer').dialog("option", "width", "$dialogWidth");
$('#viewer').scrollLeft($scrollLeft);
</script>
