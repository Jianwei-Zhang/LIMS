#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use JSON::XS; #JSON::XS is recommended to be installed for handling JSON string of big size 
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

my $svg;
my $assemblyId = param ('assemblyId') || '';
my $cookieCompanionAssemblyId = cookie("companionAssembly$assemblyId") || '';
my $companionAssemblyId = param ('companionAssemblyId') || $cookieCompanionAssemblyId;
$companionAssemblyId = '' if ($companionAssemblyId eq $assemblyId);
my $chr = param ('chr') || '0';
my $refChr = $chr % 100; #get ref chr number for polyploid genomes.
my $scrollLeft = param ('scrollLeft') || '0';
my $highlight = param ('highlight') || '';
my $assemblyChrDetails = '';
my $dialogWidth = 600;
my $textFontSize = 10;
my $textFontWidth = 7;
my $cookiePixelUnit = cookie("pixelUnit") || '500';
my $pixelUnit = param ('pixelUnit') || $cookiePixelUnit;
my %pixelUnitOptions = ( 500=>'500 bp/pixel', 5000=>'5000 bp/pixel');
my $pixelUnitList = '';
for (sort {$a <=>$b} keys %pixelUnitOptions)
{
	$pixelUnitList .= ($_ eq $pixelUnit) ? "<option value='$_' selected>$pixelUnitOptions{$_}</option>":"<option value='$_'>$pixelUnitOptions{$_}</option>";
}
my $barY = 25;
my $rulerY = 20;
my $margin = 4;
my $barHeight = 12;
my $hiddenSeqPosition = 50;
my $barSpacing = ($companionAssemblyId) ?  300: 400; #space between reference and ctgs
my $cookieTickLength = cookie("tickLength") || '10000';
my $tickLength = param ('tickLength') || $cookieTickLength;
my %tickLengthOptions = ( 10000=>'10kb', 100000=>'100kb');
my $tickLengthList;
for (sort {$a <=>$b} keys %tickLengthOptions)
{
	$tickLengthList .= ($_ eq $tickLength) ? "<option value='$_' selected>$tickLengthOptions{$_}</option>":"<option value='$_'>$tickLengthOptions{$_}</option>";
}
my $cookieAlignmentLength = cookie("alignmentLength") || '1000';
my $alignmentLength = param ('alignmentLength') || $cookieAlignmentLength;
my %alignmentLengthOptions = ( 1000=>'1kb', 10000=>'10kb', 100000=>'100kb');
my $alignmentLengthList;
for (sort {$a <=>$b} keys %alignmentLengthOptions)
{
	$alignmentLengthList .= ($_ eq $alignmentLength) ? "<option value='$_' selected>$alignmentLengthOptions{$_}</option>":"<option value='$_'>$alignmentLengthOptions{$_}</option>";
}
my $maxCol = 0;
my $totalLength;
if ($assemblyId && $chr)
{
	my $assembly=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$assembly->execute($assemblyId);
	my @assembly = $assembly->fetchrow_array();
	my $companionAssemblyList = "<option value='$assemblyId'>None</option>";
	my $companionAssembly=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assembly' AND y = ? AND id != ?");
	$companionAssembly->execute($assembly[5],$assemblyId);
	while (my @companionAssembly = $companionAssembly->fetchrow_array())
	{
		$companionAssemblyList .= ($companionAssemblyId eq $companionAssembly[0]) ? "<option value='$companionAssembly[0]' selected>$companionAssembly[2].$companionAssembly[3]</option>":"<option value='$companionAssembly[0]'>$companionAssembly[2].$companionAssembly[3]</option>";
	}
	my $refGenome=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$refGenome->execute($assembly[5]);
	my @refGenome = $refGenome->fetchrow_array();
	my $refGenomeSequence=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 99 AND x = ? AND z = ?");
	$refGenomeSequence->execute($assembly[5],$refChr);
	my @refGenomeSequence = $refGenomeSequence->fetchrow_array();

	my $refSequenceDetails = decode_json $refGenomeSequence[8];
	$refSequenceDetails->{'id'} = '' unless (exists $refSequenceDetails->{'id'});
	$refSequenceDetails->{'description'} = '' unless (exists $refSequenceDetails->{'description'});
	$refSequenceDetails->{'gapList'} = '' unless (exists $refSequenceDetails->{'gapList'});
	$refSequenceDetails->{'filter'} = '' unless (exists $refSequenceDetails->{'filter'});
	my $totalSeqs = 0;
	my $spaceForCompanion = ($companionAssemblyId) ? $barSpacing + $barHeight * 3 : 0;
	my $svgHeight = $barY + $barSpacing + $barHeight * 3 + $spaceForCompanion + $margin;

	# create an SVG object
	$svg= SVG->new(width=>'$svgWidth',height=>$svgHeight);
	$svg->rectangle(
		id    => $$,
		class   => 'hasmenuForAll',
		style => { stroke => 'white',
					fill => 'white'
					},
		x     => 0,
		y     => 0,
		width => '$svgWidth',
		height=> $svgHeight,
		'blast-url' => "blastTwoseqForm.cgi?assemblyId=$assemblyId",
		'align-url' => "alignmentCheckForm.cgi?assemblyId=$assemblyId",
		'fill-url' => "assemblyGapFillForm.cgi?assemblyId=$assemblyId"
	);

	# use explicit element constructor to generate a group element
    my $ruler=$svg->group(
        id    => 'ruler',
        style => { stroke=>'black'}
    );

	# use explicit element constructor to generate a group element
    my $assemblyRefChr=$svg->group(
        id    => 'assemblyRefChr'
    );
    my $assemblyRefChrGap=$svg->group(
        id    => 'assemblyRefChrGap'
    );
    my $assemblyRefChrFilter=$svg->group(
        id    => 'assemblyRefChrFilter'
    );

	$assemblyRefChr->rectangle(
			x     => $margin,
			y     => $barY + $spaceForCompanion,
			width => $refGenomeSequence[5] / $pixelUnit,
			height=> $barHeight,
			style => { stroke => 'black',
						fill => 'lightgreen'
					},
			onclick => "closeDialog();openDialog('seqView.cgi?seqId=$refGenomeSequence[0]')",
			id    => "assemblyRefChr$refGenomeSequence[0]"
		);

	$assemblyRefChr->text(
		id      => 'refName'.$refGenomeSequence[0],
		onclick => "closeDialog();openDialog('seqView.cgi?seqId=$refGenomeSequence[0]')",
		x       => $margin,
		y       => $barY + $spaceForCompanion + $barHeight - 2,
		style   => {
			'font-size'   => $textFontSize,
			'stroke'        =>  'black'
		}
	)->cdata("$refGenome[2]-$refGenomeSequence[2]");

	if($refSequenceDetails->{'gapList'})
	{
		foreach (split ",", $refSequenceDetails->{'gapList'})
		{
			my ($gapStart,$gapEnd) = split "-", $_;
			$assemblyRefChrGap->rectangle(
				x     => $margin + $gapStart / $pixelUnit,
				y     => $barY + $spaceForCompanion,
				width => ($gapEnd - $gapStart + 1) / $pixelUnit,
				height=> $barHeight,
				style => { stroke => 'red',
							fill => 'red'
					},
				id    => "assemblyRefChrGap$gapStart$gapEnd"
			);
		}
	}
	if($refSequenceDetails->{'filter'})
	{
		foreach (split ",", $refSequenceDetails->{'filter'})
		{
			my ($filterStart,$filterEnd) = split "-", $_;
			$assemblyRefChrFilter->rectangle(
				x     => $margin + $filterStart / $pixelUnit,
				y     => $barY + $spaceForCompanion,
				width => ($filterEnd - $filterStart + 1) / $pixelUnit,
				height=> $barHeight,
				style => { stroke => 'black',
							fill => 'lightgrey'
					},
				id    => "assemblyRefChrFilter$filterStart$filterEnd"
			);
		}
	}

    my $assemblyChrCtg=$svg->group(
        id    => 'assemblyChrCtg'
    );
    my $assemblyChrSeq=$svg->group(
        id    => 'assemblyChrSeq'
    );

    my $assemblyChrSeqAlignment=$svg->group(
        id    => 'assemblyChrSeqAlignment'
    );

    my $assemblyChrSeqHidden=$svg->group(
        id    => 'assemblyChrSeqHidden'
    );
	if($companionAssemblyId)
	{
		my $companionAssembly=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$companionAssembly->execute($companionAssemblyId);
		my @companionAssembly = $companionAssembly->fetchrow_array();
		my $companionAssemblyChrCtg=$svg->group(
			id    => 'companionAssemblyChrCtg'
		);
		my $companionAssemblyChrSeq=$svg->group(
			id    => 'companionAssemblyChrSeq'
		);

		my $companionAssemblyChrSeqAlignment=$svg->group(
			id    => 'companionAssemblyChrSeqAlignment'
		);

		my $companionAssemblyChrSeqHidden=$svg->group(
			id    => 'companionAssemblyChrSeqHidden'
		);

		my @ctgPosition;
		my $companionAssemblyCtgList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblyCtg' AND o = ? AND x = ? ORDER BY y");
		$companionAssemblyCtgList->execute($companionAssemblyId,$chr);
		while (my @companionAssemblyCtgList = $companionAssemblyCtgList->fetchrow_array())
		{
			$totalLength->{'companionAssembly'} = 0 unless (exists $totalLength->{'companionAssembly'});
			$totalLength->{'companionAssembly'} += $companionAssemblyCtgList[7];
			my $ctgRow = @ctgPosition;
			my $j = 0;
			for (@ctgPosition)
			{
				if ($_ < $margin + $companionAssemblyCtgList[6] / $pixelUnit)
				{
					$ctgRow = $j;
					last;
				}
				$j++;
			}
			$ctgPosition[$ctgRow] = $margin + ($companionAssemblyCtgList[6] + $companionAssemblyCtgList[7]) / $pixelUnit;

			my $companionAssemblySequenceName;
			my $companionAssemblySequenceId;
			my $companionAssemblySequenceData;
			my $companionAssemblySequenceType;
			my $companionAssemblySeqLength;
			my $companionAssemblySeqOrient;
			my $companionAssemblySeqStart;
			my $companionAssemblySeqEnd;
			my $companionAssemblySeqOrder;
			my $companionAssemblySeqHide;
			my $companionAssemblySeqFpcCtg;
			my $order = 0;
			my $filterLength = 0;
			foreach (split ",", $companionAssemblyCtgList[8])
			{
				next unless ($_);
				my $hide = ($_ =~ /^-/) ? 1 : 0;
				$_ =~ s/[^a-zA-Z0-9]//g;
				$companionAssemblySeqHide->{$_} = $hide;
				$order++;
				$companionAssemblySeqOrder->{$_} = $order;
				my $companionAssemblySeqList=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
				$companionAssemblySeqList->execute($_);
				my @companionAssemblySeqList = $companionAssemblySeqList->fetchrow_array();
				$companionAssemblySequenceName->{$companionAssemblySeqList[0]} = $companionAssemblySeqList[2];
				$companionAssemblySequenceId->{$companionAssemblySeqList[0]} = $companionAssemblySeqList[5];
				$companionAssemblySeqLength->{$companionAssemblySeqList[0]} = $companionAssemblySeqList[6];
				$companionAssemblySeqOrient->{$companionAssemblySeqList[0]} = ($companionAssemblySeqList[7] < 0) ? "-" : "+";
				if($companionAssemblySeqList[8])
				{
					($companionAssemblySeqStart->{$companionAssemblySeqList[0]},$companionAssemblySeqEnd->{$companionAssemblySeqList[0]}) = split ",",$companionAssemblySeqList[8];
				}
				else
				{
					$companionAssemblySeqStart->{$companionAssemblySeqList[0]} = 1;
					$companionAssemblySeqEnd->{$companionAssemblySeqList[0]} = $companionAssemblySeqList[6];
					my $updateAssemblySeq=$dbh->do("UPDATE matrix SET note = '$companionAssemblySeqStart->{$companionAssemblySeqList[0]},$companionAssemblySeqEnd->{$companionAssemblySeqList[0]}' WHERE id = $companionAssemblySeqList[0]");
				}
				my $companionAssemblySequence=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
				$companionAssemblySequence->execute($companionAssemblySeqList[5]);
				my @companionAssemblySequence = $companionAssemblySequence->fetchrow_array();
				$companionAssemblySequenceType->{$companionAssemblySeqList[0]} = $companionAssemblySequence[3];
				$companionAssemblySequenceData->{$companionAssemblySeqList[0]} = $companionAssemblySequence[8];
				$companionAssemblySeqFpcCtg->{$companionAssemblySeqList[0]} = "NA";
				my $fpcCloneList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'fpcClone' AND o = ? AND name LIKE ?");
				$fpcCloneList->execute($companionAssembly[6],$companionAssemblySeqList[2]);
				while (my @fpcCloneList = $fpcCloneList->fetchrow_array())
				{
					$fpcCloneList[8] =~ /Map "(.*)"/;
					$companionAssemblySeqFpcCtg->{$companionAssemblySeqList[0]} = $1;
				}
				my $companionSequenceDetails = decode_json $companionAssemblySequence[8];
				$companionSequenceDetails->{'filter'} = '' unless (exists $companionSequenceDetails->{'filter'});
				if ($companionSequenceDetails->{'filter'}) 
				{
					foreach (split ",", $companionSequenceDetails->{'filter'} )
					{
						my ($filterStart,$filterEnd) = split "-", $_;
						next if ($companionAssemblySeqStart->{$companionAssemblySeqList[0]} > $filterEnd);
						next if ($companionAssemblySeqEnd->{$companionAssemblySeqList[0]} < $filterStart);
						if ($companionAssemblySeqStart->{$companionAssemblySeqList[0]} >= $filterStart && $companionAssemblySeqEnd->{$companionAssemblySeqList[0]} <= $filterEnd)
						{
							$filterLength += $companionAssemblySeqEnd->{$companionAssemblySeqList[0]} - $companionAssemblySeqStart->{$companionAssemblySeqList[0]} + 1;
						}
						elsif ($companionAssemblySeqStart->{$companionAssemblySeqList[0]} >= $filterStart && $companionAssemblySeqStart->{$companionAssemblySeqList[0]} <= $filterEnd)
						{
							$filterLength += $filterEnd - $companionAssemblySeqStart->{$companionAssemblySeqList[0]} + 1;
						}
						elsif ($companionAssemblySeqEnd->{$companionAssemblySeqList[0]} >= $filterStart && $companionAssemblySeqEnd->{$companionAssemblySeqList[0]} <= $filterEnd)
						{
							$filterLength += $companionAssemblySeqEnd->{$companionAssemblySeqList[0]} - $filterStart + 1;
						}
						else
						{
							$filterLength += $filterEnd - $filterStart + 1;
						}
					}
				}
			}
			my $ctgBarX = $margin + $companionAssemblyCtgList[6] / $pixelUnit;
			my $ctgBarY = $barY + $barHeight * ($ctgRow + 1) * 3;
			my $companionCtgLabel = $companionAssembly[2].'.'.$companionAssembly[3].'-Ctg'.$companionAssemblyCtgList[2];
			my $companionCtgLabelLength = length $companionCtgLabel;
			if (($companionAssemblyCtgList[7] + $filterLength) / $pixelUnit < $companionCtgLabelLength * $textFontWidth) #no assembly name
			{
				$companionCtgLabel = 'Ctg'.$companionAssemblyCtgList[2];
				$companionCtgLabelLength = length $companionCtgLabel;
			}
			if (($companionAssemblyCtgList[7] + $filterLength) / $pixelUnit < $companionCtgLabelLength * $textFontWidth) # no prefix
			{
				$companionCtgLabel = $companionAssemblyCtgList[2];
				$companionCtgLabelLength = length $companionCtgLabel;
			}
			$companionAssemblyChrCtg->text(
				id      => 'ctgName'.$companionAssemblyCtgList[0],
				class   => 'hasmenuForCtg',
				x       => $ctgBarX,
				y       => $ctgBarY + $barHeight - 2,
				style   => {
					'font-size'   => $textFontSize,
					'stroke'        =>  'black'
				},
				'edit-url' => "assemblyCtgEdit.cgi?assemblyCtgId=$companionAssemblyCtgList[0]&assemblyId=$assembly[0]&chr=$chr&scrollLeft=$ctgBarX",
				'flip-url' => "assemblyCtgFlip.cgi?assemblyCtgId=$companionAssemblyCtgList[0]&assemblyId=$assembly[0]&chr=$chr&scrollLeft=$ctgBarX",
				'comment-url' => "comment.cgi?itemId=$companionAssemblyCtgList[0]",
				'fill-url' => "assemblyGapFillForm.cgi?assemblyCtgId=$companionAssemblyCtgList[0]"
			)->cdata($companionCtgLabel) if (($companionAssemblyCtgList[7] + $filterLength) / $pixelUnit > $companionCtgLabelLength * $textFontWidth);
			
			$companionAssemblyChrCtg->rectangle(
				id    => $companionAssemblyCtgList[0],
				class   => 'hasmenuForCtg draggable',
				x     => $ctgBarX,
				y     => $ctgBarY,
				width => ($companionAssemblyCtgList[7] + $filterLength) / $pixelUnit,
				height=> $barHeight,
				style => { stroke => 'black',
							fill => 'blue',
							'fill-opacity' => 0.3
						},
				'edit-url' => "assemblyCtgEdit.cgi?assemblyCtgId=$companionAssemblyCtgList[0]&assemblyId=$assembly[0]&chr=$chr&scrollLeft=$ctgBarX",
				'flip-url' => "assemblyCtgFlip.cgi?assemblyCtgId=$companionAssemblyCtgList[0]&assemblyId=$assembly[0]&chr=$chr&scrollLeft=$ctgBarX",
				'fill-url' => "assemblyGapFillForm.cgi?assemblyCtgId=$companionAssemblyCtgList[0]",
				'comment-url' => "comment.cgi?itemId=$companionAssemblyCtgList[0]",
				'move-url' => "assemblyCtgMove.cgi?assemblyCtgId=$companionAssemblyCtgList[0]&assemblyId=$assembly[0]&chr=$chr&assemblyPosition="
			);


			my @companionAssemblySeq = sort { $companionAssemblySeqOrder->{$a} <=> $companionAssemblySeqOrder->{$b} } keys %$companionAssemblySeqOrder;
			$totalSeqs += @companionAssemblySeq;
			my $preSeq = '';
			my $companionAssemblySeqLeftEnd;
			my $companionAssemblySeqRightEnd;
			my $hiddenSeqCount = 0;
			my @hiddenSeqPosition;
			my $companionAssemblySeqCol;
			my $seqCount = 0;
			for my $currentSeq (@companionAssemblySeq)
			{
				my $sequenceDetails = decode_json $companionAssemblySequenceData->{$currentSeq};
				$sequenceDetails->{'id'} = '' unless (exists $sequenceDetails->{'id'});
				$sequenceDetails->{'description'} = '' unless (exists $sequenceDetails->{'description'});
				$sequenceDetails->{'gapList'} = '' unless (exists $sequenceDetails->{'gapList'});
				$sequenceDetails->{'filter'} = '' unless (exists $sequenceDetails->{'filter'});
				if ($companionAssemblySeqHide->{$currentSeq} > 0)
				{
					$companionAssemblySeqLeftEnd->{$currentSeq} = ($preSeq) ? $companionAssemblySeqLeftEnd->{$preSeq} + $barHeight * ($hiddenSeqCount + 1) * $pixelUnit : $barHeight * ($hiddenSeqCount + 1) * $pixelUnit;
					my $hiddenRow = @hiddenSeqPosition;
					my $i = 0;
					for (@hiddenSeqPosition)
					{
						if ($_ < $companionAssemblySeqLeftEnd->{$currentSeq})
						{
							$hiddenRow = $i;
							last;
						}
						$i++;
					}
					my $hiddenSeqBarX = $margin + ($companionAssemblyCtgList[6] + $companionAssemblySeqLeftEnd->{$currentSeq}) / $pixelUnit;
					my $hiddenSeqBarY = $barY + $barSpacing - $hiddenSeqPosition - $barHeight * $hiddenRow * 2;
					$companionAssemblyChrSeqHidden->rectangle(
							id    => $currentSeq,
							onclick => "closeDialog();openDialog('seqView.cgi?seqId=$companionAssemblySequenceId->{$currentSeq}')",
							class   => 'hasmenuForSeq',
							x     => $hiddenSeqBarX,
							y     => $hiddenSeqBarY,
							width => $companionAssemblySeqLength->{$currentSeq} / $pixelUnit,
							height=> $barHeight,
							style => { stroke => ($currentSeq eq $highlight) ? 'red' : 'grey',
										fill => 'lightgrey',
										opacity => 0.5
										},
							'blast-url' => "blastTwoseqForm.cgi?assemblyId=$companionAssemblyId&seqOne=$companionAssemblySequenceId->{$currentSeq}",
							'align-url' => "alignmentCheckForm.cgi?assemblyId=$companionAssemblyId&seqOne=$companionAssemblySequenceId->{$currentSeq}",
							'view-url' => "seqView.cgi?seqId=$companionAssemblySequenceId->{$currentSeq}"
						);
					$scrollLeft = $companionAssemblySeqLeftEnd->{$currentSeq} / $pixelUnit if ($currentSeq eq $highlight);
					my $hiddenCompanionSeqLabel = 'Ctg'. $companionAssemblyCtgList[2] ."-". $companionAssemblySeqOrder->{$currentSeq}. "." .$companionAssemblySequenceName->{$currentSeq}."  FPC:".$companionAssemblySeqFpcCtg->{$currentSeq};
					my $hiddenCompanionSeqLabelLength = length $hiddenCompanionSeqLabel;
					if ($companionAssemblySeqLength->{$currentSeq} / $pixelUnit < $hiddenCompanionSeqLabelLength * $textFontWidth) #no FPC
					{
						$hiddenCompanionSeqLabel = 'Ctg'. $companionAssemblyCtgList[2] ."-". $companionAssemblySeqOrder->{$currentSeq}. "." .$companionAssemblySequenceName->{$currentSeq};
						$hiddenCompanionSeqLabelLength = length $hiddenCompanionSeqLabel;
					}
					if ($companionAssemblySeqLength->{$currentSeq} / $pixelUnit < $hiddenCompanionSeqLabelLength * $textFontWidth) #no Ctg
					{
						$hiddenCompanionSeqLabel = $companionAssemblySeqOrder->{$currentSeq}. "." .$companionAssemblySequenceName->{$currentSeq};
						$hiddenCompanionSeqLabelLength = length $hiddenCompanionSeqLabel;
					}
					if ($companionAssemblySeqLength->{$currentSeq} / $pixelUnit < $hiddenCompanionSeqLabelLength * $textFontWidth) #no order
					{
						$hiddenCompanionSeqLabel = $companionAssemblySequenceName->{$currentSeq};
						$hiddenCompanionSeqLabelLength = length $hiddenCompanionSeqLabel;
					}

					my $hiddenTextY = $barY + $barSpacing - $hiddenSeqPosition - $barHeight * $hiddenRow * 2 + $barHeight - 2;
					$companionAssemblyChrSeqHidden->text(
						id      => 'seqName'.$currentSeq,
						x       => $hiddenSeqBarX,
						y       => $hiddenTextY,
						style   => {
							'font-size'   => $textFontSize,
							'stroke'        => ($sequenceDetails->{'gapList'}) ? 'pink' : 'grey'
						},
						onclick => "closeDialog();openDialog('seqView.cgi?seqId=$companionAssemblySequenceId->{$currentSeq}')",
						class   => 'hasmenuForSeq',
						'blast-url' => "blastTwoseqForm.cgi?assemblyId=$companionAssemblyId&seqOne=$companionAssemblySequenceId->{$currentSeq}",
						'align-url' => "alignmentCheckForm.cgi?assemblyId=$companionAssemblyId&seqOne=$companionAssemblySequenceId->{$currentSeq}",
						'view-url' => "seqView.cgi?seqId=$companionAssemblySequenceId->{$currentSeq}"
					)->cdata($hiddenCompanionSeqLabel) if ($companionAssemblySeqLength->{$currentSeq} / $pixelUnit > $hiddenCompanionSeqLabelLength * $textFontWidth);
					$hiddenSeqPosition[$hiddenRow] = $companionAssemblySeqLeftEnd->{$currentSeq} + $companionAssemblySeqLength->{$currentSeq} - 1;
					$hiddenSeqCount++;
					next;
				}

				if($preSeq)
				{
					my $preSeqEnd;
					if($companionAssemblySeqOrient->{$preSeq} eq "+")
					{
						$preSeqEnd = $companionAssemblySeqRightEnd->{$preSeq} - ($companionAssemblySeqLength->{$preSeq} - $companionAssemblySeqEnd->{$preSeq});
					}
					else
					{
						$preSeqEnd = $companionAssemblySeqRightEnd->{$preSeq} - $companionAssemblySeqStart->{$preSeq} + 1;
					}
					if($companionAssemblySeqOrient->{$currentSeq} eq "+")
					{
						$companionAssemblySeqLeftEnd->{$currentSeq} = $preSeqEnd + 1 - $companionAssemblySeqStart->{$currentSeq} + 1;
						$companionAssemblySeqRightEnd->{$currentSeq} = $preSeqEnd + 1 - $companionAssemblySeqStart->{$currentSeq} + 1 + $companionAssemblySeqLength->{$currentSeq} - 1;
					}
					else
					{
						$companionAssemblySeqLeftEnd->{$currentSeq} = $preSeqEnd + 1 - ($companionAssemblySeqLength->{$currentSeq} - $companionAssemblySeqEnd->{$currentSeq});
						$companionAssemblySeqRightEnd->{$currentSeq} = $preSeqEnd + 1 - ($companionAssemblySeqLength->{$currentSeq} - $companionAssemblySeqEnd->{$currentSeq}) + $companionAssemblySeqLength->{$currentSeq} - 1;
					}
				}
				else
				{
					$companionAssemblySeqLeftEnd->{$currentSeq} = 1;
					$companionAssemblySeqRightEnd->{$currentSeq} = $companionAssemblySeqLength->{$currentSeq};
				}
				my $goodCol = $seqCount % 2;
				$companionAssemblySeqCol->{$currentSeq} = $goodCol;

				my $goodSequenceX = ($companionAssemblySeqOrient->{$currentSeq} eq "+") ? $margin + ($companionAssemblyCtgList[6] + $companionAssemblySeqLeftEnd->{$currentSeq} + $companionAssemblySeqStart->{$currentSeq} - 1) / $pixelUnit : $margin + ($companionAssemblyCtgList[6] + $companionAssemblySeqRightEnd->{$currentSeq} - $companionAssemblySeqEnd->{$currentSeq} + 1) / $pixelUnit ;
				$companionAssemblyChrSeq->rectangle(
						x     => $goodSequenceX,
						y     => $ctgBarY - $barHeight * (1 + $goodCol),
						width => ($companionAssemblySeqEnd->{$currentSeq} - $companionAssemblySeqStart->{$currentSeq} + 1) / $pixelUnit,
						height=> $barHeight,
						style => { stroke => 'lightgreen',
									fill => 'lightgreen'
							},
						id    => "Sequence$currentSeq$$"
					);
				if ($sequenceDetails->{'gapList'}) # draw gap
				{
					foreach (split ",", $sequenceDetails->{'gapList'} )
					{
						my ($gapStart,$gapEnd) = split "-", $_;
						my $gapX = ($companionAssemblySeqOrient->{$currentSeq} eq "+") ? $margin + ($companionAssemblyCtgList[6] + $companionAssemblySeqLeftEnd->{$currentSeq} + $gapStart - 1) / $pixelUnit : $margin + ($companionAssemblyCtgList[6] + $companionAssemblySeqRightEnd->{$currentSeq} - $gapEnd + 1) / $pixelUnit ;
						$companionAssemblyChrSeq->rectangle(
							x     => $gapX,
							y     => $ctgBarY - $barHeight * (1 + $goodCol),
							width => ($gapEnd - $gapStart + 1) / $pixelUnit,
							height=> $barHeight,
							style => { stroke => 'red',
										fill => 'red'
								},
							id    => "Gap$currentSeq$gapStart$gapEnd$$"
						);
					}
				}
				if ($sequenceDetails->{'filter'}) # draw filter
				{
					foreach (split ",", $sequenceDetails->{'filter'} )
					{
						my ($filterStart,$filterEnd) = split "-", $_;
						my $filterX = ($companionAssemblySeqOrient->{$currentSeq} eq "+") ? $margin + ($companionAssemblyCtgList[6] + $companionAssemblySeqLeftEnd->{$currentSeq} + $filterStart - 1) / $pixelUnit : $margin + ($companionAssemblyCtgList[6] + $companionAssemblySeqRightEnd->{$currentSeq} - $filterEnd + 1) / $pixelUnit ;
						$companionAssemblyChrSeq->rectangle(
							x     => $filterX,
							y     => $ctgBarY - $barHeight * (1 + $goodCol),
							width => ($filterEnd - $filterStart + 1) / $pixelUnit,
							height=> $barHeight,
							style => { stroke => 'black',
										fill => 'lightgrey'
								},
							id    => "Filter$currentSeq$filterStart$filterEnd$$"
						);
					}
				}

				my $seqBarX = $margin + ($companionAssemblyCtgList[6] + $companionAssemblySeqLeftEnd->{$currentSeq}) / $pixelUnit;
				$companionAssemblyChrSeq->rectangle(
						id    => $currentSeq,
						onclick => "closeDialog();openDialog('seqView.cgi?seqId=$companionAssemblySequenceId->{$currentSeq}')",
						class   => 'hasmenuForSeq',
						x     => $seqBarX,
						y     => $ctgBarY - $barHeight * (1 + $goodCol),
						width => $companionAssemblySeqLength->{$currentSeq} / $pixelUnit,
						height=> $barHeight,
						style => { stroke => ($currentSeq eq $highlight) ? 'red' : 'black',
									fill => 'lightgrey',
									'fill-opacity' =>  0.5
							},
						'blast-url' => "blastTwoseqForm.cgi?assemblyId=$companionAssemblyId&seqOne=$companionAssemblySequenceId->{$currentSeq}",
						'align-url' => "alignmentCheckForm.cgi?assemblyId=$companionAssemblyId&seqOne=$companionAssemblySequenceId->{$currentSeq}",
						'view-url' => "seqView.cgi?seqId=$companionAssemblySequenceId->{$currentSeq}"
					);
				$scrollLeft = $companionAssemblySeqLeftEnd->{$currentSeq} / $pixelUnit if ($currentSeq eq $highlight);
				my $companionSeqLabel = $companionAssemblySeqOrder->{$currentSeq}. "." . $companionAssemblySequenceName->{$currentSeq} . "(" .$companionAssemblySeqOrient->{$currentSeq}.")  FPC:".$companionAssemblySeqFpcCtg->{$currentSeq};
				my $companionSeqLabelLength = length $companionSeqLabel;
				if ($companionAssemblySeqLength->{$currentSeq} / $pixelUnit < $companionSeqLabelLength * $textFontWidth) # no FPC
				{
					$companionSeqLabel = $companionAssemblySeqOrder->{$currentSeq}. "." . $companionAssemblySequenceName->{$currentSeq} . "(" .$companionAssemblySeqOrient->{$currentSeq}.")";
					$companionSeqLabelLength = length $companionSeqLabel;
				}
				if ($companionAssemblySeqLength->{$currentSeq} / $pixelUnit < $companionSeqLabelLength * $textFontWidth) # no orientation
				{
					$companionSeqLabel = $companionAssemblySeqOrder->{$currentSeq}. "." . $companionAssemblySequenceName->{$currentSeq};
					$companionSeqLabelLength = length $companionSeqLabel;
				}
				if ($companionAssemblySeqLength->{$currentSeq} / $pixelUnit < $companionSeqLabelLength * $textFontWidth) # no order
				{
					$companionSeqLabel = $companionAssemblySequenceName->{$currentSeq};
					$companionSeqLabelLength = length $companionSeqLabel;
				}

				$companionAssemblyChrSeq->text(
					id      => 'seqName'.$currentSeq,
					x       => $seqBarX,
					y       => $ctgBarY - $barHeight * $goodCol - 2,
					style   => {
						'font-size'   => $textFontSize,
						'stroke'        => ($sequenceDetails->{'gapList'}) ? 'red' : 'black'
					},
					onclick => "closeDialog();openDialog('seqView.cgi?seqId=$companionAssemblySequenceId->{$currentSeq}')",
					class   => 'hasmenuForSeq',
					'blast-url' => "blastTwoseqForm.cgi?assemblyId=$companionAssemblyId&seqOne=$companionAssemblySequenceId->{$currentSeq}",
					'align-url' => "alignmentCheckForm.cgi?assemblyId=$companionAssemblyId&seqOne=$companionAssemblySequenceId->{$currentSeq}",
					'view-url' => "seqView.cgi?seqName=$currentSeq"
				)->cdata($companionSeqLabel) if ($companionAssemblySeqLength->{$currentSeq} / $pixelUnit > $companionSeqLabelLength * $textFontWidth);

				my $queryDir;
				for (my $position = 0; $position < length($companionAssemblySequenceId->{$currentSeq}); $position += 2)
				{
					$queryDir .= "/q". substr($companionAssemblySequenceId->{$currentSeq},$position,2);
				}
				my $subjectDir;
				for (my $position = 0; $position < length($refGenomeSequence[0]); $position += 2)
				{
					$subjectDir .= "/s". substr($refGenomeSequence[0],$position,2);
				}
				my $alignmentCount = 0;
				open (TBL,"$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$companionAssemblySequenceId->{$currentSeq}-$refGenomeSequence[0].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$companionAssemblySequenceId->{$currentSeq}-$refGenomeSequence[0].tbl";
				while(<TBL>)
				{
					chop;
					/^#/ and next;
					$alignmentCount++;
					my @alignments = split("\t",$_);
					next if ($alignments[12] > 0);
					next if ($alignments[3] < $alignmentLength); #skip if alignment shorter than alignmentLength
					my $seqXOne = 0;
					my $seqXTwo = 0;
					if($companionAssemblySeqOrient->{$currentSeq} eq "+")
					{
						$seqXOne = $margin + ($companionAssemblyCtgList[6] + $companionAssemblySeqLeftEnd->{$currentSeq} + $alignments[6] - 1) / $pixelUnit;
						$seqXTwo = $margin + ($companionAssemblyCtgList[6] + $companionAssemblySeqLeftEnd->{$currentSeq} + $alignments[7] - 1) / $pixelUnit;
					}
					else
					{
						$seqXOne = $margin + ($companionAssemblyCtgList[6] + $companionAssemblySeqRightEnd->{$currentSeq} - $alignments[6] + 1) / $pixelUnit;
						$seqXTwo = $margin + ($companionAssemblyCtgList[6] + $companionAssemblySeqRightEnd->{$currentSeq} - $alignments[7] + 1) / $pixelUnit;
					}
					my $refXOne = $margin + ($alignments[8] - 1 ) / $pixelUnit;;
					my $refXTwo = $margin + ($alignments[9] - 1 ) / $pixelUnit;

					my $xv = [$seqXOne,$seqXTwo,$refXTwo,$refXOne];
					my $yv = [$ctgBarY + $barHeight + 1,$ctgBarY + $barHeight + 1,$barY + $spaceForCompanion - 1,$barY + $spaceForCompanion - 1];
					my $points = $companionAssemblyChrSeqAlignment->get_path(
						x=>$xv,
						y=>$yv,
						-type=>'polygon'
					);
					$companionAssemblyChrSeqAlignment->polygon(
						%$points,
						id=>'companionAln'.$alignments[0].$alignments[1].$alignmentCount,
						onclick => "closeDialog();openDialog('alignmentView.cgi?alignmentId=$alignments[0]-$alignments[1]-$alignmentCount')",
						class=>'hasmenuForAlignment',
						style=>{ stroke => 'red',
							fill => 'yellow',
							opacity => 0.5
							},
						'hide-url' => "assemblyAlignmentHide.cgi?alignmentId=$alignments[0]-$alignments[1]-$alignmentCount&assemblyId=$companionAssemblyId&openAssemblyId=$assemblyId&chr=$chr&scrollLeft=$seqXOne",
						'view-url' => "alignmentView.cgi?alignmentId=$alignments[0]-$alignments[1]-$alignmentCount"
					);
				}
				close(TBL);
				$preSeq = $currentSeq;
				$hiddenSeqCount = 0;
				$seqCount++;
			}
		}
	}

	my @ctgPosition;
	my $assemblyCtgList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblyCtg' AND o = ? AND x = ? ORDER BY y");
	$assemblyCtgList->execute($assemblyId,$chr);
	while (my @assemblyCtgList = $assemblyCtgList->fetchrow_array())
	{
		$totalLength->{'assembly'} = 0 unless (exists $totalLength->{'assembly'});
		$totalLength->{'assembly'} += $assemblyCtgList[7];
		my $ctgRow = @ctgPosition;
		my $j = 0;
		for (@ctgPosition)
		{
			if ($_ < $margin + $assemblyCtgList[6] / $pixelUnit)
			{
				$ctgRow = $j;
				last;
			}
			$j++;
		}
		$ctgPosition[$ctgRow] = $margin + ($assemblyCtgList[6] + $assemblyCtgList[7]) / $pixelUnit;

		my $assemblySequenceName;
		my $assemblySequenceId;
		my $assemblySequenceData;
		my $assemblySequenceType;
		my $assemblySeqLength;
		my $assemblySeqOrient;
		my $assemblySeqStart;
		my $assemblySeqEnd;
		my $assemblySeqOrder;
		my $assemblySeqHide;
		my $assemblySeqFpcCtg;
		my $order = 0;
		my $filterLength = 0;
		foreach (split ",", $assemblyCtgList[8])
		{
			next unless ($_);
			my $hide = ($_ =~ /^-/) ? 1 : 0;
			$_ =~ s/[^a-zA-Z0-9]//g;
			$assemblySeqHide->{$_} = $hide;
			$order++;
			$assemblySeqOrder->{$_} = $order;
			my $assemblySeqList=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
			$assemblySeqList->execute($_);
			my @assemblySeqList = $assemblySeqList->fetchrow_array();
			$assemblySequenceName->{$assemblySeqList[0]} = $assemblySeqList[2];
			$assemblySequenceId->{$assemblySeqList[0]} = $assemblySeqList[5];
			$assemblySeqLength->{$assemblySeqList[0]} = $assemblySeqList[6];
			$assemblySeqOrient->{$assemblySeqList[0]} = ($assemblySeqList[7] < 0) ? "-" : "+";
			if($assemblySeqList[8])
			{
				($assemblySeqStart->{$assemblySeqList[0]},$assemblySeqEnd->{$assemblySeqList[0]}) = split ",",$assemblySeqList[8];
			}
			else
			{
				$assemblySeqStart->{$assemblySeqList[0]} = 1;
				$assemblySeqEnd->{$assemblySeqList[0]} = $assemblySeqList[6];
				my $updateAssemblySeq=$dbh->do("UPDATE matrix SET note = '$assemblySeqStart->{$assemblySeqList[0]},$assemblySeqEnd->{$assemblySeqList[0]}' WHERE id = $assemblySeqList[0]");
			}
			my $assemblySequence=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
			$assemblySequence->execute($assemblySeqList[5]);
			my @assemblySequence = $assemblySequence->fetchrow_array();
			$assemblySequenceType->{$assemblySeqList[0]} = $assemblySequence[3];
			$assemblySequenceData->{$assemblySeqList[0]} = $assemblySequence[8];
			$assemblySeqFpcCtg->{$assemblySeqList[0]} = "NA";
			my $fpcCloneList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'fpcClone' AND o = ? AND name LIKE ?");
			$fpcCloneList->execute($assembly[6],$assemblySeqList[2]);
			while (my @fpcCloneList = $fpcCloneList->fetchrow_array())
			{
				$fpcCloneList[8] =~ /Map "(.*)"/;
				$assemblySeqFpcCtg->{$assemblySeqList[0]} = $1;
			}
			my $assemblySequenceDetails = decode_json $assemblySequence[8];
			$assemblySequenceDetails->{'filter'} = '' unless (exists $assemblySequenceDetails->{'filter'});
			if ($assemblySequenceDetails->{'filter'}) 
			{
				foreach (split ",", $assemblySequenceDetails->{'filter'} )
				{
					my ($filterStart,$filterEnd) = split "-", $_;
					next if ($assemblySeqStart->{$assemblySeqList[0]} > $filterEnd);
					next if ($assemblySeqEnd->{$assemblySeqList[0]} < $filterStart);
					if ($assemblySeqStart->{$assemblySeqList[0]} >= $filterStart && $assemblySeqEnd->{$assemblySeqList[0]} <= $filterEnd)
					{
						$filterLength += $assemblySeqEnd->{$assemblySeqList[0]} - $assemblySeqStart->{$assemblySeqList[0]} + 1;
					}
					elsif ($assemblySeqStart->{$assemblySeqList[0]} >= $filterStart && $assemblySeqStart->{$assemblySeqList[0]} <= $filterEnd)
					{
						$filterLength += $filterEnd - $assemblySeqStart->{$assemblySeqList[0]} + 1;
					}
					elsif ($assemblySeqEnd->{$assemblySeqList[0]} >= $filterStart && $assemblySeqEnd->{$assemblySeqList[0]} <= $filterEnd)
					{
						$filterLength += $assemblySeqEnd->{$assemblySeqList[0]} - $filterStart + 1;
					}
					else
					{
						$filterLength += $filterEnd - $filterStart + 1;
					}
				}
			}
		}
		my $ctgBarX = $margin + $assemblyCtgList[6] / $pixelUnit;
		my $ctgBarY = $barY + $barSpacing - $barHeight * $ctgRow * 3 + $spaceForCompanion;

		my $ctgLabel = $assembly[2].'.'.$assembly[3].'-Ctg'.$assemblyCtgList[2];
		my $ctgLabelLength = length $ctgLabel;
		if (($assemblyCtgList[7] + $filterLength) / $pixelUnit < $ctgLabelLength * $textFontWidth) #no assembly name
		{
			$ctgLabel = 'Ctg'.$assemblyCtgList[2];
			$ctgLabelLength = length $ctgLabel;
		}
		if (($assemblyCtgList[7] + $filterLength) / $pixelUnit < $ctgLabelLength * $textFontWidth) #no prefix
		{
			$ctgLabel = $assemblyCtgList[2];
			$ctgLabelLength = length $ctgLabel;
		}
		$assemblyChrCtg->text(
			id      => 'ctgName'.$assemblyCtgList[0],
			class   => 'hasmenuForCtg',
			x       => $ctgBarX,
			y       => $ctgBarY + $barHeight - 2,
			style   => {
				'font-size'   => $textFontSize,
				'stroke'        =>  'black'
			},
			'edit-url' => "assemblyCtgEdit.cgi?assemblyCtgId=$assemblyCtgList[0]&chr=$chr&scrollLeft=$ctgBarX",
			'flip-url' => "assemblyCtgFlip.cgi?assemblyCtgId=$assemblyCtgList[0]&chr=$chr&scrollLeft=$ctgBarX",
			'comment-url' => "comment.cgi?itemId=$assemblyCtgList[0]",
			'fill-url' => "assemblyGapFillForm.cgi?assemblyCtgId=$assemblyCtgList[0]"
		)->cdata($ctgLabel) if (($assemblyCtgList[7] + $filterLength) / $pixelUnit > $ctgLabelLength * $textFontWidth);
		
		$assemblyChrCtg->rectangle(
			id    => $assemblyCtgList[0],
			class   => 'hasmenuForCtg draggable',
			x     => $ctgBarX,
			y     => $ctgBarY,
			width => ($assemblyCtgList[7] + $filterLength) / $pixelUnit,
			height=> $barHeight,
			style => { stroke => 'black',
						fill => 'blue',
						'fill-opacity' => 0.3
					},
			'edit-url' => "assemblyCtgEdit.cgi?assemblyCtgId=$assemblyCtgList[0]&chr=$chr&scrollLeft=$ctgBarX",
			'flip-url' => "assemblyCtgFlip.cgi?assemblyCtgId=$assemblyCtgList[0]&chr=$chr&scrollLeft=$ctgBarX",
			'fill-url' => "assemblyGapFillForm.cgi?assemblyCtgId=$assemblyCtgList[0]",
			'comment-url' => "comment.cgi?itemId=$assemblyCtgList[0]",
			'move-url' => "assemblyCtgMove.cgi?assemblyCtgId=$assemblyCtgList[0]&chr=$chr&assemblyPosition="
		);

		my @assemblySeq = sort { $assemblySeqOrder->{$a} <=> $assemblySeqOrder->{$b} } keys %$assemblySeqOrder;
		$totalSeqs += @assemblySeq;
		my $preSeq = '';
		my $assemblySeqLeftEnd;
		my $assemblySeqRightEnd;
		my $hiddenSeqCount = 0;
		my @hiddenSeqPosition;
		my $assemblySeqCol;
		my $seqCount = 0;

		for my $currentSeq (@assemblySeq)
		{
			my $sequenceDetails = decode_json $assemblySequenceData->{$currentSeq};
			$sequenceDetails->{'id'} = '' unless (exists $sequenceDetails->{'id'});
			$sequenceDetails->{'description'} = '' unless (exists $sequenceDetails->{'description'});
			$sequenceDetails->{'gapList'} = '' unless (exists $sequenceDetails->{'gapList'});
			$sequenceDetails->{'filter'} = '' unless (exists $sequenceDetails->{'filter'});
			if ($assemblySeqHide->{$currentSeq} > 0)
			{
				$assemblySeqLeftEnd->{$currentSeq} = ($preSeq) ? $assemblySeqLeftEnd->{$preSeq} + $barHeight * ($hiddenSeqCount + 1) * $pixelUnit : $barHeight * ($hiddenSeqCount + 1) * $pixelUnit;
				my $hiddenRow = @hiddenSeqPosition;
				my $i = 0;
				for (@hiddenSeqPosition)
				{
					if ($_ < $assemblySeqLeftEnd->{$currentSeq})
					{
						$hiddenRow = $i;
						last;
					}
					$i++;
				}
				my $hiddenSeqBarX = $margin + ($assemblyCtgList[6] + $assemblySeqLeftEnd->{$currentSeq}) / $pixelUnit;
				my $hiddenSeqBarY = $barY + $hiddenSeqPosition + $barHeight * $hiddenRow * 2 + $spaceForCompanion;
				$assemblyChrSeqHidden->rectangle(
						id    => $currentSeq,
						onclick => "closeDialog();openDialog('seqView.cgi?seqId=$assemblySequenceId->{$currentSeq}')",
						class   => 'hasmenuForSeq',
						x     => $hiddenSeqBarX,
						y     => $hiddenSeqBarY,
						width => $assemblySeqLength->{$currentSeq} / $pixelUnit,
						height=> $barHeight,
						style => { stroke => ($currentSeq eq $highlight) ? 'red' : 'grey',
									fill => 'lightgrey',
									opacity => 0.5
									},
						'blast-url' => "blastTwoseqForm.cgi?assemblyId=$assemblyId&seqOne=$assemblySequenceId->{$currentSeq}",
						'align-url' => "alignmentCheckForm.cgi?assemblyId=$assemblyId&seqOne=$assemblySequenceId->{$currentSeq}",
						'view-url' => "seqView.cgi?seqId=$assemblySequenceId->{$currentSeq}"
					);
				$scrollLeft = $assemblySeqLeftEnd->{$currentSeq} / $pixelUnit if ($currentSeq eq $highlight);
				my $hiddenSeqLabel = 'Ctg'. $assemblyCtgList[2] ."-". $assemblySeqOrder->{$currentSeq}. "." .$assemblySequenceName->{$currentSeq}."  FPC:".$assemblySeqFpcCtg->{$currentSeq};
				my $hiddenSeqLabelLength = length $hiddenSeqLabel;
				if ($assemblySeqLength->{$currentSeq} / $pixelUnit < $hiddenSeqLabelLength * $textFontWidth) #no FPC
				{
					$hiddenSeqLabel = 'Ctg'. $assemblyCtgList[2] ."-". $assemblySeqOrder->{$currentSeq}. "." .$assemblySequenceName->{$currentSeq};
					$hiddenSeqLabelLength = length $hiddenSeqLabel;
				}
				if ($assemblySeqLength->{$currentSeq} / $pixelUnit < $hiddenSeqLabelLength * $textFontWidth) #no Ctg
				{
					$hiddenSeqLabel = $assemblySeqOrder->{$currentSeq}. "." .$assemblySequenceName->{$currentSeq};
					$hiddenSeqLabelLength = length $hiddenSeqLabel;
				}
				if ($assemblySeqLength->{$currentSeq} / $pixelUnit < $hiddenSeqLabelLength * $textFontWidth) #no order
				{
					$hiddenSeqLabel = $assemblySequenceName->{$currentSeq};
					$hiddenSeqLabelLength = length $hiddenSeqLabel;
				}

				my $hiddenTextY = $barY + $hiddenSeqPosition + $barHeight * $hiddenRow * 2 + $barHeight - 2 + $spaceForCompanion;
				$assemblyChrSeqHidden->text(
					id      => 'seqName'.$currentSeq,
					x       => $hiddenSeqBarX,
					y       => $hiddenTextY,
					style   => {
						'font-size'   => $textFontSize,
						'stroke'        => ($sequenceDetails->{'gapList'}) ? 'pink' : 'grey'
					},
					onclick => "closeDialog();openDialog('seqView.cgi?seqId=$assemblySequenceId->{$currentSeq}')",
					class   => 'hasmenuForSeq',
					'blast-url' => "blastTwoseqForm.cgi?assemblyId=$assemblyId&seqOne=$assemblySequenceId->{$currentSeq}",
					'align-url' => "alignmentCheckForm.cgi?assemblyId=$assemblyId&seqOne=$assemblySequenceId->{$currentSeq}",
					'view-url' => "seqView.cgi?seqId=$assemblySequenceId->{$currentSeq}"
				)->cdata($hiddenSeqLabel) if ($assemblySeqLength->{$currentSeq} / $pixelUnit > $hiddenSeqLabelLength * $textFontWidth);
				$hiddenSeqPosition[$hiddenRow] = $assemblySeqLeftEnd->{$currentSeq} + $assemblySeqLength->{$currentSeq} - 1;
				$hiddenSeqCount++;
				next;
			}

			if($preSeq)
			{
				my $preSeqEnd;
				if($assemblySeqOrient->{$preSeq} eq "+")
				{
					$preSeqEnd = $assemblySeqRightEnd->{$preSeq} - ($assemblySeqLength->{$preSeq} - $assemblySeqEnd->{$preSeq});
				}
				else
				{
					$preSeqEnd = $assemblySeqRightEnd->{$preSeq} - $assemblySeqStart->{$preSeq} + 1;
				}
				if($assemblySeqOrient->{$currentSeq} eq "+")
				{
					$assemblySeqLeftEnd->{$currentSeq} = $preSeqEnd + 1 - $assemblySeqStart->{$currentSeq} + 1;
					$assemblySeqRightEnd->{$currentSeq} = $preSeqEnd + 1 - $assemblySeqStart->{$currentSeq} + 1 + $assemblySeqLength->{$currentSeq} - 1;
				}
				else
				{
					$assemblySeqLeftEnd->{$currentSeq} = $preSeqEnd + 1 - ($assemblySeqLength->{$currentSeq} - $assemblySeqEnd->{$currentSeq});
					$assemblySeqRightEnd->{$currentSeq} = $preSeqEnd + 1 - ($assemblySeqLength->{$currentSeq} - $assemblySeqEnd->{$currentSeq}) + $assemblySeqLength->{$currentSeq} - 1;
				}
			}
			else
			{
				$assemblySeqLeftEnd->{$currentSeq} = 1;
				$assemblySeqRightEnd->{$currentSeq} = $assemblySeqLength->{$currentSeq};
			}
			my $goodCol = $seqCount % 2;
			$assemblySeqCol->{$currentSeq} = $goodCol;

			my $goodSequenceX = ($assemblySeqOrient->{$currentSeq} eq "+") ? $margin + ($assemblyCtgList[6] + $assemblySeqLeftEnd->{$currentSeq} + $assemblySeqStart->{$currentSeq} - 1) / $pixelUnit : $margin + ($assemblyCtgList[6] + $assemblySeqRightEnd->{$currentSeq} - $assemblySeqEnd->{$currentSeq} + 1) / $pixelUnit ;
			$assemblyChrSeq->rectangle(
					x     => $goodSequenceX,
					y     => $ctgBarY + $barHeight * (1 + $goodCol),
					width => ($assemblySeqEnd->{$currentSeq} - $assemblySeqStart->{$currentSeq} + 1) / $pixelUnit,
					height=> $barHeight,
					style => { stroke => 'lightgreen',
								fill => 'lightgreen'
								},
					id    => "Sequence$currentSeq$$"
				);
			if ($sequenceDetails->{'gapList'}) # draw gap
			{
				foreach (split ",", $sequenceDetails->{'gapList'} )
				{
					my ($gapStart,$gapEnd) = split "-", $_;
					my $gapX = ($assemblySeqOrient->{$currentSeq} eq "+") ? $margin + ($assemblyCtgList[6] + $assemblySeqLeftEnd->{$currentSeq} + $gapStart - 1) / $pixelUnit : $margin + ($assemblyCtgList[6] + $assemblySeqRightEnd->{$currentSeq} - $gapEnd + 1) / $pixelUnit ;
					$assemblyChrSeq->rectangle(
						x     => $gapX,
						y     => $ctgBarY + $barHeight * (1 + $goodCol),
						width => ($gapEnd - $gapStart + 1) / $pixelUnit,
						height=> $barHeight,
						style => { stroke => 'red',
									fill => 'red'
									},
						id    => "Gap$currentSeq$gapStart$gapEnd$$"
					);
				}
			}
			if ($sequenceDetails->{'filter'}) # draw filter
			{
				foreach (split ",", $sequenceDetails->{'filter'} )
				{
					my ($filterStart,$filterEnd) = split "-", $_;
					my $filterX = ($assemblySeqOrient->{$currentSeq} eq "+") ? $margin + ($assemblyCtgList[6] + $assemblySeqLeftEnd->{$currentSeq} + $filterStart - 1) / $pixelUnit : $margin + ($assemblyCtgList[6] + $assemblySeqRightEnd->{$currentSeq} - $filterEnd + 1) / $pixelUnit ;
					$assemblyChrSeq->rectangle(
						x     => $filterX,
						y     => $ctgBarY + $barHeight * (1 + $goodCol),
						width => ($filterEnd - $filterStart + 1) / $pixelUnit,
						height=> $barHeight,
						style => { stroke => 'black',
									fill => 'lightgrey'
									},
						id    => "Filter$currentSeq$filterStart$filterEnd$$"
					);
				}
			}

			my $seqBarX = $margin + ($assemblyCtgList[6] + $assemblySeqLeftEnd->{$currentSeq}) / $pixelUnit;
			$assemblyChrSeq->rectangle(
					id    => $currentSeq,
					onclick => "closeDialog();openDialog('seqView.cgi?seqId=$assemblySequenceId->{$currentSeq}')",
					class   => 'hasmenuForSeq',
					x     => $seqBarX,
					y     => $ctgBarY + $barHeight * (1 + $goodCol),
					width => $assemblySeqLength->{$currentSeq} / $pixelUnit,
					height=> $barHeight,
					style => { stroke => ($currentSeq eq $highlight) ? 'red' : 'black',
								fill => 'lightgrey',
								'fill-opacity' =>  0.5
								},
					'blast-url' => "blastTwoseqForm.cgi?assemblyId=$assemblyId&seqOne=$assemblySequenceId->{$currentSeq}",
					'align-url' => "alignmentCheckForm.cgi?assemblyId=$assemblyId&seqOne=$assemblySequenceId->{$currentSeq}",
					'view-url' => "seqView.cgi?seqId=$assemblySequenceId->{$currentSeq}"
				);
			$scrollLeft = $assemblySeqLeftEnd->{$currentSeq} / $pixelUnit if ($currentSeq eq $highlight);
			my $seqLabel = $assemblySeqOrder->{$currentSeq}. "." . $assemblySequenceName->{$currentSeq} . "(" .$assemblySeqOrient->{$currentSeq}.")  FPC:".$assemblySeqFpcCtg->{$currentSeq};
			my $seqLabelLength = length $seqLabel;
			if ($assemblySeqLength->{$currentSeq} / $pixelUnit < $seqLabelLength * $textFontWidth) #no FPC
			{
				$seqLabel = $assemblySeqOrder->{$currentSeq}. "." . $assemblySequenceName->{$currentSeq} . "(" .$assemblySeqOrient->{$currentSeq}.")";
				$seqLabelLength = length $seqLabel;
			}
			if ($assemblySeqLength->{$currentSeq} / $pixelUnit < $seqLabelLength * $textFontWidth) #no orientation
			{
				$seqLabel = $assemblySeqOrder->{$currentSeq}. "." . $assemblySequenceName->{$currentSeq};
				$seqLabelLength = length $seqLabel;
			}
			if ($assemblySeqLength->{$currentSeq} / $pixelUnit < $seqLabelLength * $textFontWidth) #no order
			{
				$seqLabel = $assemblySequenceName->{$currentSeq};
				$seqLabelLength = length $seqLabel;
			}
			$assemblyChrSeq->text(
				id      => 'seqName'.$currentSeq,
				x       => $seqBarX,
				y       => $ctgBarY + $barHeight * (2 + $goodCol) - 2,
				style   => {
					'font-size'   => $textFontSize,
					'stroke'        => ($sequenceDetails->{'gapList'}) ? 'red' : 'black'
				},
				onclick => "closeDialog();openDialog('seqView.cgi?seqId=$assemblySequenceId->{$currentSeq}')",
				class   => 'hasmenuForSeq',
				'blast-url' => "blastTwoseqForm.cgi?assemblyId=$assemblyId&seqOne=$assemblySequenceId->{$currentSeq}",
				'align-url' => "alignmentCheckForm.cgi?assemblyId=$assemblyId&seqOne=$assemblySequenceId->{$currentSeq}",
				'view-url' => "seqView.cgi?seqId=$assemblySequenceId->{$currentSeq}"
			)->cdata($seqLabel) if ($assemblySeqLength->{$currentSeq} / $pixelUnit > $seqLabelLength * $textFontWidth);


			my $queryDir;
			for (my $position = 0; $position < length($assemblySequenceId->{$currentSeq}); $position += 2)
			{
				$queryDir .= "/q". substr($assemblySequenceId->{$currentSeq},$position,2);
			}
			my $subjectDir;
			for (my $position = 0; $position < length($refGenomeSequence[0]); $position += 2)
			{
				$subjectDir .= "/s". substr($refGenomeSequence[0],$position,2);
			}
			my $alignmentCount = 0;
			open (TBL,"$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$assemblySequenceId->{$currentSeq}-$refGenomeSequence[0].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$assemblySequenceId->{$currentSeq}-$refGenomeSequence[0].tbl";
			while(<TBL>)
			{
				chop;
				/^#/ and next;
				$alignmentCount++;
				my @alignments = split("\t",$_);
				next if ($alignments[12] > 0);
				next if ($alignments[3] < $alignmentLength); #skip if alignment shorter than alignmentLength
				my $seqXOne = 0;
				my $seqXTwo = 0;
				if($assemblySeqOrient->{$currentSeq} eq "+")
				{
					$seqXOne = $margin + ($assemblyCtgList[6] + $assemblySeqLeftEnd->{$currentSeq} + $alignments[6] - 1) / $pixelUnit;
					$seqXTwo = $margin + ($assemblyCtgList[6] + $assemblySeqLeftEnd->{$currentSeq} + $alignments[7] - 1) / $pixelUnit;
				}
				else
				{
					$seqXOne = $margin + ($assemblyCtgList[6] + $assemblySeqRightEnd->{$currentSeq} - $alignments[6] + 1) / $pixelUnit;
					$seqXTwo = $margin + ($assemblyCtgList[6] + $assemblySeqRightEnd->{$currentSeq} - $alignments[7] + 1) / $pixelUnit;
				}
				my $refXOne = $margin + ($alignments[8] - 1 ) / $pixelUnit;;
				my $refXTwo = $margin + ($alignments[9] - 1 ) / $pixelUnit;

				my $xv = [$seqXOne,$seqXTwo,$refXTwo,$refXOne];
				my $yv = [$ctgBarY - 1,$ctgBarY - 1,$barY + $barHeight + $spaceForCompanion + 1,$barY + $barHeight + $spaceForCompanion + 1];
				my $points = $assemblyChrSeqAlignment->get_path(
					x=>$xv,
					y=>$yv,
					-type=>'polygon'
				);
				$assemblyChrSeqAlignment->polygon(
					%$points,
					id=>'aln'.$alignments[0].$alignments[1].$alignmentCount,
					onclick => "closeDialog();openDialog('alignmentView.cgi?alignmentId=$alignments[0]-$alignments[1]-$alignmentCount')",
					class=>'hasmenuForAlignment',
					style=>{ stroke => 'red',
						fill => 'yellow',
						opacity => 0.5
						},
					'hide-url' => "assemblyAlignmentHide.cgi?alignmentId=$alignments[0]-$alignments[1]-$alignmentCount&assemblyId=$assemblyId&chr=$chr&scrollLeft=$seqXOne",
					'view-url' => "alignmentView.cgi?alignmentId=$alignments[0]-$alignments[1]-$alignmentCount"
				);
			}
			close(TBL);
 			$preSeq = $currentSeq;
			$hiddenSeqCount = 0;
			$seqCount++;
		}
	}
	my $longest = $refGenomeSequence[5];
	for (keys %$totalLength)
	{
		$longest = $totalLength->{$_} if ($totalLength->{$_} > $longest);
	}

	for (my $rulerNumber = 0;$rulerNumber <= $longest;$rulerNumber += $tickLength)
	{
		#dash lines
		$ruler->line(
			id    => "rulerDash$rulerNumber",
			style => {
				stroke=> ($rulerNumber % (5*$tickLength) == 0) ? 'grey' : 'lightgrey',
				'stroke-dasharray' => '3,3'
				},
			x1    => $margin + $rulerNumber / $pixelUnit,
			y1    => $rulerY,
			x2    => $margin + $rulerNumber / $pixelUnit,
			y2    => $svgHeight
		);
		#region
		$ruler->line(
			id    => "rulerRegion$rulerNumber",
			x1    => $margin + $rulerNumber / $pixelUnit,
			y1    => $rulerY,
			x2    => $margin +  ($rulerNumber + $tickLength ) / $pixelUnit,
			y2    => $rulerY
		);
		#ticks
		$ruler->line(
			id    => "rulerTick$rulerNumber",
			x1    => $margin + $rulerNumber / $pixelUnit,
			y1    => ($rulerNumber % (5*$tickLength) == 0) ? $rulerY - 5 : $rulerY - 3,
			x2    => $margin + $rulerNumber / $pixelUnit,
			y2    => $rulerY
			
		);
		if($rulerNumber % (5*$tickLength) == 0)
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
	my $graphic = $svg->xmlify;
	my $svgWidth = $longest / $pixelUnit + $margin * 2;
	$dialogWidth = ($svgWidth > 1000 ) ? 1050 : ($svgWidth < 550) ? 600 : $svgWidth + 50;
	$graphic =~ s/\$svgWidth/$svgWidth/g;
	$graphic =~ s/&/&amp;/g;

	open (SVGFILE,">$commoncfg->{TMPDIR}/GPM-$assembly[2]-$assembly[3]-Chr$chr.svg") or die "can't open file: $commoncfg->{TMPDIR}/GPM-$assembly[2]-$assembly[3]-Chr$chr.svg";
	print SVGFILE $graphic;
	close (SVGFILE);
	`gzip -f $commoncfg->{TMPDIR}/GPM-$assembly[2]-$assembly[3]-Chr$chr.svg`;

	$assemblyChrDetails ="
	<ul class='assemblyChrListMenu' style='width: 120px;float: left; margin-right: .3em; white-space: nowrap;'>
		<li><a><span class='ui-icon ui-icon-wrench'></span>You Can</a>
			<ul style='z-index: 1000;white-space: nowrap;'>
				<li><a><span class='ui-icon ui-icon-disk'></span>Download</a>
					<ul style='z-index: 1000;white-space: nowrap;'>
						<li><a href='download.cgi?assemblyId=$assemblyId&chr=$chr' target='hiddenFrame'><span class='ui-icon ui-icon-bullet'></span>Ctg PseudoMolecule</a></li>
						<li><a href='download.cgi?assemblyId=$assemblyId&chr=$chr&unit=chr' target='hiddenFrame'><span class='ui-icon ui-icon-bullet'></span>Chr PseudoMolecule</a></li>
						<li><a href='download.cgi?assemblyIdForAgp=$assemblyId&chr=$chr' target='hiddenFrame'><span class='ui-icon ui-icon-bullet'></span>Ctg-Seq AGP</a></li>
						<li><a href='download.cgi?assemblyIdForAgp=$assemblyId&chr=$chr&unit=chr&element=ctg' target='hiddenFrame'><span class='ui-icon ui-icon-bullet'></span>Chr-Ctg AGP</a></li>
						<li><a href='download.cgi?assemblyIdForAgp=$assemblyId&chr=$chr&unit=chr' target='hiddenFrame'><span class='ui-icon ui-icon-bullet'></span>Chr-Seq AGP</a></li>
						<li><a href='$commoncfg->{TMPURL}/GPM-$assembly[2]-$assembly[3]-Chr$chr.svg.gz' target='hiddenFrame'><span class='ui-icon ui-icon-bullet'></span>SVG</a></li>
					</ul>
				</li>
				<li><a onclick='printDiv(\"assemblyChrListForAlignment$assemblyId$$\")'><span class='ui-icon ui-icon-print'></span>Print</a></li>
			</ul>
		</li>
	</ul>
	<ul class='assemblyChrListMenu' style='width: 150px;float: left; margin-right: .3em; white-space: nowrap;'>
		<li><a>About Chr $chr</a>
			<ul style='z-index: 1000;white-space: nowrap;'>
				<li>Total Seqs: $totalSeqs</li>
				<li>Estimated Length: $totalLength->{'assembly'} bp (without gaps among contigs)</li>
			</ul>
		</li>
	</ul>
	<div style='float: right; margin-right: .3em; white-space: nowrap;'>
		<label for='companionAssemblyId'><b>Companion Assembly</b></label><select class='ui-widget-content ui-corner-all' name='companionAssemblyId' id='companionAssemblyId' onchange='closeViewer();openViewer(\"assemblyChrView.cgi?assemblyId=$assemblyId&chr=$chr&companionAssemblyId=\"+this.value);'>$companionAssemblyList</select>
		<label for='pixelUnit'><b>Pixel Unit</b></label><select class='ui-widget-content ui-corner-all' name='pixelUnit' id='pixelUnit' onchange='closeViewer();openViewer(\"assemblyChrView.cgi?assemblyId=$assemblyId&chr=$chr&pixelUnit=\"+this.value);'>$pixelUnitList</select>
		<label for='tickLength'><b>Tick Length</b></label><select class='ui-widget-content ui-corner-all' name='tickLength' id='tickLength' onchange='closeViewer();openViewer(\"assemblyChrView.cgi?assemblyId=$assemblyId&chr=$chr&tickLength=\"+this.value);'>$tickLengthList</select>
		<label for='alignmentLength'><b>Alignment Length</b></label><select class='ui-widget-content ui-corner-all' name='alignmentLength' id='alignmentLength' onchange='closeViewer();openViewer(\"assemblyChrView.cgi?assemblyId=$assemblyId&chr=$chr&alignmentLength=\"+this.value);'>$alignmentLengthList</select>
	</div>
	<br><br>
	<div id='assemblyChrListForAll$assemblyId$$'>
	<div id='assemblyChrListForCtg$assemblyId$$'>
	<div id='assemblyChrListForSeq$assemblyId$$'>
	<div id='assemblyChrListForAlignment$assemblyId$$'>
	$graphic
	</div>
	</div>
	</div>
	</div>";

	undef $/;# enable slurp mode
	my $html = <DATA>;
	$html =~ s/\$assemblyChrDetails/$assemblyChrDetails/g;
	$html =~ s/\$dialogWidth/$dialogWidth/g;
	$html =~ s/\$assemblyId/$assemblyId/g;
	$html =~ s/\$assemblyChr/$chr/g;
	$html =~ s/\$scrollLeft/$scrollLeft/g;
	$html =~ s/\$pixelUnit/$pixelUnit/g;
	$html =~ s/\$\$/$$/g;

	print header(-cookie=>[cookie(-name=>"companionAssembly$assemblyId",-value=>$companionAssemblyId),cookie(-name=>"pixelUnit",-value=>$pixelUnit),cookie(-name=>"tickLength",-value=>$tickLength),cookie(-name=>"alignmentLength",-value=>$alignmentLength)]);
	print $html;
}
else
{
	print header(-type=>'text/html',-status=>'402 Invalid operation');
	exit;
}


__DATA__
$assemblyChrDetails
<ul id="optionsForCtg" class="ui-helper-hidden" style='z-index: 1000;white-space: nowrap;'>
    <li data-command="edit"><a href="#"><span class="ui-icon ui-icon-pencil"></span>Edit Contig</a></li>
    <li data-command="flip"><a href="#"><span class="ui-icon ui-icon-arrowreturnthick-1-w"></span>Flip Contig</a></li>
    <li data-command="fill"><a href="#"><span class="ui-icon ui-icon-link"></span>Gap Filler</a></li>
    <li data-command="comment"><a href="#"><span class="ui-icon ui-icon-comment"></span>Comment</a></li>
</ul>
<ul id="optionsForSeq" class="ui-helper-hidden" style='z-index: 1000;white-space: nowrap;'>
    <li><a><span class="ui-icon ui-icon-transfer-e-w"></span>Alignment</a>
		<ul>
		    <li data-command="blast"><a href="#"><span class="ui-icon ui-icon-bullet"></span>BLAST2SEQ</a>
			<li data-command="align"><a href="#"><span class="ui-icon ui-icon-bullet"></span>Alignment Checker</a></li>
		</ul>
    </li>
    <li data-command="view"><a href="#"><span class="ui-icon ui-icon-clipboard"></span>View Details</a></li>
</ul>
<ul id="optionsForAlignment" class="ui-helper-hidden" style='z-index: 1000;white-space: nowrap;'>
    <li data-command="hide"><a href="#"><span class="ui-icon ui-icon-cancel"></span>Hide Alignment</a></li>
    <li data-command="view"><a href="#"><span class="ui-icon ui-icon-clipboard"></span>View Alignment</a></li>
</ul>
<ul id="optionsForAll" class="ui-helper-hidden" style='z-index: 1000;white-space: nowrap;'>
    <li><a><span class="ui-icon ui-icon-transfer-e-w"></span>Alignment</a>
		<ul>
		    <li data-command="blast"><a href="#"><span class="ui-icon ui-icon-bullet"></span>BLAST2SEQ</a>
			<li data-command="align"><a href="#"><span class="ui-icon ui-icon-bullet"></span>Alignment Checker</a></li>
		</ul>
    </li>
    <li data-command="fill"><a href="#"><span class="ui-icon ui-icon-link"></span>Gap Filler</a></li>
</ul>
<script>
buttonInit();
$( ".assemblyChrListMenu" ).menu();
$("#assemblyChrListForCtg$assemblyId$$").contextmenu({
    delegate: ".hasmenuForCtg",
	menu: "#optionsForCtg",
	position: {my: "left bottom", at: "right top"},
    select: function(event, ui) {
		if (ui.cmd == 'edit')
		{
			openDialog(ui.target.attr('edit-url'));
		}
		if (ui.cmd == 'flip')
		{
			var h=confirm("Are you sure to flip this?");
			if (h==true)
			{
				loaddiv('hiddenDiv',ui.target.attr('flip-url'));
			}
		}
		if (ui.cmd == 'fill')
		{
			openDialog(ui.target.attr('fill-url'));
		}
		if (ui.cmd == 'comment')
		{
			openDialog(ui.target.attr('comment-url'));
		}
    }
});
$( ".draggable" )
	.draggable({ axis: "x",
		stop: function( event, ui ) {
			var position = parseInt (ui.position.left * $pixelUnit);
			loaddiv('hiddenDiv',this.getAttribute('move-url') + position + '&scrollLeft=' + ui.position.left);
			}
		})
	.bind('mousedown', function(event, ui){
		// bring target to front
		$(event.target.parentElement).append( event.target );
		})
	.bind('drag', function(event, ui){
		// update coordinates manually, since top/left style props don't work on SVG
		event.target.setAttribute('x', ui.position.left);
		event.target.setAttribute('y', ui.position.top);
		});
$("#assemblyChrListForSeq$assemblyId$$").contextmenu({
    delegate: ".hasmenuForSeq",
	menu: "#optionsForSeq",
	position: {my: "left bottom", at: "right top"},
    select: function(event, ui) {
		if (ui.cmd == 'blast')
		{
			openDialog(ui.target.attr('blast-url'));
		}
 		if (ui.cmd == 'align')
		{
			openDialog(ui.target.attr('align-url'));
		}
		if (ui.cmd == 'view')
		{
			openDialog(ui.target.attr('view-url'));
		}
    }
});
$("#assemblyChrListForAlignment$assemblyId$$").contextmenu({
    delegate: ".hasmenuForAlignment",
	menu: "#optionsForAlignment",
    select: function(event, ui) {
		if (ui.cmd == 'hide')
		{
			var h=confirm("Are you sure to hide this?");
			if (h==true)
			{
				loaddiv('hiddenDiv',ui.target.attr('hide-url'));
			}
		}
		if (ui.cmd == 'view')
		{
			openDialog(ui.target.attr('view-url'));
		}
    }
});
$("#assemblyChrListForAll$assemblyId$$").contextmenu({
    delegate: ".hasmenuForAll",
	menu: "#optionsForAll",
    select: function(event, ui) {
		if (ui.cmd == 'blast')
		{
			openDialog(ui.target.attr('blast-url'));
		}
 		if (ui.cmd == 'align')
		{
			openDialog(ui.target.attr('align-url'));
		}
 		if (ui.cmd == 'fill')
		{
			openDialog(ui.target.attr('fill-url'));
		}
   }
});
$('#viewer').dialog("option", "title", "GPM Viewer: Chromosome $assemblyChr");
$('#viewer').dialog("option", "width", "$dialogWidth");
$('#viewer').scrollLeft($scrollLeft);
</script>
