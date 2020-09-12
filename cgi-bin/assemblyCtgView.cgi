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
my $assemblyCtgId = param ('assemblyCtgId') || '';
my $scrollLeft = param ('scrollLeft') || '0';
my $highlight = param ('highlight') || '';
my $assemblyCtgDetails = '';
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
my $barSpacing = 300; #space between two bars
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
my $gapLength = 100;
if ($assemblyCtgId)
{
	my $assemblyCtg=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$assemblyCtg->execute($assemblyCtgId);
	my @assemblyCtg = $assemblyCtg->fetchrow_array();
	my $assemblyPreCtgButton = '';
	my $assemblyCurrentCtgButton = "{ text: 'Reload - Ctg$assemblyCtg[2]', click: function() { closeViewer();openViewer('assemblyCtgView.cgi?assemblyCtgId=$assemblyCtg[0]');} }";
	my $assemblyNextCtgButton = '';
	my $assemblyCtgList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblyCtg' AND o = ? AND x = ? ORDER BY y");
	$assemblyCtgList->execute($assemblyCtg[3],$assemblyCtg[4]);
	while (my @assemblyCtgList = $assemblyCtgList->fetchrow_array())
	{
		$assemblyPreCtgButton = "{ text: 'Previous - Ctg$assemblyCtgList[2]', click: function() { closeViewer();openViewer('assemblyCtgView.cgi?assemblyCtgId=$assemblyCtgList[0]');} }" if ($assemblyCtgList[5] < $assemblyCtg[5]);
		$assemblyNextCtgButton = "{ text: 'Next - Ctg$assemblyCtgList[2]', click: function() { closeViewer();openViewer('assemblyCtgView.cgi?assemblyCtgId=$assemblyCtgList[0]');} }" if ($assemblyCtgList[5] > $assemblyCtg[5]);
		last if ($assemblyNextCtgButton);
	}
	my $assemblySwitchCtgButton = '';
	if ($assemblyPreCtgButton) 
	{
		if($assemblyNextCtgButton)
		{
			$assemblySwitchCtgButton = "$assemblyPreCtgButton,$assemblyCurrentCtgButton,$assemblyNextCtgButton";
		}
		else
		{
			$assemblySwitchCtgButton = "$assemblyPreCtgButton,$assemblyCurrentCtgButton";
		}
	}
	else
	{
		if($assemblyNextCtgButton)
		{
			$assemblySwitchCtgButton = "$assemblyCurrentCtgButton,$assemblyNextCtgButton";
		}
		else
		{
			$assemblySwitchCtgButton = "$assemblyCurrentCtgButton";
		}
	}	
	my $assembly=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$assembly->execute($assemblyCtg[3]);
	my @assembly = $assembly->fetchrow_array();
	my $assemblySequenceName;
	my $assemblySequenceGapType;
	my $assemblySequenceId;
	my $assemblySequenceData;
	my $assemblySequenceType;
	my $assemblySeqLength;
	my $assemblySeqOrient;
	my $assemblySeqStart;
	my $assemblySeqEnd;
	my $assemblySeqOrder;
	my $firstAssemblySeq = "";
	my $lastAssemblySeq = "";	
	my $assemblySeqHide;
	my $assemblySeqFpcCtg;
	my $order = 0;
	foreach (split ",", $assemblyCtg[8])
	{
		next unless ($_);
		my $hide = ($_ =~ /^-/) ? 1 : 0;
		$_ =~ s/[^a-zA-Z0-9]//g;
		$assemblySeqHide->{$_} = $hide;
		$firstAssemblySeq = $_ if ($hide eq 0 && $firstAssemblySeq eq "");
		$lastAssemblySeq = $_ if ($hide eq 0);
		$order++;
		$assemblySeqOrder->{$_} = $order;
		my $assemblySeqList=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$assemblySeqList->execute($_);
		my @assemblySeqList = $assemblySeqList->fetchrow_array();
		$assemblySequenceName->{$assemblySeqList[0]} = $assemblySeqList[2];
		$assemblySequenceGapType->{$assemblySeqList[0]} = $assemblySeqList[4];		
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
	}

	# create an SVG object
	$svg= SVG->new(width=>'$svgWidth',height=>'$svgHeight'); # set width and height after knowing the size
	$svg->rectangle(
		id    => $$,
		class   => 'hasmenuForAll',
		style => { stroke => 'white',
					fill => 'white'
					},
		x     => 0,
		y     => 0,
		width => '$svgWidth',
		height=> '$svgHeight',
		'edit-url' => "assemblyCtgEdit.cgi?assemblyCtgId=$assemblyCtgId",
		'flip-url' => "assemblyCtgFlip.cgi?assemblyCtgId=$assemblyCtgId",
		'filter-url' => "assemblyCtgOverlap.cgi?assemblyCtgId=$assemblyCtgId",
		'reset-url' => "assemblyCtgReset.cgi?assemblyCtgId=$assemblyCtgId",
		'blast-url' => "blastTwoseqForm.cgi?assemblyId=$assemblyCtg[3]",		
		'align-url' => "alignmentCheckForm.cgi?assemblyId=$assemblyCtg[3]"		
	);

	# use explicit element constructor to generate a group element
    my $ruler=$svg->group(
        id    => 'ruler',
        style => { stroke=>'black'}
    );
	# use explicit element constructor to generate a group element
    my $assemblyCtgSeq=$svg->group(
        id    => 'assemblyCtgSeq'
    );

	# use explicit element constructor to generate a group element
    my $assemblyCtgSeqAlignment=$svg->group(
        id    => 'assemblyCtgSeqAlignment'
    );

	# use explicit element constructor to generate a group element
    my $assemblyCtgSeqHidden=$svg->group(
        id    => 'assemblyCtgSeqHidden'
    );

	my @assemblySeq = sort { $assemblySeqOrder->{$a} <=> $assemblySeqOrder->{$b} } keys %$assemblySeqOrder;
	my $totalSeqs = @assemblySeq;
    my $preSeq = '';
    my $assemblySeqLeftEnd;
    my $assemblySeqRightEnd;
	my $assemblySeqMaxEnd = 0;
	my $assemblyCtgLength = 0;
	my $hiddenSeqCount = 0;
	my $filterLength=0;
	my @hiddenSeqPosition;
	my $assemblySeqCol;
    my $seqCount = 0;
	my $lastComponentType = '';

    for my $currentSeq (@assemblySeq)
    {
		my $sequenceDetails = decode_json $assemblySequenceData->{$currentSeq};
		$sequenceDetails->{'id'} = '' unless (exists $sequenceDetails->{'id'});
		$sequenceDetails->{'description'} = '' unless (exists $sequenceDetails->{'description'});
		$sequenceDetails->{'sequence'} = '' unless (exists $sequenceDetails->{'sequence'});
		$sequenceDetails->{'sequence'} =~ tr/a-zA-Z/N/c; #replace nonword characters.;
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
			my $hiddenSeqBarX = $margin + $assemblySeqLeftEnd->{$currentSeq} / $pixelUnit;
			$assemblyCtgSeqHidden->rectangle(
					id    => $currentSeq,
					onclick => "closeDialog();openDialog('seqView.cgi?seqId=$assemblySequenceId->{$currentSeq}')",
					class   => 'hasmenuForSeq',
					x     => $hiddenSeqBarX,
					y     => $barY + $hiddenSeqPosition + $barHeight * $hiddenRow * 2,
					width => $assemblySeqLength->{$currentSeq} / $pixelUnit,
					height=> $barHeight,
					style => { stroke => ($currentSeq eq $highlight) ? 'red' : 'grey',
								fill => 'lightgrey',
								opacity => 0.5
								},
					'flip-url' => "assemblySeqFlip.cgi?assemblySeqId=$currentSeq&scrollLeft=$hiddenSeqBarX",
					'hide-url' => "assemblySeqHide.cgi?assemblySeqId=$currentSeq&scrollLeft=$hiddenSeqBarX",
					'edit-url' => "assemblySeqEditForm.cgi?assemblySeqId=$currentSeq&scrollLeft=$hiddenSeqBarX",
					'break-url'=> "assemblyCtgBreakForm.cgi?assemblyCtgId=$assemblyCtgId&assemblySeqId=$currentSeq&scrollLeft=$hiddenSeqBarX",
					'move-url' => "assemblyCtgEdit.cgi?assemblyCtgId=$assemblyCtgId&scrollLeft=$hiddenSeqBarX",
					'delete-url' => "itemDelete.cgi?items=$currentSeq&scrollLeft=$hiddenSeqBarX",
					'blast-url' => "blastTwoseqForm.cgi?assemblyId=$assemblyCtg[3]&seqOne=$assemblySequenceId->{$currentSeq}",		
					'align-url' => "alignmentCheckForm.cgi?assemblyId=$assemblyCtg[3]&seqOne=$assemblySequenceId->{$currentSeq}",	
					'view-url' => "seqView.cgi?seqId=$assemblySequenceId->{$currentSeq}",
					'colNumber'=> 0
				);
			$scrollLeft = $assemblySeqLeftEnd->{$currentSeq} / $pixelUnit if ($currentSeq eq $highlight);
			my $hiddenTextY = $barY + $hiddenSeqPosition + $barHeight * $hiddenRow * 2 + $barHeight - 2;
			$assemblyCtgSeqHidden->text(
				id      => 'seqName'.$currentSeq,
				onclick => "closeDialog();openDialog('seqView.cgi?seqId=$assemblySequenceId->{$currentSeq}')",
				class   => 'hasmenuForSeq',
				x       => $hiddenSeqBarX,
				y       => $hiddenTextY,
				style   => {
					'font-size'   => $textFontSize,
					'stroke'        => ($sequenceDetails->{'gapList'}) ? 'pink' : 'grey'
				},
				'flip-url' => "assemblySeqFlip.cgi?assemblySeqId=$currentSeq&scrollLeft=$hiddenSeqBarX",
				'hide-url' => "assemblySeqHide.cgi?assemblySeqId=$currentSeq&scrollLeft=$hiddenSeqBarX",
				'edit-url' => "assemblySeqEditForm.cgi?assemblySeqId=$currentSeq&scrollLeft=$hiddenSeqBarX",
				'break-url'=> "assemblyCtgBreakForm.cgi?assemblySeqId=$currentSeq&scrollLeft=$hiddenSeqBarX",
				'move-url' => "assemblyCtgEdit.cgi?assemblyCtgId=$assemblyCtgId&scrollLeft=$hiddenSeqBarX",
				'delete-url' => "itemDelete.cgi?items=$currentSeq&scrollLeft=$hiddenSeqBarX",
				'blast-url' => "blastTwoseqForm.cgi?assemblyId=$assemblyCtg[3]&seqOne=$assemblySequenceId->{$currentSeq}",		
				'align-url' => "alignmentCheckForm.cgi?assemblyId=$assemblyCtg[3]&seqOne=$assemblySequenceId->{$currentSeq}",	
				'view-url' => "seqView.cgi?seqId=$assemblySequenceId->{$currentSeq}",
				'colNumber'=> 0
			)->cdata($assemblySeqOrder->{$currentSeq}. "." .$assemblySequenceName->{$currentSeq}."  FPC:".$assemblySeqFpcCtg->{$currentSeq});
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

		my $goodSequenceX = ($assemblySeqOrient->{$currentSeq} eq "+") ? $margin + ($assemblySeqLeftEnd->{$currentSeq} + $assemblySeqStart->{$currentSeq} - 1) / $pixelUnit : $margin + ($assemblySeqRightEnd->{$currentSeq} - $assemblySeqEnd->{$currentSeq} + 1) / $pixelUnit ;
		$assemblyCtgSeq->rectangle(
				x     => $goodSequenceX,
				y     => $barY + ($barHeight + $barSpacing) * $goodCol,
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
				my $gapX = ($assemblySeqOrient->{$currentSeq} eq "+") ? $margin + ($assemblySeqLeftEnd->{$currentSeq} + $gapStart - 1) / $pixelUnit : $margin + ($assemblySeqRightEnd->{$currentSeq} - $gapEnd + 1) / $pixelUnit ;
				$assemblyCtgSeq->rectangle(
						x     => $gapX,
						y     => $barY + ($barHeight + $barSpacing) * $goodCol,
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
				my $filterX = ($assemblySeqOrient->{$currentSeq} eq "+") ? $margin + ($assemblySeqLeftEnd->{$currentSeq} + $filterStart - 1) / $pixelUnit : $margin + ($assemblySeqRightEnd->{$currentSeq} - $filterEnd + 1) / $pixelUnit ;
				$assemblyCtgSeq->rectangle(
						x     => $filterX,
						y     => $barY + ($barHeight + $barSpacing) * $goodCol,
						width => ($filterEnd - $filterStart + 1) / $pixelUnit,
						height=> $barHeight,
						style => { stroke => 'black',
									fill => 'lightgrey'
									},
						id    => "Filter$currentSeq$filterStart$filterEnd$$"
					);
				next if ($assemblySeqStart->{$currentSeq} > $filterEnd);
				next if ($assemblySeqEnd->{$currentSeq} < $filterStart);
				if ($assemblySeqStart->{$currentSeq} >= $filterStart && $assemblySeqEnd->{$currentSeq} <= $filterEnd)
				{
					$filterLength += $assemblySeqEnd->{$currentSeq} - $assemblySeqStart->{$currentSeq} + 1;
				}
				elsif ($assemblySeqStart->{$currentSeq} >= $filterStart && $assemblySeqStart->{$currentSeq} <= $filterEnd)
				{
					$filterLength += $filterEnd - $assemblySeqStart->{$currentSeq} + 1;

				}
				elsif ($assemblySeqEnd->{$currentSeq} >= $filterStart && $assemblySeqEnd->{$currentSeq} <= $filterEnd)
				{
					$filterLength += $assemblySeqEnd->{$currentSeq} - $filterStart + 1;
				}
				else
				{
					$filterLength += $filterEnd - $filterStart + 1;
				}
			}
		}

		my $seqBarX = $margin + $assemblySeqLeftEnd->{$currentSeq} / $pixelUnit;
		$assemblyCtgSeq->rectangle(
				id    => $currentSeq,
				onclick => "closeDialog();openDialog('seqView.cgi?seqId=$assemblySequenceId->{$currentSeq}')",
				class   => 'hasmenuForSeq',
				x     => $seqBarX,
				y     => $barY + ($barHeight + $barSpacing) * $goodCol,
				width => $assemblySeqLength->{$currentSeq} / $pixelUnit,
				height=> $barHeight,
				style => { stroke => ($currentSeq eq $highlight) ? 'red' : 'black',
							fill => 'lightgrey',
							'fill-opacity' =>  0.5
							},
				'flip-url' => "assemblySeqFlip.cgi?assemblySeqId=$currentSeq&scrollLeft=$seqBarX",
				'hide-url' => "assemblySeqHide.cgi?assemblySeqId=$currentSeq&scrollLeft=$seqBarX",
				'edit-url' => "assemblySeqEditForm.cgi?assemblySeqId=$currentSeq&scrollLeft=$seqBarX",
				'break-url'=> "assemblyCtgBreakForm.cgi?assemblySeqId=$currentSeq&scrollLeft=$seqBarX",
				'move-url' => "assemblyCtgEdit.cgi?assemblyCtgId=$assemblyCtgId&scrollLeft=$seqBarX",
				'delete-url' => "itemDelete.cgi?items=$currentSeq&scrollLeft=$seqBarX",
				'blast-url' => "blastTwoseqForm.cgi?assemblyId=$assemblyCtg[3]&seqOne=$assemblySequenceId->{$currentSeq}",		
				'align-url' => "alignmentCheckForm.cgi?assemblyId=$assemblyCtg[3]&seqOne=$assemblySequenceId->{$currentSeq}",	
				'view-url' => "seqView.cgi?seqId=$assemblySequenceId->{$currentSeq}",
				'colNumber'=> $goodCol
			);
		$scrollLeft = $assemblySeqLeftEnd->{$currentSeq} / $pixelUnit if ($currentSeq eq $highlight);

		$assemblyCtgSeq->text(
			id      => 'seqName'.$currentSeq,
			x       => $seqBarX,
			y       => $barY + ($barHeight + $barSpacing) * $goodCol + $barHeight - 2,
			style   => {
				'font-size'   => $textFontSize,
				'stroke'        => ($sequenceDetails->{'gapList'}) ? 'red' : 'black'
			},
			onclick => "closeDialog();openDialog('seqView.cgi?seqId=$assemblySequenceId->{$currentSeq}')",
			class   => 'hasmenuForSeq',
			'flip-url' => "assemblySeqFlip.cgi?assemblySeqId=$currentSeq&scrollLeft=$seqBarX",
			'hide-url' => "assemblySeqHide.cgi?assemblySeqId=$currentSeq&scrollLeft=$seqBarX",
			'edit-url' => "assemblySeqEditForm.cgi?assemblySeqId=$currentSeq&scrollLeft=$seqBarX",
			'break-url'=> "assemblyCtgBreakForm.cgi?assemblySeqId=$currentSeq&scrollLeft=$seqBarX",
			'move-url' => "assemblyCtgEdit.cgi?assemblyCtgId=$assemblyCtgId&scrollLeft=$seqBarX",
			'delete-url' => "itemDelete.cgi?items=$currentSeq&scrollLeft=$seqBarX",
			'blast-url' => "blastTwoseqForm.cgi?assemblyId=$assemblyCtg[3]&seqOne=$assemblySequenceId->{$currentSeq}",		
			'align-url' => "alignmentCheckForm.cgi?assemblyId=$assemblyCtg[3]&seqOne=$assemblySequenceId->{$currentSeq}",	
			'view-url' => "seqView.cgi?seqId=$assemblySequenceId->{$currentSeq}",
			'colNumber'=> $goodCol
		)->cdata($assemblySeqOrder->{$currentSeq}. "." . $assemblySequenceName->{$currentSeq} . "(" .$assemblySeqOrient->{$currentSeq}.")  FPC:".$assemblySeqFpcCtg->{$currentSeq});

    	if($preSeq) #alignment
    	{
			my $preY = 0;
			my $currentY = 0;
			if($assemblySeqCol->{$preSeq} < $assemblySeqCol->{$currentSeq})
			{
				$preY = $barY + ($barHeight + $barSpacing) * $assemblySeqCol->{$preSeq} + $barHeight + 1;
				$currentY =  $barY + ($barHeight + $barSpacing) * $assemblySeqCol->{$currentSeq} - 1;
			}
			else
			{
				$preY = $barY + ($barHeight + $barSpacing) * $assemblySeqCol->{$preSeq} - 1;
				$currentY =  $barY + ($barHeight + $barSpacing) * $assemblySeqCol->{$currentSeq} + $barHeight + 1;
			}


			my $queryDir;
			for (my $position = 0; $position < length($assemblySequenceId->{$preSeq}); $position += 2)
			{
				$queryDir .= "/q". substr($assemblySequenceId->{$preSeq},$position,2);
			}
			my $subjectDir;
			for (my $position = 0; $position < length($assemblySequenceId->{$currentSeq}); $position += 2)
			{
				$subjectDir .= "/s". substr($assemblySequenceId->{$currentSeq},$position,2);
			}
			if(-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$assemblySequenceId->{$preSeq}-$assemblySequenceId->{$currentSeq}.tbl")
			{
				my $alignmentCount = 0;
				open (TBL,"$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$assemblySequenceId->{$preSeq}-$assemblySequenceId->{$currentSeq}.tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$assemblySequenceId->{$preSeq}-$assemblySequenceId->{$currentSeq}.tbl";
				while(<TBL>)
				{
					chop;
					/^#/ and next;
					$alignmentCount++;
					my @alignments = split("\t",$_);
					next if ($alignments[12] > 0);
					next if ($alignments[3] < $alignmentLength); #skip if alignment shorter than alignmentLength
					my $preXOne = 0;
					my $preXTwo = 0;
					if($assemblySeqOrient->{$preSeq} eq "+")
					{
						$preXOne = $margin + ($assemblySeqLeftEnd->{$preSeq} + $alignments[6] - 1) / $pixelUnit;
						$preXTwo = $margin + ($assemblySeqLeftEnd->{$preSeq} + $alignments[7] - 1) / $pixelUnit;
					}
					else
					{
						$preXOne = $margin + ($assemblySeqRightEnd->{$preSeq} - $alignments[6] + 1) / $pixelUnit;
						$preXTwo = $margin + ($assemblySeqRightEnd->{$preSeq} - $alignments[7] + 1) / $pixelUnit;
					}
					my $currentXOne = 0;
					my $currentXTwo = 0;
					if($assemblySeqOrient->{$currentSeq} eq "+")
					{
						$currentXOne = $margin + ($assemblySeqLeftEnd->{$currentSeq} + $alignments[8] - 1) / $pixelUnit;
						$currentXTwo = $margin + ($assemblySeqLeftEnd->{$currentSeq} + $alignments[9] - 1) / $pixelUnit;
					}
					else
					{
						$currentXOne = $margin + ($assemblySeqRightEnd->{$currentSeq} - $alignments[8] + 1) / $pixelUnit;
						$currentXTwo = $margin + ($assemblySeqRightEnd->{$currentSeq} - $alignments[9] + 1) / $pixelUnit;
					}		
					my $xv = [$preXOne,$preXTwo,$currentXTwo,$currentXOne];
					my $yv = [$preY,$preY,$currentY,$currentY];
					my $points = $assemblyCtgSeqAlignment->get_path(
						x=>$xv,
						y=>$yv,
						-type=>'polygon'
					);
					$assemblyCtgSeqAlignment->polygon(
						%$points,
						id=>$alignments[0].$alignments[1].$alignmentCount,
						onclick => "closeDialog();openDialog('alignmentView.cgi?alignmentId=$alignments[0]-$alignments[1]-$alignmentCount')",
						class=>'hasmenuForAlignment',
						style=>{ stroke => 'red',
							fill => 'yellow',
							opacity => 0.5
							},
						'overlap-url' => "assemblyCtgOverlap.cgi?assemblyCtgId=$assemblyCtgId&alignmentId=$alignments[0]-$alignments[1]-$alignmentCount&scrollLeft=$preXOne",
						'hide-url' => "assemblyAlignmentHide.cgi?alignmentId=$alignments[0]-$alignments[1]-$alignmentCount&assemblyCtgId=$assemblyCtgId&scrollLeft=$preXOne",
						'view-url' => "alignmentView.cgi?alignmentId=$alignments[0]-$alignments[1]-$alignmentCount",
						'move-url' => "assemblyCtgEdit.cgi?assemblyCtgId=$assemblyCtgId&scrollLeft=$preXOne"
					);
				}
				close(TBL);
			}
		}
    	$assemblySeqMaxEnd = $assemblySeqRightEnd->{$currentSeq} if($assemblySeqMaxEnd < $assemblySeqRightEnd->{$currentSeq});
 
     	if($lastComponentType ne 'U')
    	{
			if($currentSeq ne $firstAssemblySeq)
			{
				if ($assemblySequenceGapType->{$currentSeq} eq 1 || $assemblySequenceGapType->{$currentSeq} eq 3 || $assemblySequenceGapType->{$currentSeq} eq 4 || $assemblySequenceGapType->{$currentSeq} eq 6 || $assemblySequenceGapType->{$currentSeq} eq 7 || $assemblySequenceGapType->{$currentSeq} eq 8) # add 5' 100 Ns
				{
					$assemblyCtgLength += $gapLength;
					$lastComponentType = 'U';
					#graphic to be added
				}
			}
    	}

    	$assemblyCtgLength += $assemblySeqEnd->{$currentSeq} - $assemblySeqStart->{$currentSeq} + 1;
    	$lastComponentType = 'D';

    	if($lastComponentType ne 'U')
    	{
			if($currentSeq ne $lastAssemblySeq)
			{
				if ($assemblySequenceGapType->{$currentSeq} eq 2 || $assemblySequenceGapType->{$currentSeq} eq 3 || $assemblySequenceGapType->{$currentSeq} eq 5 || $assemblySequenceGapType->{$currentSeq} eq 6 || $assemblySequenceGapType->{$currentSeq} eq 7 || $assemblySequenceGapType->{$currentSeq} eq 8) # add 5' 100 Ns
				{
					$assemblyCtgLength += $gapLength;
					$lastComponentType = 'U';
					#graphic to be added
				}
			}
    	}
    	$preSeq = $currentSeq;
    	$hiddenSeqCount = 0;
    	$seqCount++;
	}

	for (my $rulerNumber = 0;$rulerNumber <= $assemblySeqMaxEnd;$rulerNumber += $tickLength)
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
			y2    => '$svgHeight'
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

	
    $assemblyCtgLength = $assemblyCtgLength - $filterLength;
	my $updateAssemblyCtgLength=$dbh->do("UPDATE matrix SET barcode = $assemblyCtgLength WHERE id = $assemblyCtgId");
	$assemblyCtgLength = commify($assemblyCtgLength);

	# now render the SVG object, implicitly use svg namespace
	my $graphic = $svg->xmlify;
	my $svgWidth = $assemblySeqMaxEnd / $pixelUnit + $margin * 2;
	my $svgHeight = $barY + ($barHeight + $barSpacing) * 2 - $barSpacing + $barHeight + $margin;
	$graphic =~ s/\$svgWidth/$svgWidth/g;
	$graphic =~ s/\$svgHeight/$svgHeight/g;
	$graphic =~ s/&/&amp;/g;

	open (SVGFILE,">$commoncfg->{TMPDIR}/GPM-Ctg$assemblyCtg[2]-$assemblyCtgId.svg") or die "can't open file: $commoncfg->{TMPDIR}/GPM-Ctg$assemblyCtg[2]-$assemblyCtgId.svg";
	print SVGFILE $graphic;
	close (SVGFILE);
	`gzip -f $commoncfg->{TMPDIR}/GPM-Ctg$assemblyCtg[2]-$assemblyCtgId.svg`;

	$dialogWidth = ($svgWidth > 1000 ) ? 1050 : ($svgWidth < 550) ? 600 : $svgWidth + 50;
	$assemblyCtgDetails ="
	<ul class='assemblyCtgMenu' style='width: 100px;float: left; margin-right: .3em; white-space: nowrap;'>
		<li><a><span class='ui-icon ui-icon-wrench'></span>You Can</a>
			<ul style='z-index: 1000;white-space: nowrap;'>
				<li><a onclick='openDialog(\"assemblyCtgEdit.cgi?assemblyCtgId=$assemblyCtgId\");'><span class='ui-icon ui-icon-pencil'></span>Edit Ctg$assemblyCtg[2]</a></li>
				<li><a onclick='loaddiv(\"hiddenDiv\",\"assemblyCtgFlip.cgi?assemblyCtgId=$assemblyCtgId\")'><span class='ui-icon ui-icon-arrowreturnthick-1-w'></span>Flip Ctg$assemblyCtg[2]</a></li>
				<li><a onclick='deleteItem($assemblyCtgId);'><span class='ui-icon ui-icon-trash'></span>Delete Ctg$assemblyCtg[2]</a></li>
				<li><a><span class='ui-icon ui-icon-gripsmall-diagonal-se'></span>Redundancy</a>
					<ul style='z-index: 1000;white-space: nowrap;'>
						<li><a onclick='loaddiv(\"hiddenDiv\",\"assemblyCtgOverlap.cgi?assemblyCtgId=$assemblyCtgId\")'><span class='ui-icon ui-icon-bullet'></span>Filter</a></li>
						<li><a onclick='loaddiv(\"hiddenDiv\",\"assemblyCtgReset.cgi?assemblyCtgId=$assemblyCtgId\")'><span class='ui-icon ui-icon-bullet'></span>Reset</a></li>
					</ul>
				</li>
				<li><a><span class='ui-icon ui-icon-disk'></span>Download</a>
					<ul style='z-index: 1000;white-space: nowrap;'>
						<li><a href='download.cgi?assemblyCtgIdForAgp=$assemblyCtgId' target='hiddenFrame'><span class='ui-icon ui-icon-bullet'></span>AGP</a></li>
						<li><a href='download.cgi?assemblyCtgId=$assemblyCtgId' target='hiddenFrame'><span class='ui-icon ui-icon-bullet'></span>Raw Sequence</a></li>
						<li><a href='download.cgi?assemblyCtgId=$assemblyCtgId&pseudo=1' target='hiddenFrame'><span class='ui-icon ui-icon-bullet'></span>PseudoMolecule</a></li>
						<li><a href='$commoncfg->{TMPURL}/GPM-Ctg$assemblyCtg[2]-$assemblyCtgId.svg.gz' target='hiddenFrame'><span class='ui-icon ui-icon-bullet'></span>SVG</a></li>
					</ul>
				</li>
				<li><a onclick='printDiv(\"assemblyCtgForAlignment$assemblyCtgId$$\")'><span class='ui-icon ui-icon-print'></span>Print</a></li>
				<li><a>Learn about Ctg$assemblyCtg[2]</a>
					<ul style='z-index: 1000;white-space: nowrap;'>
						<li>Total Seqs: $totalSeqs</li>
						<li>Estimated Length: $assemblyCtgLength bp</li>
						<li>Chromosome: $assemblyCtg[4]</li>
					</ul>
				</li>
			</ul>
		</li>
	</ul>
	<div style='float: right; margin-right: .3em; white-space: nowrap;'>
		<label for='pixelUnit'><b>Pixel Unit</b></label><select class='ui-widget-content ui-corner-all' name='pixelUnit' id='pixelUnit' onchange='closeViewer();openViewer(\"assemblyCtgView.cgi?assemblyCtgId=$assemblyCtg[0]&pixelUnit=\"+this.value);'>$pixelUnitList</select>
		<label for='tickLength'><b>Tick Length</b></label><select class='ui-widget-content ui-corner-all' name='tickLength' id='tickLength' onchange='closeViewer();openViewer(\"assemblyCtgView.cgi?assemblyCtgId=$assemblyCtg[0]&tickLength=\"+this.value);'>$tickLengthList</select>
		<label for='alignmentLength'><b>Alignment Length</b></label><select class='ui-widget-content ui-corner-all' name='alignmentLength' id='alignmentLength' onchange='closeViewer();openViewer(\"assemblyCtgView.cgi?assemblyCtgId=$assemblyCtg[0]&alignmentLength=\"+this.value);'>$alignmentLengthList</select>
		<a onclick='openDialog(\"comment.cgi?itemId=$assemblyCtgId\");'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-comment'></span>Comments</a>
	</div>
	<br><br>
	<div id='assemblyCtgForAll$assemblyCtgId$$'>
	<div id='assemblyCtgForSeq$assemblyCtgId$$'>
	<div id='assemblyCtgForAlignment$assemblyCtgId$$'>
	$graphic
	</div>
	</div>
	</div>
	";

	undef $/;# enable slurp mode
	my $html = <DATA>;
	$html =~ s/\$assemblyCtgDetails/$assemblyCtgDetails/g;
	$html =~ s/\$dialogWidth/$dialogWidth/g;
	$html =~ s/\$assemblyCtgId/$assemblyCtgId/g;
	$html =~ s/\$assemblyCtg/$assemblyCtg[2]/g;
	$html =~ s/\$assemblySwitchCtgButton/$assemblySwitchCtgButton/g;
	$html =~ s/\$scrollLeft/$scrollLeft/g;
	$html =~ s/\$\$/$$/g;

	print header(-cookie=>[cookie(-name=>'assemblyCtg',-value=>$assemblyCtgId),cookie(-name=>"pixelUnit",-value=>$pixelUnit),cookie(-name=>"tickLength",-value=>$tickLength),cookie(-name=>"alignmentLength",-value=>$alignmentLength)]);
	print $html;
}
else
{
	print header(-type=>'text/html',-status=>'402 Invalid operation');
	exit;
}

__DATA__
$assemblyCtgDetails
<ul id="optionsForSeq" class="ui-helper-hidden" style='z-index: 1000;white-space: nowrap;'>
    <li data-command="flip"><a href="#"><span class="ui-icon ui-icon-arrowreturnthick-1-w"></span>Flip Sequence</a></li>
    <li data-command="hide"><a href="#"><span class="ui-icon ui-icon-lightbulb"></span>Hide/Show</a></li>
    <li data-command="edit"><a href="#"><span class="ui-icon ui-icon-carat-2-e-w"></span>Edit Sequence</a></li>
    <li data-command="break"><a href="#"><span class="ui-icon ui-icon-scissors"></span>Break Contig</a></li>
    <li data-command="move"><a href="#"><span class="ui-icon ui-icon-arrow-2-e-w"></span>Move Seq</a></li>
    <li data-command="delete"><a href="#"><span class="ui-icon ui-icon-trash"></span>Delete Seq</a></li>
    <li><a><span class="ui-icon ui-icon-transfer-e-w"></span>Alignment</a>
		<ul>
		    <li data-command="blast"><a href="#"><span class="ui-icon ui-icon-bullet"></span>BLAST2SEQ</a>
			<li data-command="align"><a href="#"><span class="ui-icon ui-icon-bullet"></span>Alignment Checker</a></li>
		</ul>
    </li>
    <li data-command="view"><a href="#"><span class="ui-icon ui-icon-clipboard"></span>View Details</a></li>
</ul>
<ul id="optionsForAlignment" class="ui-helper-hidden" style='z-index: 1000;white-space: nowrap;'>
    <li data-command="overlap"><a href="#"><span class="ui-icon ui-icon-gripsmall-diagonal-se"></span>Smart Redundancy Filter</a></li>
    <li data-command="hide"><a href="#"><span class="ui-icon ui-icon-cancel"></span>Hide Alignment</a></li>
    <li data-command="view"><a href="#"><span class="ui-icon ui-icon-clipboard"></span>View Alignment</a></li>
    <li data-command="move"><a href="#"><span class="ui-icon ui-icon-arrow-2-e-w"></span>Move Seq</a></li>
</ul>
<ul id="optionsForAll" class="ui-helper-hidden" style='z-index: 1000;white-space: nowrap;'>
    <li data-command="edit"><a href="#"><span class="ui-icon ui-icon-pencil"></span>Edit Ctg$assemblyCtg</a></li>
    <li data-command="flip"><a href="#"><span class="ui-icon ui-icon-arrowreturnthick-1-w"></span>Flip Ctg$assemblyCtg</a></li>
	<li><a><span class='ui-icon ui-icon-gripsmall-diagonal-se'></span>Redundancy</a>
		<ul>
			<li data-command="filter"><a href="#"><span class='ui-icon ui-icon-bullet'></span>Filter</a></li>
			<li data-command="reset"><a href="#"><span class='ui-icon ui-icon-bullet'></span>Reset</a></li>
		</ul>
	</li>
    <li><a><span class="ui-icon ui-icon-transfer-e-w"></span>Alignment</a>
		<ul>
		    <li data-command="blast"><a href="#"><span class="ui-icon ui-icon-bullet"></span>BLAST2SEQ</a>
			<li data-command="align"><a href="#"><span class="ui-icon ui-icon-bullet"></span>Alignment Checker</a></li>
		</ul>
    </li>
</ul>
<script>
buttonInit();
$( ".assemblyCtgMenu" ).menu();
$("#assemblyCtgForSeq$assemblyCtgId$$").contextmenu({
    delegate: ".hasmenuForSeq",
	menu: "#optionsForSeq",
	position: function(event, ui){
			if(ui.target.attr('colNumber') > 0)
			{
				return {my: "left bottom", at: "right top"};
			}
			else
			{
				return {my: "left top", at: "right bottom"};
			}
		},
    select: function(event, ui) {
		if (ui.cmd == 'flip')
		{
			var h=confirm("Are you sure to flip this?");
			if (h==true)
			{
				loaddiv('hiddenDiv',ui.target.attr('flip-url'));
			}
		}
		if (ui.cmd == 'hide')
		{
			loaddiv('hiddenDiv',ui.target.attr('hide-url'));
		}
		if (ui.cmd == 'edit')
		{
			openDialog(ui.target.attr('edit-url'));
		}
		if (ui.cmd == 'break')
		{
			openDialog(ui.target.attr('break-url'));
		}
		if (ui.cmd == 'move')
		{
			openDialog(ui.target.attr('move-url'));
		}
		if (ui.cmd == 'delete')
		{
			var h=confirm("Are you sure to delete this?");
			if (h==true)
			{
				loaddiv('hiddenDiv',ui.target.attr('delete-url'));
			}
		}
  		if(ui.cmd == 'blast')
		{
			openDialog(ui.target.attr('blast-url'));
		}
 		if(ui.cmd == 'align')
		{
			openDialog(ui.target.attr('align-url'));
		}
		if (ui.cmd == 'view')
		{
			openDialog(ui.target.attr('view-url'));
		}
    }
});
$("#assemblyCtgForAlignment$assemblyCtgId$$").contextmenu({
    delegate: ".hasmenuForAlignment",
	menu: "#optionsForAlignment",
    select: function(event, ui) {
		if (ui.cmd == 'overlap')
		{
			loaddiv('hiddenDiv',ui.target.attr('overlap-url'));
		}
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
		if (ui.cmd == 'move')
		{
			openDialog(ui.target.attr('move-url'));
		}
    }
});
$("#assemblyCtgForAll$assemblyCtgId$$").contextmenu({
    delegate: ".hasmenuForAll",
	menu: "#optionsForAll",
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
		if(ui.cmd == 'filter')
		{
			var h=confirm("Are you sure to filter this redundancy?");
			if (h==true)
			{
		 		loaddiv('hiddenDiv',ui.target.attr('filter-url'));
		 	}
		}
 		if(ui.cmd == 'reset')
		{
			var h=confirm("Are you sure to reset this redundancy?");
			if (h==true)
			{
		 		loaddiv('hiddenDiv',ui.target.attr('reset-url'));
		 	}
		}
 		if(ui.cmd == 'blast')
		{
			openDialog(ui.target.attr('blast-url'));
		}
 		if(ui.cmd == 'align')
		{
			openDialog(ui.target.attr('align-url'));
		}
   }
});
$('#viewer').dialog("option", "title", "GPM Viewer: Ctg$assemblyCtg");
$('#viewer').dialog("option", "width", "$dialogWidth");
$("#viewer" ).dialog( "option", "buttons", [ $assemblySwitchCtgButton ] );
$('#viewer').scrollLeft($scrollLeft);
</script>
