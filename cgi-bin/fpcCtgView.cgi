#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
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
my $fpcCtgId = param ('fpcCtgId') || '';
my $highlight = param ('highlight') || '';
my $scrollLeft = param ('scrollLeft') || '0';
my $fpcCtgDetails = '';
my $dialogWidth = 600;
if ($fpcCtgId)
{
	my $fpcCtgList=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$fpcCtgList->execute($fpcCtgId);
	my @fpcCtgList = $fpcCtgList->fetchrow_array();
	$fpcCtgList[8] =~ s/\n/<br>/g;
    my $textFontSize = 10;
    my $textFontWidth = 7;
	my $pixelUnit = 5;
	my $barY = 25;
	my $rulerY = 20;
	my $margin = 4;
	my $barHeight = 12;
	my $barSpacing = 5; #space between two bars
	my $unitLength = 10;
    my $textLength;
	my $fpcClone;
	my $fpcCloneSequenced;
	my $fpcCloneMTP;
	my $fpcCloneHighlighted;
	my $fpcCloneLeftEnd;
	my $fpcCloneRightEnd;
	my $fpcCloneMaxEnd = 0;
	my $fpcCloneList=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'fpcClone' AND o = ? AND MATCH (note) AGAINST (?)");
	$fpcCloneList->execute($fpcCtgList[3],$fpcCtgList[2]);
	while (my @fpcCloneList = $fpcCloneList->fetchrow_array())
	{
		$fpcClone->{$fpcCloneList[2]} = $fpcCloneList[8];
		$fpcCloneSequenced->{$fpcCloneList[2]} = $fpcCloneList[4];
		$fpcCloneMTP->{$fpcCloneList[2]} = $fpcCloneList[5];
		$fpcCloneHighlighted->{$fpcCloneList[2]} = $fpcCloneList[6];
		$fpcCloneList[8] =~ /Ends Left (\d*)/;
		$fpcCloneLeftEnd->{$fpcCloneList[2]} = $1;
		$fpcCloneList[8] =~ /Ends Right (\d*)/;
		$fpcCloneRightEnd->{$fpcCloneList[2]} = $1;
		$fpcCloneMaxEnd = $1 if ($1 > $fpcCloneMaxEnd);
		$textLength->{$fpcCloneList[2]} = length $fpcCloneList[2];
		$fpcCloneMaxEnd = $fpcCloneLeftEnd->{$fpcCloneList[2]}  + $textLength->{$fpcCloneList[2]} * $textFontWidth / $pixelUnit if($fpcCloneLeftEnd->{$fpcCloneList[2]}  + $textLength->{$fpcCloneList[2]} * $textFontWidth / $pixelUnit > $fpcCloneMaxEnd);
	}
	my @fpcClone = sort { $fpcCloneLeftEnd->{$a} <=> $fpcCloneLeftEnd->{$b} } keys %$fpcCloneLeftEnd;

	# create an SVG object
	$svg= SVG->new(width=>'$svgWidth',height=>'$svgHeight'); # set width and height after knowing the size

	# use explicit element constructor to generate a group element
    my $ruler=$svg->group(
        id    => 'ruler',
        style => { stroke=>'black'}
    );
	for (my $rulerNumber = 0;$rulerNumber < $fpcCloneMaxEnd;$rulerNumber += $unitLength)
	{
		#dash lines
		$ruler->line(
			id    => "rulerDash$rulerNumber",
			style => {
				stroke=> ($rulerNumber % (5*$unitLength) == 0) ? 'grey' : 'lightgrey',
				'stroke-dasharray' => '3,3'
				},
			x1    => $margin + $rulerNumber * $pixelUnit,
			y1    => $rulerY,
			x2    => $margin + $rulerNumber * $pixelUnit,
			y2    => '$svgHeight'
		);
		#region
		$ruler->line(
			id    => "rulerRegion$rulerNumber",
			x1    => $margin + $rulerNumber * $pixelUnit,
			y1    => $rulerY,
			x2    => $margin +  ($rulerNumber + $unitLength ) * $pixelUnit,
			y2    => $rulerY
		);
		#ticks
		$ruler->line(
			id    => "rulerTick$rulerNumber",
			x1    => $margin + $rulerNumber * $pixelUnit,
			y1    => ($rulerNumber % (5*$unitLength) == 0) ? $rulerY - 5 : $rulerY - 3,
			x2    => $margin + $rulerNumber * $pixelUnit,
			y2    => $rulerY
			
		);
		my $textX =$margin + $rulerNumber * $pixelUnit;
		$ruler->text(
			id      => "rulerNumber$rulerNumber",
			x       => $textX,
			y       => $rulerY - 2,
			style   => {
				'font-size'   => ($rulerNumber % (5*$unitLength) == 0) ? 11 : 9,
				'stroke'        => ($rulerNumber % (5*$unitLength) == 0) ? 'black' : 'grey'
			}
		)->cdata($rulerNumber);
	}

	# use explicit element constructor to generate a group element
    my $fpcCtgClone=$svg->group(
        id    => 'fpcCtgClone'
    );
	my $colCount=0;
	my $sequenced = 0;
	my @colPostion;
	push @colPostion, -1;
    for my $currentClone (@fpcClone)
    {
    	my $col = 0;
    	my $goodCol = @colPostion;
    	for(@colPostion)
    	{
    		if($_ < $fpcCloneLeftEnd->{$currentClone})
    		{
    			$goodCol = $col;
    			last;
    		}
    		$col++;
    	}
		$fpcCtgClone->anchor(
			id      => $currentClone,
			onclick => "closeDialog();openDialog('cloneView.cgi?cloneName=$currentClone')",
			title   => $currentClone
			)->rectangle(
				x     => $margin + $fpcCloneLeftEnd->{$currentClone} * $pixelUnit,
				y     => $barY + ($barHeight + $barSpacing) * $goodCol + $barSpacing,
				width => ($fpcCloneRightEnd->{$currentClone} - $fpcCloneLeftEnd->{$currentClone}) * $pixelUnit,
				height=> $barHeight,
				style => ($fpcCloneHighlighted->{$currentClone}) ? 
						{ stroke =>'red',
							fill => ($fpcCloneSequenced->{$currentClone}) ? 'yellow' : 'white',
							'stroke-dasharray' => ($fpcCloneMTP->{$currentClone}) ? '' : '3,1'
							} :
						{ stroke =>'grey',
							fill => ($fpcCloneSequenced->{$currentClone}) ? 'yellow' : 'white',
							'stroke-dasharray' => ($fpcCloneMTP->{$currentClone}) ? '' : '3,1'
							},
				id    => "$currentClone$$"
			);
		$scrollLeft = $fpcCloneLeftEnd->{$currentClone} * $pixelUnit if ($currentClone eq $highlight);
    	my $textX = $margin + $fpcCloneLeftEnd->{$currentClone} * $pixelUnit + 2;
    	my $textY = $barY + ($barHeight + $barSpacing) * ($goodCol + 1) - 2;
    	my $textEnd = $fpcCloneLeftEnd->{$currentClone}  + $textLength->{$currentClone} * $textFontWidth / $pixelUnit;
		$fpcCtgClone->text(
			id      => 'cloneName'.$currentClone,
			onclick => "closeDialog();openDialog('cloneView.cgi?cloneName=$currentClone')",
			x       => $textX,
			y       => $textY,
			style   => {
				'font-size'   => $textFontSize,
				'stroke'        => ($currentClone eq $highlight) ? 'red' : 'grey'
			}
		)->cdata($currentClone);
    	$colPostion[$goodCol] = ($textEnd > $fpcCloneRightEnd->{$currentClone}) ? $textEnd : $fpcCloneRightEnd->{$currentClone};
		$sequenced++ if($fpcCloneSequenced->{$currentClone});
		$colCount++;
    }

	# now render the SVG object, implicitly use svg namespace
	my $graphic = $svg->xmlify;
	my $svgWidth = $fpcCloneMaxEnd * $pixelUnit + $margin * 2;
	my $svgHeight = $barY + ($barHeight + $barSpacing) * @colPostion + $margin;
	$graphic =~ s/\$svgWidth/$svgWidth/g;
	$graphic =~ s/\$svgHeight/$svgHeight/g;
	$graphic =~ s/&/&amp;/g;

	open (SVGFILE,">$commoncfg->{TMPDIR}/FPC-$fpcCtgList[2]-$fpcCtgList[0].svg") or die "can't open file: $commoncfg->{TMPDIR}/FPC-$fpcCtgList[2]-$fpcCtgList[0].svg";
	print SVGFILE $graphic;
	close (SVGFILE);
	`gzip -f $commoncfg->{TMPDIR}/FPC-$fpcCtgList[2]-$fpcCtgList[0].svg`;

	$dialogWidth = ($svgWidth > 1000 ) ? 1050 : ($svgWidth < 550) ? 600 : $svgWidth + 50;

	$fpcCtgDetails ="
	<ul class='fpcCtgListMenu' style='width: 100px;float: left; margin-right: .3em; white-space: nowrap;'>
		<li><a><span class='ui-icon ui-icon-wrench'></span>You Can</a>
			<ul style='z-index: 1000;white-space: nowrap;'>
				<li><a onclick='openDialog(\"fpcCtgList.cgi?fpcCtgId=$fpcCtgId\");'><span class='ui-icon ui-icon-document'></span>View Clone List</a></li>
				<li><a href='$commoncfg->{TMPURL}/FPC-$fpcCtgList[2]-$fpcCtgList[0].svg.gz' target='hiddenFrame'><span class='ui-icon ui-icon-disk'></span>Download SVG</a></li>
				<li><a>Legends</a>
					<ul>
						<li>Solid bolder: MTP Clones</li>
						<li>Red bolder: Highlighted Clones</li>
						<li>Yellow filling: Sequenced Clones</li>
					</ul>
				</li>
				<li><a>Learn about</a>
					<ul>
						<li><b>$fpcCtgList[2]</b></li>
						<li>Total Clones: $colCount ($sequenced sequenced)</li>
						<li>$fpcCtgList[8]</li>
					</ul>
				</li>
			</ul>
		</li>
	</ul>
	<br><br>
	<div id='fpcCtgList$fpcCtgId$$'>
	$graphic
	</div>";

	undef $/;# enable slurp mode
	my $html = <DATA>;
	$html =~ s/\$fpcCtgDetails/$fpcCtgDetails/g;
	$html =~ s/\$dialogWidth/$dialogWidth/g;
	$html =~ s/\$fpcCtgId/$fpcCtgId/g;
	$html =~ s/\$fpcCtg/$fpcCtgList[2]/g;
	$html =~ s/\$scrollLeft/$scrollLeft/g;
	$html =~ s/\$\$/$$/g;

	print header(-cookie=>cookie(-name=>'fpcCtg',-value=>$fpcCtgId));
	print $html;
}
else
{
	print header(-type=>'text/html',-status=>'402 Invalid operation');
	exit;
}

__DATA__
$fpcCtgDetails
<script>
buttonInit();
$( ".fpcCtgListMenu" ).menu();
$('#viewer').dialog("option", "title", "View FPC $fpcCtg");
$('#viewer').dialog("option", "width", "$dialogWidth");
$('#viewer').scrollLeft($scrollLeft);
</script>