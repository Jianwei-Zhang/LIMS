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
use user;
use userCookie;

my $userCookie = new userCookie;
my $userId = (cookie('cid')) ? $userCookie->checkCookie(cookie('cid')) : 0;
exit if (!$userId);

my $sortCtgByLength = 1; #temporally use fixed number, but will change to use cookie in the future
#my $sortCtgByLength = cookie('sortCtgByLength') || '0';

my $user = new user;
my $userDetail = $user->getAllFieldsWithUserId($userId);
my $userName = $userDetail->{"userName"};

my $commoncfg = readConfig("main.conf");
my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});

my $assemblyId = param ('assemblyId') || '';
my $assemblyStatus;
$assemblyStatus->{2} = "Frozen";
$assemblyStatus->{1} = "Assembled";
$assemblyStatus->{0} = "Initialized";
$assemblyStatus->{-1} = "Assembling...";
$assemblyStatus->{-2} = "SeqToSeq Aligning...";
$assemblyStatus->{-3} = "Physical Reference-based Assembling...";
$assemblyStatus->{-4} = "SeqToGenome Aligning...";
$assemblyStatus->{-5} = "Assigning Chr#...";
$assemblyStatus->{-6} = "Filtering Seqs...";
$assemblyStatus->{-7} = "Orienting Seqs...";
$assemblyStatus->{-8} = "Orienting Contigs...";
$assemblyStatus->{-9} = "EndToEnd Merging...";
$assemblyStatus->{-10} = "Filtering Overlaps...";
$assemblyStatus->{-11} = "Estimating Length...";
$assemblyStatus->{-12} = "BesToSeq Aligning...";
my $autoSeqSearchUrl = 'autoSeqSearch.cgi';
if($assemblyId)
{
	$autoSeqSearchUrl .= "?assemblyId=$assemblyId";
}

my $assemblyList = '';
my $assemblySortableStyle = '';
my $assemblySortableJs = '';
my $lengthList;
my $totalAssembled;
my $totalLength;
my $maxLength;
my $minLength;

if ($assemblyId)
{
	my $assembly=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$assembly->execute($assemblyId);
	my @assembly = $assembly->fetchrow_array();

	my $assemblyDetails = decode_json $assembly[8];
	$assemblyDetails->{'autoCheckNewSeq'} = 0 if (!exists $assemblyDetails->{'autoCheckNewSeq'});
	$assemblyDetails->{'log'} = 'None.' if (!exists $assemblyDetails->{'log'});
	$assemblyDetails->{'log'} = 'None.' unless ($assemblyDetails->{'log'}); #give a value to $assemblyDetails->{'log'}
	$assemblyDetails->{'log'} = escapeHTML($assemblyDetails->{'log'});
	$assemblyDetails->{'log'} =~ s/\n/<br>/g;
	$assemblyDetails->{'description'} = 'None.' if (!exists $assemblyDetails->{'description'});
	$assemblyDetails->{'description'} = 'None.' unless ($assemblyDetails->{'description'}); #give a value to $assemblyDetails->{'description'}
	$assemblyDetails->{'description'} = escapeHTML($assemblyDetails->{'description'});
	$assemblyDetails->{'description'} =~ s/\n/<br>/g;

	my $target=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$target->execute($assembly[4]);
	my @target = $target->fetchrow_array();

	my $gapFillerViewer = '';
	$gapFillerViewer = "<li><a onclick='openDialog(\"$commoncfg->{TMPURL}/$assemblyId.gapfill.html\")'><span class='ui-icon ui-icon-bullet'></span>View</a></li>" if (-e "$commoncfg->{TMPDIR}/$assemblyId.gapfill.html");
	my $gapFillerDownload = '';
	$gapFillerDownload = "<li><a href='$commoncfg->{TMPURL}/$assemblyId.gapfill.txt.gz' target='hiddenFrame'><span class='ui-icon ui-icon-bullet'></span>Download</a></li>" if (-e "$commoncfg->{TMPDIR}/$assemblyId.gapfill.txt.gz");
	$assemblyList .= "<ul class='assemblyMenu' style='left: 0px;top: 0px;display:inline-block;width: 200px;'>
				<li><a><span class='ui-icon ui-icon-wrench'></span>Assembly Tools</a>
					<ul style='z-index: 1000;white-space: nowrap;'>
						<li><a onclick='openDialog(\"assemblyEdit.cgi?assemblyId=$assemblyId\")'><span class='ui-icon ui-icon-pencil'></span>Edit This Assembly</a></li>
						<li><a onclick='openDialog(\"assemblyView.cgi?assemblyId=$assemblyId\")'><span class='ui-icon ui-icon-copy'></span>Copy This Assembly</a></li>
						<li><a onclick='openDialog(\"assemblyRunForm.cgi?assemblyId=$assemblyId\")'><span class='ui-icon ui-icon-shuffle'></span>Run Assembly</a></li>
						<li><a><span class='ui-icon ui-icon-transfer-e-w'></span>Alignment</a>
							<ul style='z-index: 1000;white-space: nowrap;'>
								<li><a onclick='openDialog(\"alignmentForm.cgi?assemblyId=$assemblyId\")'><span class='ui-icon ui-icon-bullet'></span>New Alignment</a></li>
								<li><a onclick='openDialog(\"alignmentLoadForm.cgi?assemblyId=$assemblyId\")' title='BLAST Tabular Format'><span class='ui-icon ui-icon-bullet'></span>Load Alignment in Tabular</a></li>";
	$assemblyList .= ($target[1] eq 'library') ? "<li><a onclick='openDialog(\"besToSeqForm.cgi?libraryId=$assembly[4]&targetId=$assembly[5]\")'><span class='ui-icon ui-icon-bullet'></span>BES to Seq</a></li>" 
						: ($target[6] > 0) ? "<li><a onclick='openDialog(\"besToSeqForm.cgi?libraryId=$target[6]&targetId=$assembly[5]\")'><span class='ui-icon ui-icon-bullet'></span>BES to Seq</a></li>" 
						: "";
	$assemblyList .= "</ul>
						</li>
						<li><a onclick='openDialog(\"assemblyAssign.cgi?assemblyId=$assemblyId\")'><span class='ui-icon ui-icon-arrow-2-e-w'></span>Assign Ctg To Chromosome</a></li>
						<li><a onclick='openDialog(\"assemblyCtgResetForm.cgi?assemblyId=$assemblyId\")'><span class='ui-icon ui-icon-bullet'></span>Reset Redundancy</a></li>
						<li><a onclick='openDialog(\"assemblyValidation.cgi?assemblyId=$assemblyId\")'><span class='ui-icon ui-icon-search'></span>Assembly Validation</a></li>
						<li><a><span class='ui-icon ui-icon-link'></span>Gap Filler</a>
							<ul style='z-index: 1000;white-space: nowrap;'>
								<li><a onclick='openDialog(\"assemblyGapFillForm.cgi?assemblyId=$assemblyId\")'><span class='ui-icon ui-icon-bullet'></span>Run</a></li>
								$gapFillerViewer
								$gapFillerDownload
							</ul>
						</li>
						<li><a><span class='ui-icon ui-icon-disk'></span>Download</a>
							<ul style='z-index: 1000;white-space: nowrap;'>";
	$assemblyList .= ($target[1] eq 'library') ? "<li><a href='download.cgi?libraryId=$assembly[4]' target='hiddenFrame'><span class='ui-icon ui-icon-bullet'></span>BAC Sequence</a></li>
						<li><a href='download.cgi?besLibraryId=$assembly[4]' target='hiddenFrame'><span class='ui-icon ui-icon-bullet'></span>BAC End Sequence</a></li>"
						: "";

	$assemblyList .= "<li><a href='download.cgi?assemblyId=$assemblyId&chr=all' target='hiddenFrame'><span class='ui-icon ui-icon-bullet'></span>Ctg PseudoMolecules</a></li>
						<li><a href='download.cgi?assemblyId=$assemblyId&chr=all&unit=chr' target='hiddenFrame' title='100 Ns will be added to connect two contigs'><span class='ui-icon ui-icon-bullet'></span>Chr PseudoMolecules</a></li>
						<li><a href='download.cgi?assemblyIdForAgp=$assemblyId&chr=all' target='hiddenFrame'><span class='ui-icon ui-icon-bullet'></span>Ctg-Seq AGP</a></li>
						<li><a href='download.cgi?assemblyIdForAgp=$assemblyId&chr=all&unit=chr&element=ctg' target='hiddenFrame' title='100 Ns will be added to connect two contigs'><span class='ui-icon ui-icon-bullet'></span>Chr-Ctg AGP</a></li>
						<li><a href='download.cgi?assemblyIdForAgp=$assemblyId&chr=all&unit=chr' target='hiddenFrame' title='100 Ns will be added to connect sequences at edges of two contigs'><span class='ui-icon ui-icon-bullet'></span>Chr-Seq AGP</a></li>
							</ul>
						</li>
						<li><a onclick='printDiv(\"assemblyTab$assemblyId$$\")'><span class='ui-icon ui-icon-print'></span>Print</a></li>
					</ul>
				</li>
			</ul>";
			
	$assemblyList .= ($assembly[7] < 1) ? "<button style='float: right;' onclick='refresh(\"menu\")' title='Click to Refresh'>Status: $assemblyStatus->{$assembly[7]}</button>" : 
							($assembly[7] == 1) ? "<button style='float: right;' onclick='loaddiv(\"hiddenDiv\", \"assemblyFreeze.cgi?assemblyId=$assemblyId\")' title='Click to Freeze'>$assemblyStatus->{$assembly[7]}</button>
							<button style='float: right;' onclick='refresh(\"menu\")'>Refresh</button>" :
							"<button style='float: right;' onclick='loaddiv(\"hiddenDiv\", \"assemblyFreeze.cgi?assemblyId=$assemblyId\")' title='Click to Unfreeze'>$assemblyStatus->{$assembly[7]}</button>
							<button style='float: right;' onclick='refresh(\"menu\")'>Refresh</button>";
	$assemblyList .= "<input style='float: right;'  class='ui-widget-content ui-corner-all' name='seqName' id='assemblySearchSeqName$$' size='16' type='text' maxlength='32' VALUE='' placeholder='Search Seq' />
			<hr>";
	my $assemblySeq;
	my $assemblyCtg;
	my $assemblyCtgNumber = 0;
	my $assemblyCtgChr;
	my $assemblyCtgChrOrder;
	my $assemblyCtgLength;
	my $assemblyCtgSeqNumber;
	my $inCtgSeq;
	my $assemblyCtgs = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblyCtg' AND o = ?");
	$assemblyCtgs->execute($assembly[0]);
	while(my @assemblyCtgs = $assemblyCtgs->fetchrow_array())
	{
		$assemblyCtg->{$assemblyCtgs[0]} = $assemblyCtgs[2];
		$assemblyCtgNumber = $assemblyCtgs[2] if($assemblyCtgs[2] > $assemblyCtgNumber);
		$assemblyCtgChr->{$assemblyCtgs[0]} = $assemblyCtgs[4];
		$assemblyCtgChrOrder->{$assemblyCtgs[0]} = ($assemblyCtgs[4]) ? $assemblyCtgs[5] : ($sortCtgByLength) ? "-$assemblyCtgs[7]" : $assemblyCtgs[2];
		$assemblyCtgLength->{$assemblyCtgs[0]} = $assemblyCtgs[7];

		unless (exists $totalLength->{'All'})
		{
			$totalAssembled->{'All'} = 0;
			$totalLength->{'All'} = 0;
			$maxLength->{'All'} = 0;
			$minLength->{'All'} = 999999999;
		}
		push @{$lengthList->{'All'}},$assemblyCtgs[7];
		$totalAssembled->{'All'}++;
		$totalLength->{'All'} += $assemblyCtgs[7];
		$maxLength->{'All'} = $assemblyCtgs[7] if ($assemblyCtgs[7] > $maxLength->{'All'});
		$minLength->{'All'} = $assemblyCtgs[7] if ($assemblyCtgs[7] < $minLength->{'All'});

		if($assemblyCtgs[4])
		{
			unless (exists $totalLength->{'All Anchored'})
			{
				$totalAssembled->{'All Anchored'} = 0;
				$totalLength->{'All Anchored'} = 0;
				$maxLength->{'All Anchored'} = 0;
				$minLength->{'All Anchored'} = 999999999;
			}
			push @{$lengthList->{'All Anchored'}},$assemblyCtgs[7];
			$totalAssembled->{'All Anchored'}++;
			$totalLength->{'All Anchored'} += $assemblyCtgs[7];
			$maxLength->{'All Anchored'} = $assemblyCtgs[7] if ($assemblyCtgs[7] > $maxLength->{'All Anchored'});
			$minLength->{'All Anchored'} = $assemblyCtgs[7] if ($assemblyCtgs[7] < $minLength->{'All Anchored'});

			if ($assemblyCtgs[4] % 100 == 98)
			{
				unless (exists $totalLength->{'Chloroplast'})
				{
					$totalAssembled->{'Chloroplast'} = 0;
					$totalLength->{'Chloroplast'} = 0;
					$maxLength->{'Chloroplast'} = 0;
					$minLength->{'Chloroplast'} = 999999999;
				}
				push @{$lengthList->{'Chloroplast'}},$assemblyCtgs[7];
				$totalAssembled->{'Chloroplast'}++;
				$totalLength->{'Chloroplast'} += $assemblyCtgs[7];
				$maxLength->{'Chloroplast'} = $assemblyCtgs[7] if ($assemblyCtgs[7] > $maxLength->{'Chloroplast'});
				$minLength->{'Chloroplast'} = $assemblyCtgs[7] if ($assemblyCtgs[7] < $minLength->{'Chloroplast'});
			}
			elsif($assemblyCtgs[4] % 100 == 99)
			{
				unless (exists $totalLength->{'Mitochondrion'})
				{
					$totalAssembled->{'Mitochondrion'} = 0;
					$totalLength->{'Mitochondrion'} = 0;
					$maxLength->{'Mitochondrion'} = 0;
					$minLength->{'Mitochondrion'} = 999999999;
				}
				push @{$lengthList->{'Mitochondrion'}},$assemblyCtgs[7];
				$totalAssembled->{'Mitochondrion'}++;
				$totalLength->{'Mitochondrion'} += $assemblyCtgs[7];
				$maxLength->{'Mitochondrion'} = $assemblyCtgs[7] if ($assemblyCtgs[7] > $maxLength->{'Mitochondrion'});
				$minLength->{'Mitochondrion'} = $assemblyCtgs[7] if ($assemblyCtgs[7] < $minLength->{'Mitochondrion'});
			}
			elsif($assemblyCtgs[4] % 100 == 0)
			{
				unless (exists $totalLength->{'Contamination'})
				{
					$totalAssembled->{'Contamination'} = 0;
					$totalLength->{'Contamination'} = 0;
					$maxLength->{'Contamination'} = 0;
					$minLength->{'Contamination'} = 999999999;
				}
				push @{$lengthList->{'Contamination'}},$assemblyCtgs[7];
				$totalAssembled->{'Contamination'}++;
				$totalLength->{'Contamination'} += $assemblyCtgs[7];
				$maxLength->{'Contamination'} = $assemblyCtgs[7] if ($assemblyCtgs[7] > $maxLength->{'Contamination'});
				$minLength->{'Contamination'} = $assemblyCtgs[7] if ($assemblyCtgs[7] < $minLength->{'Contamination'});
			}
			else
			{
				unless (exists $totalLength->{'All Chromosome'})
				{
					$totalAssembled->{'All Chromosome'} = 0;
					$totalLength->{'All Chromosome'} = 0;
					$maxLength->{'All Chromosome'} = 0;
					$minLength->{'All Chromosome'} = 999999999;
				}
				push @{$lengthList->{'All Chromosome'}},$assemblyCtgs[7];
				$totalAssembled->{'All Chromosome'}++;
				$totalLength->{'All Chromosome'} += $assemblyCtgs[7];
				$maxLength->{'All Chromosome'} = $assemblyCtgs[7] if ($assemblyCtgs[7] > $maxLength->{'All Chromosome'});
				$minLength->{'All Chromosome'} = $assemblyCtgs[7] if ($assemblyCtgs[7] < $minLength->{'All Chromosome'});
				if($assemblyCtgs[4] > 100)
				{
					my $subGenomeNumber = substr ($assemblyCtgs[4], 0, -2);
					my $subGenomeName = "Subgenome-$subGenomeNumber Chromosome";

					unless (exists $totalLength->{$subGenomeName})
					{
						$totalAssembled->{$subGenomeName} = 0;
						$totalLength->{$subGenomeName} = 0;
						$maxLength->{$subGenomeName} = 0;
						$minLength->{$subGenomeName} = 999999999;
					}
					push @{$lengthList->{$subGenomeName}},$assemblyCtgs[7];
					$totalAssembled->{$subGenomeName}++;
					$totalLength->{$subGenomeName} += $assemblyCtgs[7];
					$maxLength->{$subGenomeName} = $assemblyCtgs[7] if ($assemblyCtgs[7] > $maxLength->{$subGenomeName});
					$minLength->{$subGenomeName} = $assemblyCtgs[7] if ($assemblyCtgs[7] < $minLength->{$subGenomeName});
				}
			}
		}
		else
		{
			unless (exists $totalLength->{'Unplaced'})
			{
				$totalAssembled->{'Unplaced'} = 0;
				$totalLength->{'Unplaced'} = 0;
				$maxLength->{'Unplaced'} = 0;
				$minLength->{'Unplaced'} = 999999999;
			}
			push @{$lengthList->{'Unplaced'}},$assemblyCtgs[7];
			$totalAssembled->{'Unplaced'}++;
			$totalLength->{'Unplaced'} += $assemblyCtgs[7];
			$maxLength->{'Unplaced'} = $assemblyCtgs[7] if ($assemblyCtgs[7] > $maxLength->{'Unplaced'});
			$minLength->{'Unplaced'} = $assemblyCtgs[7] if ($assemblyCtgs[7] < $minLength->{'Unplaced'});
		}
		my @seqs = split ",", $assemblyCtgs[8];
		$assemblyCtgSeqNumber->{$assemblyCtgs[0]} = @seqs;
		for (@seqs)
		{
			next unless ($_);
			$_ =~ s/[^a-zA-Z0-9]//g;
			$inCtgSeq->{$_} = 1;
		}
	}
	my $totalInAssembly = keys %$inCtgSeq;

	my $assemblySeqs = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblySeq' AND o = ?");
	$assemblySeqs->execute($assembly[0]);
	while(my @assemblySeqs = $assemblySeqs->fetchrow_array())
	{
		$assemblySeq->{$assemblySeqs[5]} = $assemblySeqs[4];
	}

	my $margin = 20;
	my $chartWidth = 600;
	my $chartHeight = 400;
	my $svgWidth = $chartWidth + $margin*2;
	my $svgHeight = $chartHeight + $margin*2;
	my $assemblyStats = '';
	foreach my $contigType (sort keys %$totalAssembled)
	{
		@{$lengthList->{$contigType}} = sort {$b <=> $a} @{$lengthList->{$contigType}};
		my $meanLength = int ($totalLength->{$contigType}/$totalAssembled->{$contigType});
		my $medianLength = 0;
		if($totalAssembled->{$contigType} % 2 == 0)
		{
			my $median = $totalAssembled->{$contigType}/2;
			$medianLength = ($lengthList->{$contigType}[$median]+$lengthList->{$contigType}[$median-1])/2;
		}
		else
		{
			my $median = int ($totalAssembled->{$contigType}/2);
			$medianLength = $lengthList->{$contigType}[$median];
		}
		my $n50Length = 0;
		my $subtotal = 0;
		my $widthUnit = int ($maxLength->{$contigType}/ $chartWidth) || 1; 
		my $heightUnit = ($totalAssembled->{$contigType} > 0) ? $chartHeight/$totalAssembled->{$contigType} : 0; 
		my $assemblyCtgLengthCount;
		my %countByLength = ('halfK'  => '0',
			'oneK'    => '0',
			'tenK'   => '0',
			'hundredK'  => '0',
			'oneM' => '0');
		my %sumLength = ('halfK'  => '0',
			'oneK'    => '0',
			'tenK'   => '0',
			'hundredK'  => '0',
			'oneM' => '0');
		my $lFifty=0;
		my $lFiftyCounter=0;
		foreach (@{$lengthList->{$contigType}})
		{
			$subtotal += $_;
			$lFiftyCounter++;
			if($subtotal > $totalLength->{$contigType}/2 && $n50Length == 0)
			{
				$n50Length = $_;
				$lFifty = $lFiftyCounter;
			}
			for (my $i = 1; $i <= $_; $i += $widthUnit)
			{
				$assemblyCtgLengthCount->{$i} = 0 unless (exists $assemblyCtgLengthCount->{$i});
				$assemblyCtgLengthCount->{$i}++;
			}
			$countByLength{'halfK'}++ if ($_ > 500);
			$countByLength{'oneK'}++ if ($_ > 1000);
			$countByLength{'tenK'}++ if ($_ > 10000);
			$countByLength{'hundredK'}++ if ($_ > 100000);
			$countByLength{'oneM'}++ if ($_ > 1000000);
			$sumLength{'halfK'} += $_ if ($_ > 500);
			$sumLength{'oneK'} += $_ if ($_ > 1000);
			$sumLength{'tenK'} += $_ if ($_ > 10000);
			$sumLength{'hundredK'} += $_ if ($_ > 100000);
			$sumLength{'oneM'} += $_ if ($_ > 1000000);
		}
		my $graphic = '';
		if ($totalLength->{$contigType} > 0)
		{
			# create an SVG object
			my $svg= SVG->new(width=>$svgWidth,height=>$svgHeight); # set width and height after knowing the size
			# use explicit element constructor to generate a group element
			my $length=$svg->group(
				id    => 'length'.$contigType,
				style => { stroke=>'black',
					fill =>'white'
				}
			);
			$length->rectangle(
					x=> $margin, y=> $margin,
					width=>$chartWidth, height=>$chartHeight,
					id=>'rect_1'
				);
			my $stringOfContigNumber = ($totalAssembled->{$contigType} > 0) ? "Total ". commify($totalAssembled->{$contigType}) ." contigs" : "Total ". commify($totalAssembled->{$contigType}) ." contig";
			my $lengthOfContigNumber = length ($stringOfContigNumber);
			my $xPositionOfContigNumber = $margin - 5;
			my $yPositionOfContigNumber = $lengthOfContigNumber * 7;
			$length->text(
					id      => 'totalContigs'.$contigType,
					x       => $xPositionOfContigNumber,
					y       => $yPositionOfContigNumber,
					style=>{
						stroke=>'black',
					fill =>'black'
					},
					transform => "rotate(-90,$xPositionOfContigNumber,$yPositionOfContigNumber)"
				)->cdata($stringOfContigNumber);

			$length->text(
					id      => 'minLength'.$contigType,
					x       => 0,
					y       => $svgHeight - 5,
					style=>{
						stroke=>'black',
					fill =>'black'
					}
				)->cdata("Shortest: ". commify($minLength->{$contigType}) . " bp");

			my $lengthOfString = length ("Longest: ". commify($maxLength->{$contigType}) . " bp");
			$length->text(
					id      => 'maxLength'.$contigType,
					x       => $svgWidth - $lengthOfString * 8,
					y       => $svgHeight - 5,
					style=>{
						stroke=>'black',
					fill =>'black'
					}
				)->cdata("Longest: ". commify($maxLength->{$contigType}) . " bp");

			my $n50y = 0;
			my $mediany = 0;
			my @xvArray;
			my @yvArray;
			for (sort {$a <=> $b} keys %$assemblyCtgLengthCount)
			{
				push @xvArray, $_ / $widthUnit + $margin;
				push @yvArray, $chartHeight - $heightUnit * $assemblyCtgLengthCount->{$_}  + $margin;
				$n50y  = $chartHeight - $heightUnit * $assemblyCtgLengthCount->{$_}  + $margin if( $n50Length > $_);
				$mediany  = $chartHeight - $heightUnit * $assemblyCtgLengthCount->{$_} + $margin if( $medianLength > $_);
			}
			my $xv = [@xvArray];
			my $yv = [@yvArray];
			my $points = $length->get_path(
				x=>$xv, y=>$yv,
				-type=>'polyline',
				-closed=>'false' #specify that the polyline is closed.
			);
			$length->polyline (
				%$points,
				id=>'pline'.$contigType,
				style=>{
					'fill-opacity'=>0,
					'stroke-color'=>'rgb(250,123,23)'
				}
			);

			$length->line(
					id=>'mediana'.$contigType,
					style => {
						stroke=> 'lightgrey',
						'stroke-dasharray' => '3,3'
						},
					x1=>$medianLength / $widthUnit  + $margin, y1=>$mediany,
					x2=>$medianLength / $widthUnit  + $margin, y2=>$chartHeight + $margin
				);
	# 		$length->line(
	# 				id=>'medianb'.$contigType,
	# 				x1=> $margin, y1=>$mediany,
	# 				x2=>$medianLength / $widthUnit + $margin, y2=>$mediany
	# 			);
			$length->text(
					id      => 'mediantext'.$contigType,
					x       => $medianLength / $widthUnit  + $margin,
					y       => $mediany,
					style=>{
						stroke=>'black',
					fill =>'black'
					}
				)->cdata("Median: ". commify($medianLength) . " bp");

			$length->line(
					id=>'n50a'.$contigType,
					style => {
						stroke=> 'lightgrey',
						'stroke-dasharray' => '3,3'
						},
					x1=>$n50Length / $widthUnit + $margin, y1=>$n50y,
					x2=>$n50Length / $widthUnit + $margin, y2=>$chartHeight + $margin
				);
	# 		$length->line(
	# 				id=>'n50b'.$contigType,
	# 				x1=> $margin, y1=>$n50y,
	# 				x2=>$n50Length / $widthUnit + $margin, y2=>$n50y
	# 			);
			$length->text(
					id      => 'n50text'.$contigType,
					x       => $n50Length / $widthUnit + $margin,
					y       => $n50y,
					style=>{
						stroke=>'black',
					fill =>'black'
					}
				)->cdata("N50: ". commify($n50Length) . " bp");

			# now render the SVG object, implicitly use svg namespace
			$graphic = $svg->xmlify;
		}

		if ($totalAssembled->{$contigType} > 0)
		{
			$assemblyStats .= "<hr style='page-break-before: always;'><div><h3>$contigType Contigs</h3>
			<div style='display:inline-block;'>$graphic</div>
			<div style='display:inline-block;'>Number of sequences: " . commify($totalAssembled->{$contigType}). ".<br>
			Total size of sequences: " . commify($totalLength->{$contigType}). " bp.<br>
			Longest sequence: " . commify($maxLength->{$contigType}). " bp.<br>
			Shortest sequence: " . commify($minLength->{$contigType}). " bp.<br>
			Number of sequences > 500 nt: $countByLength{'halfK'} (" . commify($sumLength{'halfK'}) ." bp).<br>
			Number of sequences > 1k nt: $countByLength{'oneK'} (" . commify($sumLength{'oneK'}) ." bp).<br>
			Number of sequences > 10k nt: $countByLength{'tenK'} (" . commify($sumLength{'tenK'}) ." bp).<br>
			Number of sequences > 100k nt: $countByLength{'hundredK'} (" . commify($sumLength{'hundredK'}) ." bp).<br>
			Number of sequences > 1M nt: $countByLength{'oneM'} (" . commify($sumLength{'oneM'}) ." bp).<br>
			Mean sequence size: " . commify($meanLength). " bp.<br>
			Median sequence size: " . commify($medianLength). " bp.<br>
			N50 sequence length: " . commify($n50Length). " bp.<br>
			L50 sequence count: " . commify($lFifty). ".</div>
			</div>";
		}
	}

	my @assemblyCtg = sort { $assemblyCtgChrOrder->{$a} <=> $assemblyCtgChrOrder->{$b} } keys %$assemblyCtgChrOrder;
	my $assembledCtgByChr;
	my $totalAssembledContigByChr;
	my $totalAssembledSeqNumberByChr;
	my $assemblySortableJsByChr;
	my $assemblySortableStyleByChr;
	my $assemblyCtgOrders;
    for (@assemblyCtg)
    {
		my %commentType = (
			0=>"<a style='float: right;' onclick='closeDialog();openDialog(\"comment.cgi?itemId=$_\")'><span style=' left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-comment' title='Comments'></span></a>",
			1=>"<a style='float: right;' class='ui-state-error-text' onclick='closeDialog();openDialog(\"comment.cgi?itemId=$_\")'><span style=' left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-alert' title='Comments'></span></a>",
			2=>"<a style='float: right;' class='ui-state-error-text' onclick='closeDialog();openDialog(\"comment.cgi?itemId=$_\")'><span style=' left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-heart' title='Comments'></span></a>"
			);
		my $comment = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'comment' AND o = ?");
		$comment->execute($_);
		my $commentIcon;
		if ($comment->rows > 0)
		{
			my @comment = $comment->fetchrow_array();
			$commentIcon = $commentType{$comment[4]};
		}
		else
		{
			$commentIcon = "<a style='float: right;' class='ui-state-disabled' onclick='closeDialog();openDialog(\"comment.cgi?itemId=$_\")'><span style=' left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-comment' title='Add Comments'></span></a>";
		}
    	$assembledCtgByChr->{$assemblyCtgChr->{$_}} = "<form id='assemblyChr$assemblyId$assemblyCtgChr->{$_}' name='assemblyChr$assemblyId$assemblyCtgChr->{$_}' action='assemblyCtgOrder.cgi' enctype='multipart/form-data' method='post' target='hiddenFrame'><input name='assemblyId' id='assemblyId' type='hidden' value='$assemblyId' /><ul id='sortableCtgList$assemblyId$assemblyCtgChr->{$_}'>" unless (exists $assembledCtgByChr->{$assemblyCtgChr->{$_}});
		$assembledCtgByChr->{$assemblyCtgChr->{$_}} .= ($assemblyCtgSeqNumber->{$_} > 1) ?
				"<li class='ui-state-default' id='$_'>$commentIcon<a onclick='closeViewer();openViewer(\"assemblyCtgView.cgi?assemblyCtgId=$_\")'>Ctg$assemblyCtg->{$_} ($assemblyCtgSeqNumber->{$_})</a><br>" . commify($assemblyCtgLength->{$_}). " bp</li>"
				: "<li class='ui-state-default' id='$_'>$commentIcon<a onclick='closeViewer();openViewer(\"assemblyCtgView.cgi?assemblyCtgId=$_\")'>Ctg$assemblyCtg->{$_}</a><br>" . commify($assemblyCtgLength->{$_}). " bp</li>";
		
		$totalLength->{$assemblyCtgChr->{$_}} = 0 unless (exists $totalLength->{$assemblyCtgChr->{$_}});
		$totalLength->{$assemblyCtgChr->{$_}} += $assemblyCtgLength->{$_};
		$maxLength->{$assemblyCtgChr->{$_}} = 0 unless (exists $maxLength->{$assemblyCtgChr->{$_}});
		$maxLength->{$assemblyCtgChr->{$_}} = $assemblyCtgLength->{$_} if ($assemblyCtgLength->{$_} > $maxLength->{$assemblyCtgChr->{$_}});
		$minLength->{$assemblyCtgChr->{$_}} = 999999 unless (exists $minLength->{$assemblyCtgChr->{$_}});
		$minLength->{$assemblyCtgChr->{$_}} = $assemblyCtgLength->{$_} if ($assemblyCtgLength->{$_} < $minLength->{$assemblyCtgChr->{$_}});
    	$totalAssembledContigByChr->{$assemblyCtgChr->{$_}} = 0 unless (exists $totalAssembledContigByChr->{$assemblyCtgChr->{$_}});
		$totalAssembledSeqNumberByChr->{$assemblyCtgChr->{$_}} = 0 unless (exists $totalAssembledSeqNumberByChr->{$assemblyCtgChr->{$_}});
		$totalAssembledContigByChr->{$assemblyCtgChr->{$_}}++;
		$totalAssembledSeqNumberByChr->{$assemblyCtgChr->{$_}} += $assemblyCtgSeqNumber->{$_};

    	$assemblyCtgOrders->{$assemblyCtgChr->{$_}} = '' unless (exists $assemblyCtgOrders->{$assemblyCtgChr->{$_}});
    	$assemblyCtgOrders->{$assemblyCtgChr->{$_}} .= ($assemblyCtgOrders->{$assemblyCtgChr->{$_}}) ? ",$_" : $_;

    	$assemblySortableJsByChr->{$assemblyCtgChr->{$_}} = <<END;
			\$( "#sortableCtgList$assemblyId$assemblyCtgChr->{$_}" ).sortable({
				  placeholder: "ui-state-highlight",
				  forcePlaceholderSize: true,
				  stop: function(event, ui) {
					var oldSortedIds = \$('#assemblyCtgOrders$assemblyId$assemblyCtgChr->{$_}').val();
					var newSortedIds = \$( "#sortableCtgList$assemblyId$assemblyCtgChr->{$_}" ).sortable( "toArray" );
					if(oldSortedIds != newSortedIds)
					{
						savingShow();
						\$('#assemblyCtgOrders$assemblyId$assemblyCtgChr->{$_}').val(newSortedIds);
						submitForm('assemblyChr$assemblyId$assemblyCtgChr->{$_}');
					}
				}
			});
			\$( "#sortableCtgList$assemblyId$assemblyCtgChr->{$_}" ).disableSelection();
END
    	$assemblySortableStyleByChr->{$assemblyCtgChr->{$_}} = <<END;
		  <style>
		  #sortableCtgList$assemblyId$assemblyCtgChr->{$_} { list-style-type: none; display:inline-block;margin: 0; padding: 0; width: 100%; }
		  #sortableCtgList$assemblyId$assemblyCtgChr->{$_} li { margin: 3px 3px 3px 0; padding: 1px; float: left; width: 120px; text-align: center; }
		  </style>
END
    }
    my $assembledCtgDetails = '';
	my $assemblyStatsSummary = "";
    for (sort {($b % 100 == $a % 100) ? ($b <=> $a) : ($b % 100 <=> $a % 100) } keys %$assembledCtgByChr)
    {
		$assembledCtgByChr->{$_} .= "</ul><input name='assemblyCtgOrders' id='assemblyCtgOrders$assemblyId$_' type='hidden' value='$assemblyCtgOrders->{$_}' /></form>";
		my $headerByChr;
		my $formattedChr = ($_ == 0) ? "Unplaced"
			: ($_ % 100 == 0) ? "Contamination" :
				($_ % 100 == 98) ? "Chloroplast" : 
					($_ % 100 == 99) ? "Mitochondrion" : 
						($_ > 100) ? "Subgenome-" . substr ($_, 0, -2) . " Chromosome " . substr ($_, -2) : "Chromosome $_";
		$headerByChr->{$_} = ($_ % 100 > 0 && $assembly[5] > 0) ?
			"<h3><a onclick='closeViewer();openViewer(\"assemblyChrView.cgi?assemblyId=$assemblyId&chr=$_\")'>$formattedChr</a>"
			: "<h3><a class='ui-state-error-text'>$formattedChr</a>";

		$headerByChr->{$_} .= ($totalAssembledContigByChr->{$_} > 1) ?
			" ($totalAssembledContigByChr->{$_} Contigs)"
			: " ($totalAssembledContigByChr->{$_} Contig)";

		$headerByChr->{$_} .= ($totalAssembledSeqNumberByChr->{$_} > 1) ?
			" - $totalAssembledSeqNumberByChr->{$_} Seqs"
			: " - $totalAssembledSeqNumberByChr->{$_} Seq";

		$headerByChr->{$_} .= ($totalAssembledContigByChr->{$_} > 1) ?
			" <sup><a title = 'total length with Ns between assemblySeqs but excluding Ns (like for artificial telomeres) at Ctg ends'>" . commify($totalLength->{$_}). "</a> bp <ul class='assemblyMenu' style='width: 150px;display:inline-block; white-space: nowrap;'><li><a><span class='ui-icon ui-icon-disk'></span>Download</a><ul style='z-index: 1000;white-space: nowrap;'><li><a href='download.cgi?assemblyId=$assemblyId&chr=$_' target='hiddenFrame' title='Click to Download'><span class='ui-icon ui-icon-bullet'></span>Ctg Pseudomolecules</a></li>"
			: " <sup><a title = 'total length with Ns between assemblySeqs but excluding Ns (like for artificial telomeres) at Ctg ends'>" . commify($totalLength->{$_}). "</a> bp <ul class='assemblyMenu' style='width: 150px;display:inline-block; white-space: nowrap;'><li><a><span class='ui-icon ui-icon-disk'></span>Download</a><ul style='z-index: 1000;white-space: nowrap;'><li><a href='download.cgi?assemblyId=$assemblyId&chr=$_' target='hiddenFrame' title='Click to Download'><span class='ui-icon ui-icon-bullet'></span>Ctg Pseudomolecule</a></li>";

		$headerByChr->{$_} .= ($_ > 0) ?
			($_ % 100 == 0) ? "<li><a href='download.cgi?assemblyId=$assemblyId&chr=$_&unit=chr' target='hiddenFrame'><span class='ui-icon ui-icon-bullet'></span>Chr Pseudomolecule</a></li><li><a href='download.cgi?assemblyIdForAgp=$assemblyId&chr=$_' target='hiddenFrame' title='Click to Download'><span class='ui-icon ui-icon-bullet'></span>Ctg-Seq AGP</a></li><li><a href='download.cgi?assemblyIdForAgp=$assemblyId&chr=$_&unit=chr&element=ctg' target='hiddenFrame'><span class='ui-icon ui-icon-bullet'></span>Chr-Ctg AGP</a></li></ul></li><li><a class='ui-state-error-text' onclick='deleteItem($assemblyId,\"contamination\")' title='Delete Contamination Contigs'><span class='ui-icon ui-icon-trash'></span>Delete Ctgs</a></li></ul></sup></h3>"
			: "<li><a href='download.cgi?assemblyId=$assemblyId&chr=$_&unit=chr' target='hiddenFrame' title='100 Ns will be added to connect two contigs'><span class='ui-icon ui-icon-bullet'></span>Chr Pseudomolecule</a></li><li><a href='download.cgi?assemblyIdForAgp=$assemblyId&chr=$_' target='hiddenFrame' title='Click to Download'><span class='ui-icon ui-icon-bullet'></span>Ctg-Seq AGP</a></li><li><a href='download.cgi?assemblyIdForAgp=$assemblyId&chr=$_&unit=chr&element=ctg' target='hiddenFrame' title='100 Ns will be added to connect two contigs'><span class='ui-icon ui-icon-bullet'></span>Chr-Ctg AGP</a></li><li><a href='download.cgi?assemblyIdForAgp=$assemblyId&chr=$_&unit=chr' target='hiddenFrame' title='100 Ns will be added to connect sequences at edges of two contigs'><span class='ui-icon ui-icon-bullet'></span>Chr-Seq AGP</a></li></ul></li></ul></sup></h3>"
			: "<li><a href='download.cgi?assemblyId=$assemblyId&chr=$_&unit=chr' target='hiddenFrame'><span class='ui-icon ui-icon-bullet'></span>Chr Pseudomolecule</a></li><li><a href='download.cgi?assemblyIdForAgp=$assemblyId&chr=$_' target='hiddenFrame' title='Click to Download'><span class='ui-icon ui-icon-bullet'></span>Ctg-Seq AGP</a></li><li><a href='download.cgi?assemblyIdForAgp=$assemblyId&chr=$_&unit=chr&element=ctg' target='hiddenFrame'><span class='ui-icon ui-icon-bullet'></span>Chr-Ctg AGP</a></li></ul></li><li><a class='ui-state-error-text' onclick='deleteItem($assemblyId,\"unplaced\")' title='Delete Unplaced Contigs'><span class='ui-icon ui-icon-trash'></span>Delete Ctgs</a></li></ul></sup></h3>";
		$assemblyStatsSummary = ($_ % 100 > 0) ? "<tr><td>$formattedChr</td><td>$totalAssembledContigByChr->{$_}</td><td>" . commify($totalLength->{$_}). "</td><td>$totalAssembledSeqNumberByChr->{$_}</td></tr>" . $assemblyStatsSummary
			: $assemblyStatsSummary . "<tr><td>$formattedChr</td><td>$totalAssembledContigByChr->{$_}</td><td>" . commify($totalLength->{$_}). "</td><td>$totalAssembledSeqNumberByChr->{$_}</td></tr>";

		$assembledCtgDetails = ($_ % 100 > 0) ? $headerByChr->{$_}.$assembledCtgByChr->{$_}.$assembledCtgDetails : $assembledCtgDetails.$headerByChr->{$_}.$assembledCtgByChr->{$_}; #list chromosome-assigned first.
		$assemblySortableStyle .= $assemblySortableStyleByChr->{$_};
		$assemblySortableJs .= $assemblySortableJsByChr->{$_};
    }
	$assemblyStatsSummary = "<table id='assemblyStatsSummary$$' class='display'><thead><tr><th><b>Assignment</b></th><th><b>Number of Contigs</b></th><th><b>Length (bp)</b></th><th><b>Component Seqs</b></th></tr></thead><tbody>$assemblyStatsSummary</tbody><tfoot><tr><td><b>Total</b></td><td>$totalAssembled->{'All'}</td><td>" . commify($totalLength->{'All'}). "</td><td>$totalInAssembly</td></tr></tfoot></table>";

	$assemblyList .= "<div id='assemblyTab$assemblyId$$'><ul>";
	$assemblyList .= ($totalAssembled->{'All'} > 1) ? "<li><a href='#tabs-assembled$assemblyId$$'>Assembly ($totalAssembled->{'All'} Contigs)" : "<li><a href='#tabs-assembled$assemblyId$$'>Assembly ($totalAssembled->{'All'} Contig)";
	$assemblyList .= ($totalInAssembly > 1) ? " - $totalInAssembly Sequences</a></li>" : " - $totalInAssembly Sequence</a></li>";
	$assemblyList .= ($totalAssembled->{'All'} > 0) ? "<li><a href='assemblyCtgList.cgi?assemblyId=$assemblyId'>Contig List</a></li>" : "";
	$assemblyList .= ($totalAssembled->{'All'} > 0) ? "<li><a href='#tabs-assemblyStats$assemblyId$$'>Stats</a></li>" : "";

	if($assemblyDetails->{'autoCheckNewSeq'})
	{
		my $totalSequenced = 0;
		if($target[1] eq 'library')
		{
			my $getClones = $dbh->prepare("SELECT * FROM clones WHERE libraryId = ? AND sequenced > 0");
			$getClones->execute($assembly[4]);
			while(my @getClones= $getClones->fetchrow_array())
			{
				my $getSeqs = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o < 50 AND name LIKE ?");
				$getSeqs->execute($getClones[1]);
				while(my @getSeqs = $getSeqs->fetchrow_array())
				{
					next if (exists $assemblySeq->{$getSeqs[0]});
					$totalSequenced++;
				}
				last if ($totalSequenced > 1);
			}
		}
		if($target[1] eq 'genome')
		{
			my $getSeqs = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 99 AND x = ?");
			$getSeqs->execute($assembly[4]);
			while(my @getSeqs = $getSeqs->fetchrow_array())
			{
				next if (exists $assemblySeq->{$getSeqs[0]});
				$totalSequenced++;
				last if ($totalSequenced > 1);
			}
		}
		if($totalSequenced)
		{
			$assemblyList .= ($totalSequenced > 1) ? "<li><a href='assemblyCheckSeq.cgi?assemblyId=$assemblyId'>New Sequences Available</a></li>" : "<li><a href='assemblyCheckSeq.cgi?assemblyId=$assemblyId'>New Sequence Available</a></li>";
		}
	}
	else
	{
		$assemblyList .= "<li><a href='assemblyCheckSeq.cgi?assemblyId=$assemblyId'>Check New Sequences</a></li>";
	}

	$assemblyList .= "<li><a href='#tabs-assemblyAbout$assemblyId$$'>About</a></li>";

	my $assembledCtg = ($totalAssembled->{'All'} > 0) ? "<div id='tabs-assembled$assemblyId$$'>$assembledCtgDetails</div><div id='tabs-assemblyStats$assemblyId$$'>$assemblyStatsSummary$assemblyStats</div>"
			: "<div id='tabs-assembled$assemblyId$$'>No assembly! Please run assembly!</div>";

	my $assemblyStatsSummaryJs = '';
	if ($totalAssembled->{'All'} > 0)
	{
		$assemblyStatsSummaryJs = <<END;
\$( "#assemblyStatsSummary$$" ).dataTable({
	"dom": 'lfrtip',
	"scrollCollapse": true,
	"paging": false,
	"searching": false,
	"ordering": false,
	"info": false
});
END
	}

	my $fpcOrAgpId = '';
	if($assembly[6])
	{
		my $fpcOrAgpList=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$fpcOrAgpList->execute($assembly[6]);
		my @fpcOrAgpList = $fpcOrAgpList->fetchrow_array();
		if($fpcOrAgpList[1] eq 'fpc')
		{
			$fpcOrAgpId .= "<tr><td style='text-align:right'><b>Physical Reference</b></td><td title='$fpcOrAgpList[8]'>FPC: $fpcOrAgpList[2] v.$fpcOrAgpList[3]</td></tr>";
		}
		elsif($fpcOrAgpList[1]  eq 'agp')
		{
			$fpcOrAgpId .= "<tr><td style='text-align:right'><b>Physical Reference</b></td><td>AGP: $fpcOrAgpList[2] v.$fpcOrAgpList[3]</td></tr>";
		}
	}

	my $refGenomeId = '';
	if($assembly[5])
	{
		my $genomeList=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
		$genomeList->execute($assembly[5]);
		my @genomeList = $genomeList->fetchrow_array();
		$refGenomeId = "<tr><td style='text-align:right'><b>Reference Genome</b></td><td title='$genomeList[8]'><a onclick='openDialog(\"genomeView.cgi?genomeId=$genomeList[0]\")' title='View'>$genomeList[2]<span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-arrow-1-ne' title='View'></span></a></td></tr>";
	}

	my $assemblyExtraIds = "<tr><td style='text-align:right'><b>Extra Genome</b></td><td>";
	my $checkAsbGenome = $dbh->prepare("SELECT * FROM matrix,link WHERE link.type LIKE 'asbGenome' AND link.child = matrix.id AND link.parent = ?");
	$checkAsbGenome->execute($assemblyId);
	if ($checkAsbGenome->rows > 0)
	{
		while(my @checkAsbGenome=$checkAsbGenome->fetchrow_array())
		{
			$assemblyExtraIds .= "<a onclick='closeDialog();openDialog(\"genomeView.cgi?genomeId=$checkAsbGenome[0]\")' title='$checkAsbGenome[8]'>$checkAsbGenome[2]<span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-arrow-1-ne'></span></a><br>";
		}
		$assemblyExtraIds .= "</td></tr>";
	}
	else
	{
		$assemblyExtraIds .= "None.</td></tr>";
	}
	
	if($target[1] eq 'library')
	{
		$target[2] = "<a onclick='openDialog(\"libraryView.cgi?libraryId=$target[0]\")' title='View'>$target[2]<span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-arrow-1-ne' title='View'></span></a>";
	}
	elsif($target[1] eq 'genome')
	{
		$target[2] = "<a onclick='openDialog(\"genomeView.cgi?genomeId=$target[0]\")' title='$target[8]'>$target[2]<span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-arrow-1-ne' title='View'></span></a>";
	}
	$target[1] = ucfirst ($target[1]);
	$assemblyList .= <<END;
	</ul>$assembledCtg <div id='tabs-assemblyAbout$assemblyId$$'><h3>About '$assembly[2]'</h3><table>
	<tr><td style='text-align:right'><b>Assembly Name</b></td><td>$assembly[2]<br>Version $assembly[3] <sup class='ui-state-disabled'>by $assembly[9] on $assembly[10]</sup></td></tr>
	<tr><td style='text-align:right'><b>Source $target[1]</b></td><td>$target[2]</td></tr>
	$fpcOrAgpId
	$refGenomeId
	$assemblyExtraIds
	<tr><td style='text-align:right'><b>Description</b></td><td>$assemblyDetails->{'description'}</td></tr>
	</table><hr>
	<h3>Run Logs</h3>$assemblyDetails->{'log'}</div> $assemblySortableStyle</div>
END

	undef $/;# enable slurp mode
	my $html = <DATA>;
	$html =~ s/\$assemblyList/$assemblyList/g;
	$html =~ s/\$assemblyId/$assemblyId/g;
	$html =~ s/\$assemblySortableJs/$assemblySortableJs/g;
	$html =~ s/\$assemblyStatsSummaryJs/$assemblyStatsSummaryJs/g;
	$html =~ s/\$autoSeqSearchUrl/$autoSeqSearchUrl/g;
	$html =~ s/\$commoncfg->{HTDOCS}/$commoncfg->{HTDOCS}/g;
	$html =~ s/\$\$/$$/g;

	print header(-cookie=>cookie(-name=>'assembly',-value=>$assemblyId));
	print $html;
}
else
{
	print header(-type=>'text/html',-status=>'402 Invalid operation');
	exit;
}

__DATA__
$assemblyList
<script>
buttonInit();
$("#assemblyTab$assemblyId$$").tabs({
	// loading spinner
	beforeLoad: function(event, ui) {
		ui.panel.html('<img src="$commoncfg->{HTDOCS}/css/images/loading.gif" style="vertical-align:middle;"> Loading...');
	},
	active: 0
});
$( ".assemblyMenu" ).menu();
$( "#assemblySearchSeqName$$" ).autocomplete({
	source: "$autoSeqSearchUrl",
	minLength: 3,
	select: function( event, ui ) {
		openDialog("seqView.cgi?seqId=" + ui.item.id);
	}
});
$assemblySortableJs
$assemblyStatsSummaryJs
loadingHide();

</script>
