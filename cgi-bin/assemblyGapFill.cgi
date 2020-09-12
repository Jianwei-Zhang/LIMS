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
my $assemblyId = param ('assemblyId') || '';
my $ctgOne = param ('ctgOne') || '';
my $ctgTwo = param ('ctgTwo') || '';
$ctgOne =~ s/\D//g;
$ctgTwo =~ s/\D//g;

unless (-e $commoncfg->{TMPDIR})
{
	mkdir $commoncfg->{TMPDIR};
}
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
	0 => "NA",
	1 => "f",
	2 => "r"
	);

print header;

if($assemblyId)
{
	#connect to the mysql server
	my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});
	my $assembly=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$assembly->execute($assemblyId);
	my @assembly = $assembly->fetchrow_array();

	my $target=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$target->execute($assembly[4]);
	my @target = $target->fetchrow_array();
	if($target[1] eq 'genome')
	{
		print <<END;
<script>
parent.errorPop("Sorry...this function is not ready yet!");
parent.loadingHide();
</script>	
END
		exit;	
	}
	elsif($target[1] eq 'library')
	{
		if($ctgOne && $ctgTwo)
		{
			my $assemblyCtgOne=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblyCtg' AND o = ? AND name LIKE ?");
			$assemblyCtgOne->execute($assemblyId,$ctgOne);
			if($assemblyCtgOne->rows < 1)
			{
				print <<END;
<script>
parent.errorPop("Contig-1 number NOT found!");
parent.loadingHide();
</script>	
END
				exit;
			}
			my @assemblyCtgOne = $assemblyCtgOne->fetchrow_array();
			my $assemblyCtgOneSeedSeq;
			foreach (split ",", $assemblyCtgOne[8])
			{
				next if(/^-/);
				$_ =~ s/[^a-zA-Z0-9]//g;
				$assemblyCtgOneSeedSeq = $_;
			}
			my $assemblySeqOne=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
			$assemblySeqOne->execute($assemblyCtgOneSeedSeq);
			my @assemblySeqOne = $assemblySeqOne->fetchrow_array();
			my $assemblySequenceOne=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
			$assemblySequenceOne->execute($assemblySeqOne[5]);
			my @assemblySequenceOne = $assemblySequenceOne->fetchrow_array();
			my $fpcCtgOne;
			my $fpcCtgOneLeft;
			my $fpcCtgOneRight;
			my $getFpcSeqOne = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'fpcClone' AND o = ? AND name LIKE ?");
			$getFpcSeqOne->execute($assembly[6],$assemblySeqOne[2]);
			while (my @getFpcSeqOne = $getFpcSeqOne->fetchrow_array())
			{
				if ($getFpcSeqOne[8] =~ /Map "(.*)" Ends Left/)
				{
					$fpcCtgOne = ucfirst ($1);
				}
				if ($getFpcSeqOne[8] =~ /Ends Left (\d*)/)
				{
					$fpcCtgOneLeft = $1;
				}
				if ($getFpcSeqOne[8] =~ /Ends Right (\d*)/)
				{
					$fpcCtgOneRight = $1;
				}
			}
			my $matchedSeqOne;
			my $matchedSeqOneByTag;
			my $getTagsOne = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'tag' AND x = ? AND name LIKE ?");
			$getTagsOne->execute($assembly[4],$assemblySeqOne[2]);
			my $wgpTagOne = $getTagsOne->rows;
			while(my @getTagsOne = $getTagsOne->fetchrow_array())
			{
				my $getTagsOneMatched = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'tag' AND x = ? AND o = ? AND name NOT LIKE ?");
				$getTagsOneMatched->execute($assembly[4],$getTagsOne[3],$assemblySeqOne[2]);
				while(my @getTagsOneMatched = $getTagsOneMatched->fetchrow_array())
				{
					$matchedSeqOne->{$getTagsOneMatched[2]} = 1;
					$matchedSeqOneByTag->{$getTagsOneMatched[2]} = 0 if (!exists $matchedSeqOneByTag->{$getTagsOneMatched[2]});
					$matchedSeqOneByTag->{$getTagsOneMatched[2]}++;
				}
			}

			my $queryIdOnSeqOne;
			my $queryDirOne;
			for (my $position = 0; $position < length($assemblySeqOne[5]); $position += 2)
			{
				$queryDirOne .= "/q". substr($assemblySeqOne[5],$position,2);
			}
			if(-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirOne")
			{
				my @dirs;
				push @dirs, "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirOne";

				do{
					my $dir = shift @dirs;
					opendir(DIR, $dir) or die "can't opendir $dir: $!";
					my @files = readdir(DIR);
					closedir DIR;
					foreach my $file (sort @files)
					{
						next if ($file =~ /^\./);
						if(-f "$dir/$file")
						{
							if ($file =~ /(\d+)-(\d+).tbl$/)
							{
								my $subjectId = $2;
								open (TBL, "$dir/$file") or die "can't open file: $dir/$file";
								while(<TBL>)
								{
									chop;
									/^#/ and next;
									my @hit = split("\t",$_);
									next if ($hit[12] > 0);
									$queryIdOnSeqOne->{$subjectId} = $hit[3] if (!exists $queryIdOnSeqOne->{$subjectId});
									$queryIdOnSeqOne->{$subjectId} = $hit[3] if ($queryIdOnSeqOne->{$subjectId} < $hit[3]);
								}
								close(TBL);
							}
						}
						elsif(-d "$dir/$file")
						{
							unshift @dirs, "$dir/$file";
						}
					}
				}while (@dirs);
			}

			my $matchedSeqOneByBes;
			my $getBesOne = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 98 AND x = ?");
			$getBesOne->execute($assembly[4]);
			while(my @getBesOne = $getBesOne->fetchrow_array())
			{
				if(exists $queryIdOnSeqOne->{$getBesOne[0]})
				{
					$matchedSeqOne->{$getBesOne[2]} = 1;
					$matchedSeqOneByBes->{$getBesOne[2]} = "$seqDir{$getBesOne[6]}-end";
				}
			}

			my $matchedSeqOneByGenome;
			my $getGenomeOne = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 99 AND x = ?");
			$getGenomeOne->execute($assembly[4]);
			while(my @getGenomeOne = $getGenomeOne->fetchrow_array())
			{
				if(exists $queryIdOnSeqOne->{$getGenomeOne[0]})
				{
					$matchedSeqOne->{$getGenomeOne[2]} = 1;
					$matchedSeqOneByGenome->{$getGenomeOne[2]} = 0 if (!exists $matchedSeqOneByGenome->{$getGenomeOne[2]});
					$matchedSeqOneByGenome->{$getGenomeOne[2]} = $queryIdOnSeqOne->{$getGenomeOne[0]} if($matchedSeqOneByGenome->{$getGenomeOne[2]} < $queryIdOnSeqOne->{$getGenomeOne[0]});
				}
			}

			my $assemblyCtgTwo=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblyCtg' AND o = ? AND name LIKE ?");
			$assemblyCtgTwo->execute($assemblyId,$ctgTwo);
			if($assemblyCtgTwo->rows < 1)
			{
				print <<END;
<script>
parent.errorPop("Contig-2 number NOT found!");
parent.loadingHide();
</script>	
END
				exit;
			}
			my @assemblyCtgTwo = $assemblyCtgTwo->fetchrow_array();
			my $assemblyCtgTwoSeedSeq;
			foreach (split ",", $assemblyCtgTwo[8])
			{
				next if(/^-/);
				$_ =~ s/[^a-zA-Z0-9]//g;
				$assemblyCtgTwoSeedSeq = $_;
				last;
			}
			my $assemblySeqTwo=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
			$assemblySeqTwo->execute($assemblyCtgTwoSeedSeq);
			my @assemblySeqTwo = $assemblySeqTwo->fetchrow_array();
			my $assemblySequenceTwo=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
			$assemblySequenceTwo->execute($assemblySeqTwo[5]);
			my @assemblySequenceTwo = $assemblySequenceTwo->fetchrow_array();
			my $fpcCtgTwo = 'NA';
			my $fpcCtgTwoLeft = '0';
			my $fpcCtgTwoRight = '0';
			my $getFpcSeqTwo = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'fpcClone' AND o = ? AND name LIKE ?");
			$getFpcSeqTwo->execute($assembly[6],$assemblySeqTwo[2]);
			while (my @getFpcSeqTwo = $getFpcSeqTwo->fetchrow_array())
			{
				if ($getFpcSeqTwo[8] =~ /Map "(.*)" Ends Left/)
				{
					$fpcCtgTwo = ucfirst ($1);
				}
				if ($getFpcSeqTwo[8] =~ /Ends Left (\d*)/)
				{
					$fpcCtgTwoLeft = $1;
				}
				if ($getFpcSeqTwo[8] =~ /Ends Right (\d*)/)
				{
					$fpcCtgTwoRight = $1;
				}
			}
			my $matchedSeqTwoByTag;
			my $getTagsTwo = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'tag' AND x = ? AND name LIKE ?");
			$getTagsTwo->execute($assembly[4],$assemblySeqTwo[2]);
			my $wgpTagTwo = $getTagsTwo->rows;
			while(my @getTagsTwo = $getTagsTwo->fetchrow_array())
			{
				my $getTagsTwoMatched = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'tag' AND x = ? AND o = ? AND name NOT LIKE ?");
				$getTagsTwoMatched->execute($assembly[4],$getTagsTwo[3],$assemblySeqTwo[2]);
				while(my @getTagsTwoMatched = $getTagsTwoMatched->fetchrow_array())
				{
					$matchedSeqTwoByTag->{$getTagsTwoMatched[2]} = 0 if (!exists $matchedSeqTwoByTag->{$getTagsTwoMatched[2]});
					$matchedSeqTwoByTag->{$getTagsTwoMatched[2]}++;
				}
			}

			my $queryIdOnSeqTwo;
			my $queryDirTwo;
			for (my $position = 0; $position < length($assemblySeqTwo[5]); $position += 2)
			{
				$queryDirTwo .= "/q". substr($assemblySeqTwo[5],$position,2);
			}
			if(-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirTwo")
			{
				my @dirs;
				push @dirs, "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirTwo";

				do{
					my $dir = shift @dirs;
					opendir(DIR, $dir) or die "can't opendir $dir: $!";
					my @files = readdir(DIR);
					closedir DIR;
					foreach my $file (sort @files)
					{
						next if ($file =~ /^\./);
						if(-f "$dir/$file")
						{
							if ($file =~ /(\d+)-(\d+).tbl$/)
							{
								my $subjectId = $2;
								open (TBL, "$dir/$file") or die "can't open file: $dir/$file";
								while(<TBL>)
								{
									chop;
									/^#/ and next;
									my @hit = split("\t",$_);
									next if ($hit[12] > 0);
									$queryIdOnSeqTwo->{$subjectId} = $hit[3] if (!exists $queryIdOnSeqTwo->{$subjectId});
									$queryIdOnSeqTwo->{$subjectId} = $hit[3] if ($queryIdOnSeqTwo->{$subjectId} < $hit[3]);
								}
								close(TBL);
							}
						}
						elsif(-d "$dir/$file")
						{
							unshift @dirs, "$dir/$file";
						}
					}
				}while (@dirs);
			}

			my $matchedSeqTwoByBes;
			my $getBesTwo = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 98 AND x = ?");
			$getBesTwo->execute($assembly[4]);
			while(my @getBesTwo = $getBesTwo->fetchrow_array())
			{
				if(exists $queryIdOnSeqTwo->{$getBesTwo[0]})
				{
					$matchedSeqTwoByBes->{$getBesTwo[2]} = "$seqDir{$getBesTwo[6]}-end";
				}
			}

			my $matchedSeqTwoByGenome;
			my $getGenomeTwo = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 99 AND x = ?");
			$getGenomeTwo->execute($assembly[4]);
			while(my @getGenomeTwo = $getGenomeTwo->fetchrow_array())
			{
				if(exists $queryIdOnSeqTwo->{$getGenomeTwo[0]})
				{
					$matchedSeqTwoByGenome->{$getGenomeTwo[2]} = 0 if (!exists $matchedSeqTwoByGenome->{$getGenomeTwo[2]});
					$matchedSeqTwoByGenome->{$getGenomeTwo[2]} = $queryIdOnSeqTwo->{$getGenomeTwo[0]} if ($matchedSeqTwoByGenome->{$getGenomeTwo[2]} < $queryIdOnSeqTwo->{$getGenomeTwo[0]});
				}
			}


			open (GAPFILL,">$commoncfg->{TMPDIR}/$assemblyId-$ctgOne-$ctgTwo.gapfill.html") or die "can't open file: $commoncfg->{TMPDIR}/$assemblyId-$ctgOne-$ctgTwo.gapfill.html";
			print GAPFILL <<END;
<div id='gapFillerDiv$$' name='gapFillerDiv$$'>
<table id='gapFiller$$' class='display'>
<thead>
<tr><th colspan='2'>Filler Candidates</th>
<th colspan='4'><b>Ctg$ctgOne</b> R-end: <a title='$seqType{$assemblySequenceOne[3]}'>$assemblySeqOne[2]($assemblyCtgOneSeedSeq)</a> ($bacAssignType{$assemblySequenceOne[7]}: <a title='mapped tags'>$assemblySequenceOne[6]</a>)</th>
<tr><th>Clone</th><th>Score</th><th>Tag: $wgpTagOne</th><th>BES</th><th>Genome</th><th>FPC $fpcCtgOne ($fpcCtgOneLeft-$fpcCtgOneRight)</th></tr>
</thead>
<tfoot>
<tr><th>Clone</th><th>Score</th><th>Tag: $wgpTagTwo</th><th>BES</th><th>Genome</th><th>FPC $fpcCtgTwo ($fpcCtgTwoLeft-$fpcCtgTwoRight)</th></tr>
<tr><th colspan='2'>Filler Candidates</th>
<th colspan='4'><b>Ctg$ctgTwo</b> L-end: <a title='$seqType{$assemblySequenceTwo[3]}'>$assemblySeqTwo[2]($assemblyCtgTwoSeedSeq)</a> ($bacAssignType{$assemblySequenceTwo[7]}: <a title='mapped tags'>$assemblySequenceTwo[6]</a>)</th></tr>
</tfoot>
<tbody>
END
			for (keys %$matchedSeqOne)
			{
				next unless (exists $matchedSeqTwoByTag->{$_} || exists $matchedSeqTwoByBes->{$_} ||  exists $matchedSeqTwoByGenome->{$_});
				my $score = 0;
				if(exists $matchedSeqOneByTag->{$_} && exists $matchedSeqTwoByTag->{$_})
				{
					$score++;
				}
				if(exists $matchedSeqOneByBes->{$_} && exists $matchedSeqTwoByBes->{$_} && ($matchedSeqOneByBes->{$_} ne $matchedSeqTwoByBes->{$_}))
				{
					$score++;
				}
				if(exists $matchedSeqTwoByGenome->{$_} && exists $matchedSeqTwoByGenome->{$_})
				{
					$score++;
				}
				my $fpcCtgMatch = 'NA';
				my $fpcCtgMatchLeft = '0';
				my $fpcCtgMatchRight = '0';
				my $fpcCtgMTPClone = 0;
				my $fpcCtgSequencedClone = '';
				my $getFpcCloneMatch = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'fpcClone' AND o = ? AND name LIKE ?");
				$getFpcCloneMatch->execute($assembly[6],$_);
				while (my @getFpcCloneMatch = $getFpcCloneMatch->fetchrow_array())
				{
					if ($getFpcCloneMatch[8] =~ /Map "(.*)" Ends Left/)
					{
						$fpcCtgMatch = ucfirst ($1);
					}
					if ($getFpcCloneMatch[8] =~ /Ends Left (\d*)/)
					{
						$fpcCtgMatchLeft = $1;
					}
					if ($getFpcCloneMatch[8] =~ /Ends Right (\d*)/)
					{
						$fpcCtgMatchRight = $1;
					}
					$fpcCtgMTPClone = $getFpcCloneMatch[5];
					$fpcCtgSequencedClone = ($getFpcCloneMatch[4] > 0) ? "<span class='ui-icon ui-icon-check' style='float: left; margin-right: 0;' title='Sequence available'></span>" : '';
				}
				if($fpcCtgMTPClone > 0)
				{
					print GAPFILL "<tr><td><span class='ui-icon ui-icon-flag' style='float: left; margin-right: 0;' title='MTP Clone'></span>$fpcCtgSequencedClone $_ </td><td>$score</td>";
				}
				else
				{
					print GAPFILL "<tr><td>$fpcCtgSequencedClone $_</td><td>$score</td>";
				}

				if(exists $matchedSeqOneByTag->{$_})
				{
					print GAPFILL "<td>$matchedSeqOneByTag->{$_}";
				}
				else
				{
					print GAPFILL "<td>-";
				}
				if(exists $matchedSeqTwoByTag->{$_})
				{
					print GAPFILL "/$matchedSeqTwoByTag->{$_}</td>";
				}
				else
				{
					print GAPFILL "/-</td>";
				}
				if(exists $matchedSeqOneByBes->{$_})
				{
					print GAPFILL "<td>$matchedSeqOneByBes->{$_}";
				}
				else
				{
					print GAPFILL "<td>-";
				}
				if(exists $matchedSeqTwoByBes->{$_})
				{
					print GAPFILL "/$matchedSeqTwoByBes->{$_}</td>";
				}
				else
				{
					print GAPFILL "/-</td>";
				}
				if(exists $matchedSeqOneByGenome->{$_})
				{
					print GAPFILL "<td>$matchedSeqOneByGenome->{$_}";
				}
				else
				{
					print GAPFILL "<td>-";
				}
				if(exists $matchedSeqTwoByGenome->{$_})
				{
					print GAPFILL "/$matchedSeqTwoByGenome->{$_}</td>";
				}
				else
				{
					print GAPFILL "/-</td>";
				}
				print GAPFILL "<td>$fpcCtgMatch ($fpcCtgMatchLeft-$fpcCtgMatchRight)</td></tr>";
			}

			print GAPFILL <<END;
</tbody></table>
</div>
<script>
\$('#dialog').dialog("option", "title", "Gap Filler for Ctg$ctgOne & Ctg$ctgTwo");
\$( "#dialog" ).dialog( "option", "buttons", [ { text: "Fill Another Gap", click: function() {  closeDialog();openDialog('assemblyGapFillForm.cgi?assemblyId=$assemblyId'); } }, { text: "Print", click: function() {printDiv('gapFillerDiv$$'); } },{ text: "OK", click: function() {closeDialog(); } } ] );
\$('#dialog').dialog("option", "width", 1000);
\$( "#gapFiller$$" ).dataTable({
	"scrollY": "500px",
	"scrollCollapse": true,
	"paging": false,
	"searching": false
});
</script>
END
			close (GAPFILL);
			print <<END;
		<script>
		parent.loadingHide();
		parent.closeDialog();
		parent.openDialog('$commoncfg->{TMPURL}/$assemblyId-$ctgOne-$ctgTwo.gapfill.html');
		</script>
END
		}
		else
		{
			open (GAPFILLTEXT,">$commoncfg->{TMPDIR}/$assemblyId.gapfill.txt") or die "can't open file: $commoncfg->{TMPDIR}/$assemblyId.gapfill.txt";
			open (GAPFILL,">$commoncfg->{TMPDIR}/$assemblyId.gapfill.html") or die "can't open file: $commoncfg->{TMPDIR}/$assemblyId.gapfill.html";
			print GAPFILLTEXT "#Gap Filler for Assembly $assembly[2]\n#Contig#\tClone\tTag\tBES\tFPC\tScore\n";
			print GAPFILL <<END;
<a style='float: right;' href='$commoncfg->{TMPURL}/$assemblyId.gapfill.txt.gz' target='hiddenFrame'><span style='left: 0px;top: 0px;display:inline-block;' class='ui-icon ui-icon-disk' title='Save'></span>Download Table</a>
<div id='gapFillerDiv$$' name='gapFillerDiv$$'>
<table id='gapFiller$$' class='display'>
<thead>
<tr><th><b>Contig#</b></th><th>Clone</th><th>Tag</th><th>BES</th><th>FPC</th><th>Score</th></tr>
</thead>
<tbody>
END

			my $lastAssemblyCtgChr = '';
			my $lastAssemblyCtgNumber = '';
			my $lastAssemblyCtgSeqs = '';
			my $chrStart = 0;
			my $assemblyCtgs = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblyCtg' AND o = ? ORDER BY x,z");
			$assemblyCtgs->execute($assemblyId);
			while(my @assemblyCtgs = $assemblyCtgs->fetchrow_array())
			{
				if($assemblyCtgs[4] ne $lastAssemblyCtgChr)
				{
					$lastAssemblyCtgChr = $assemblyCtgs[4];
					$lastAssemblyCtgNumber = $assemblyCtgs[2];
					$lastAssemblyCtgSeqs = $assemblyCtgs[8];
					$chrStart = 1;
					next;
				}
				#Contig A
				my $assemblyCtgOneSeedSeq;
				foreach (split ",", $lastAssemblyCtgSeqs)
				{
					next if(/^-/);
					$_ =~ s/[^a-zA-Z0-9]//g;
					$assemblyCtgOneSeedSeq = $_;
				}
				my $assemblySeqOne=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
				$assemblySeqOne->execute($assemblyCtgOneSeedSeq);
				my @assemblySeqOne = $assemblySeqOne->fetchrow_array();
				my $assemblySequenceOne=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
				$assemblySequenceOne->execute($assemblySeqOne[5]);
				my @assemblySequenceOne = $assemblySequenceOne->fetchrow_array();
				my $fpcCtgOne = 'NA';
				my $fpcCtgOneLeft = 0;
				my $fpcCtgOneRight = 0;
				my $getFpcSeqOne = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'fpcClone' AND o = ? AND name LIKE ?");
				$getFpcSeqOne->execute($assembly[6],$assemblySeqOne[2]);
				while (my @getFpcSeqOne = $getFpcSeqOne->fetchrow_array())
				{
					if ($getFpcSeqOne[8] =~ /Map "(.*)" Ends Left/)
					{
						$fpcCtgOne = ucfirst ($1);
					}
					if ($getFpcSeqOne[8] =~ /Ends Left (\d*)/)
					{
						$fpcCtgOneLeft = $1;
					}
					if ($getFpcSeqOne[8] =~ /Ends Right (\d*)/)
					{
						$fpcCtgOneRight = $1;
					}
				}
				my $matchedSeqOne;
				my $matchedSeqOneByTag;
				my $getTagsOne = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'tag' AND x = ? AND name LIKE ?");
				$getTagsOne->execute($assembly[4],$assemblySeqOne[2]);
				my $wgpTagOne = $getTagsOne->rows;
				while(my @getTagsOne = $getTagsOne->fetchrow_array())
				{
					my $getTagsOneMatched = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'tag' AND x = ? AND o = ? AND name NOT LIKE ?");
					$getTagsOneMatched->execute($assembly[4],$getTagsOne[3],$assemblySeqOne[2]);
					while(my @getTagsOneMatched = $getTagsOneMatched->fetchrow_array())
					{
						$matchedSeqOne->{$getTagsOneMatched[2]} = 1;
						$matchedSeqOneByTag->{$getTagsOneMatched[2]} = 0 if (!exists $matchedSeqOneByTag->{$getTagsOneMatched[2]});
						$matchedSeqOneByTag->{$getTagsOneMatched[2]}++;
					}
				}

				my $queryIdOnSeqOne;
				my $queryDirOne;
				for (my $position = 0; $position < length($assemblySeqOne[5]); $position += 2)
				{
					$queryDirOne .= "/q". substr($assemblySeqOne[5],$position,2);
				}
				if(-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirOne")
				{
					my @dirs;
					push @dirs, "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirOne";

					do{
						my $dir = shift @dirs;
						opendir(DIR, $dir) or die "can't opendir $dir: $!";
						my @files = readdir(DIR);
						closedir DIR;
						foreach my $file (sort @files)
						{
							next if ($file =~ /^\./);
							if(-f "$dir/$file")
							{
								if ($file =~ /(\d+)-(\d+).tbl$/)
								{
									my $subjectId = $2;
									open (TBL, "$dir/$file") or die "can't open file: $dir/$file";
									while(<TBL>)
									{
										chop;
										/^#/ and next;
										my @hit = split("\t",$_);
										next if ($hit[12] > 0);
										$queryIdOnSeqOne->{$subjectId} = $hit[3] if (!exists $queryIdOnSeqOne->{$subjectId});
										$queryIdOnSeqOne->{$subjectId} = $hit[3] if ($queryIdOnSeqOne->{$subjectId} < $hit[3]);
									}
									close(TBL);
								}
							}
							elsif(-d "$dir/$file")
							{
								unshift @dirs, "$dir/$file";
							}
						}
					}while (@dirs);
				}

				my $matchedSeqOneByBes;
				my $getBesOne = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 98 AND x = ?");
				$getBesOne->execute($assembly[4]);
				while(my @getBesOne = $getBesOne->fetchrow_array())
				{
					if(exists $queryIdOnSeqOne->{$getBesOne[0]})
					{
						$matchedSeqOne->{$getBesOne[2]} = 1;
						$matchedSeqOneByBes->{$getBesOne[2]} = "$seqDir{$getBesOne[6]}-end";
					}
				}

				my $matchedSeqOneByGenome;
				my $getGenomeOne = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 99 AND x = ?");
				$getGenomeOne->execute($assembly[4]);
				while(my @getGenomeOne = $getGenomeOne->fetchrow_array())
				{
					if(exists $queryIdOnSeqOne->{$getGenomeOne[0]})
					{
						$matchedSeqOne->{$getGenomeOne[2]} = 1;
						$matchedSeqOneByGenome->{$getGenomeOne[2]} = 0 if (!exists $matchedSeqOneByGenome->{$getGenomeOne[2]});
						$matchedSeqOneByGenome->{$getGenomeOne[2]} = $queryIdOnSeqOne->{$getGenomeOne[0]} if($matchedSeqOneByGenome->{$getGenomeOne[2]} < $queryIdOnSeqOne->{$getGenomeOne[0]});
					}
				}

				#Contig B
				my $assemblyCtgTwoSeedSeq;
				foreach (split ",", $assemblyCtgs[8])
				{
					next if(/^-/);
					$_ =~ s/[^a-zA-Z0-9]//g;
					$assemblyCtgTwoSeedSeq = $_;
					last;
				}
				my $assemblySeqTwo=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
				$assemblySeqTwo->execute($assemblyCtgTwoSeedSeq);
				my @assemblySeqTwo = $assemblySeqTwo->fetchrow_array();
				my $assemblySequenceTwo=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
				$assemblySequenceTwo->execute($assemblySeqTwo[5]);
				my @assemblySequenceTwo = $assemblySequenceTwo->fetchrow_array();
				my $fpcCtgTwo = 'NA';
				my $fpcCtgTwoLeft = '0';
				my $fpcCtgTwoRight = '0';
				my $getFpcSeqTwo = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'fpcClone' AND o = ? AND name LIKE ?");
				$getFpcSeqTwo->execute($assembly[6],$assemblySeqTwo[2]);
				while (my @getFpcSeqTwo = $getFpcSeqTwo->fetchrow_array())
				{
					if ($getFpcSeqTwo[8] =~ /Map "(.*)" Ends Left/)
					{
						$fpcCtgTwo = ucfirst ($1);
					}
					if ($getFpcSeqTwo[8] =~ /Ends Left (\d*)/)
					{
						$fpcCtgTwoLeft = $1;
					}
					if ($getFpcSeqTwo[8] =~ /Ends Right (\d*)/)
					{
						$fpcCtgTwoRight = $1;
					}
				}
				my $matchedSeqTwoByTag;
				my $getTagsTwo = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'tag' AND x = ? AND name LIKE ?");
				$getTagsTwo->execute($assembly[4],$assemblySeqTwo[2]);
				my $wgpTagTwo = $getTagsTwo->rows;
				while(my @getTagsTwo = $getTagsTwo->fetchrow_array())
				{
					my $getTagsTwoMatched = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'tag' AND x = ? AND o = ? AND name NOT LIKE ?");
					$getTagsTwoMatched->execute($assembly[4],$getTagsTwo[3],$assemblySeqTwo[2]);
					while(my @getTagsTwoMatched = $getTagsTwoMatched->fetchrow_array())
					{
						$matchedSeqTwoByTag->{$getTagsTwoMatched[2]} = 0 if (!exists $matchedSeqTwoByTag->{$getTagsTwoMatched[2]});
						$matchedSeqTwoByTag->{$getTagsTwoMatched[2]}++;
					}
				}

				my $queryIdOnSeqTwo;
				my $queryDirTwo;
				for (my $position = 0; $position < length($assemblySeqTwo[5]); $position += 2)
				{
					$queryDirTwo .= "/q". substr($assemblySeqTwo[5],$position,2);
				}
				if(-e "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirTwo")
				{
					my @dirs;
					push @dirs, "$commoncfg->{DATADIR}/alignments/seqToSeq$queryDirTwo";

					do{
						my $dir = shift @dirs;
						opendir(DIR, $dir) or die "can't opendir $dir: $!";
						my @files = readdir(DIR);
						closedir DIR;
						foreach my $file (sort @files)
						{
							next if ($file =~ /^\./);
							if(-f "$dir/$file")
							{
								if ($file =~ /(\d+)-(\d+).tbl$/)
								{
									my $subjectId = $2;
									open (TBL, "$dir/$file") or die "can't open file: $dir/$file";
									while(<TBL>)
									{
										chop;
										/^#/ and next;
										my @hit = split("\t",$_);
										next if ($hit[12] > 0);
										$queryIdOnSeqTwo->{$subjectId} = $hit[3] if (!exists $queryIdOnSeqTwo->{$subjectId});
										$queryIdOnSeqTwo->{$subjectId} = $hit[3] if ($queryIdOnSeqTwo->{$subjectId} < $hit[3]);
									}
									close(TBL);
								}
							}
							elsif(-d "$dir/$file")
							{
								unshift @dirs, "$dir/$file";
							}
						}
					}while (@dirs);
				}

				my $matchedSeqTwoByBes;
				my $getBesTwo = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 98 AND x = ?");
				$getBesTwo->execute($assembly[4]);
				while(my @getBesTwo = $getBesTwo->fetchrow_array())
				{
					if(exists $queryIdOnSeqTwo->{$getBesTwo[0]})
					{
						$matchedSeqTwoByBes->{$getBesTwo[2]} = "$seqDir{$getBesTwo[6]}-end";
					}
				}

				my $matchedSeqTwoByGenome;
				my $getGenomeTwo = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 99 AND x = ?");
				$getGenomeTwo->execute($assembly[4]);
				while(my @getGenomeTwo = $getGenomeTwo->fetchrow_array())
				{
					if(exists $queryIdOnSeqTwo->{$getGenomeTwo[0]})
					{
						$matchedSeqTwoByGenome->{$getGenomeTwo[2]} = 0 if (!exists $matchedSeqTwoByGenome->{$getGenomeTwo[2]});
						$matchedSeqTwoByGenome->{$getGenomeTwo[2]} = $queryIdOnSeqTwo->{$getGenomeTwo[0]} if ($matchedSeqTwoByGenome->{$getGenomeTwo[2]} < $queryIdOnSeqTwo->{$getGenomeTwo[0]});
					}
				}

				if($chrStart)
				{
					print GAPFILLTEXT ">Chromosome $assemblyCtgs[4]\t\t\t\t\t\t\n";
					print GAPFILLTEXT "Ctg$lastAssemblyCtgNumber\tR-end: $assemblySeqOne[2]($assemblyCtgOneSeedSeq) ($seqType{$assemblySequenceOne[3]},$bacAssignType{$assemblySequenceOne[7]}: $assemblySequenceOne[6])\t$wgpTagOne\tBES\tGenome\t$fpcCtgOne ($fpcCtgOneLeft-$fpcCtgOneRight)\t\n";
					print GAPFILL <<END;
<tr><td class="ui-state-highlight" style='white-space: nowrap;'><span class='ui-icon ui-icon ui-icon-carat-1-e' style='float: left; margin-right: 0;'></span><b>Chromosome $assemblyCtgs[4]</b></td><td><hr class="ui-state-highlight ui-corner-all" /></td><td></td><td></td><td></td><td></td><td></td></tr>
<tr><td style='text-align:right' title='Chromosome $assemblyCtgs[4]'><b>Ctg$lastAssemblyCtgNumber</b><span class='ui-icon ui-icon-carat-1-sw' style='float: right; margin-right: 0;' title='R-end'></span></td><td>R-end: <a title='$seqType{$assemblySequenceOne[3]}'>$assemblySeqOne[2]($assemblyCtgOneSeedSeq)</a> ($bacAssignType{$assemblySequenceOne[7]}: <a title='mapped tags'>$assemblySequenceOne[6]</a>)</td><td>$wgpTagOne</td><td>BES</td><td>Genome</td><td>$fpcCtgOne ($fpcCtgOneLeft-$fpcCtgOneRight)</td><td></td></tr>
END
				}
				else
				{
					print GAPFILLTEXT "\tR-end: $assemblySeqOne[2]($assemblyCtgOneSeedSeq) ($seqType{$assemblySequenceOne[3]},$bacAssignType{$assemblySequenceOne[7]}: $assemblySequenceOne[6])\t$wgpTagOne\tBES\tGenome\t$fpcCtgOne ($fpcCtgOneLeft-$fpcCtgOneRight)\t\n";
					print GAPFILL <<END;
<tr><td style='text-align:right' title='Chromosome $assemblyCtgs[4]'><span class='ui-icon ui-icon-carat-1-sw' style='float: right; margin-right: 0;' title='R-end'></span></td><td>R-end: <a title='$seqType{$assemblySequenceOne[3]}'>$assemblySeqOne[2]($assemblyCtgOneSeedSeq)</a> ($bacAssignType{$assemblySequenceOne[7]}: <a title='mapped tags'>$assemblySequenceOne[6]</a>)</td><td>$wgpTagOne</td><td>BES</td><td>Genome</td><td>$fpcCtgOne ($fpcCtgOneLeft-$fpcCtgOneRight)</td><td></td></tr>
END
				}
				my $candidates = 0;
				for (keys %$matchedSeqOne)
				{
					next unless (exists $matchedSeqTwoByTag->{$_} || exists $matchedSeqTwoByBes->{$_}  || exists $matchedSeqTwoByGenome->{$_});
					$candidates++;
					my $score = 0;
					if(exists $matchedSeqOneByTag->{$_} && exists $matchedSeqTwoByTag->{$_})
					{
						$score++;
					}
					if(exists $matchedSeqOneByBes->{$_} && exists $matchedSeqTwoByBes->{$_} && ($matchedSeqOneByBes->{$_} ne $matchedSeqTwoByBes->{$_}))
					{
						$score++;
					}
					if(exists $matchedSeqTwoByGenome->{$_} && exists $matchedSeqTwoByGenome->{$_})
					{
						$score++;
					}
					my $fpcCtgMatch = 'NA';
					my $fpcCtgMatchLeft = '0';
					my $fpcCtgMatchRight = '0';
					my $fpcCtgMTPClone = 0;
					my $fpcCtgSequencedClone = '';
					my $fpcCtgSequencedCloneMark = '';
					my $getFpcCloneMatch = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'fpcClone' AND o = ? AND name LIKE ?");
					$getFpcCloneMatch->execute($assembly[6],$_);
					while (my @getFpcCloneMatch = $getFpcCloneMatch->fetchrow_array())
					{
						if ($getFpcCloneMatch[8] =~ /Map "(.*)" Ends Left/)
						{
							$fpcCtgMatch = ucfirst ($1);
						}
						if ($getFpcCloneMatch[8] =~ /Ends Left (\d*)/)
						{
							$fpcCtgMatchLeft = $1;
						}
						if ($getFpcCloneMatch[8] =~ /Ends Right (\d*)/)
						{
							$fpcCtgMatchRight = $1;
						}
						$fpcCtgMTPClone = $getFpcCloneMatch[5];
						$fpcCtgSequencedClone = ($getFpcCloneMatch[4] > 0) ? "<span class='ui-icon ui-icon-check' style='float: left; margin-right: 0;' title='Sequence available'></span>" : '';
						$fpcCtgSequencedCloneMark = ($getFpcCloneMatch[4] > 0) ? "+" : '';
					}
					if($fpcCtgMTPClone > 0)
					{
						print GAPFILLTEXT "*|\t$fpcCtgSequencedCloneMark$_";
						print GAPFILL "<tr><td style='text-align:right' title='Chromosome $assemblyCtgs[4]'><span class='ui-icon ui-icon-grip-solid-vertical' style='float: right; margin-right: 0;' title='Filler'></span><span class='ui-icon ui-icon-flag' style='float: right; margin-right: 0;' title='MTP Clone'></span></td><td>$fpcCtgSequencedClone $_</td>";
					}
					else
					{
						print GAPFILLTEXT "|\t$fpcCtgSequencedCloneMark$_";
						print GAPFILL "<tr><td style='text-align:right' title='Chromosome $assemblyCtgs[4]'><span class='ui-icon ui-icon-grip-solid-vertical' style='float: right; margin-right: 0;' title='Filler'></span></td><td>$fpcCtgSequencedClone $_</td>";
					}

					if(exists $matchedSeqOneByTag->{$_})
					{
						print GAPFILLTEXT "\t$matchedSeqOneByTag->{$_}";
						print GAPFILL "<td>$matchedSeqOneByTag->{$_}";
					}
					else
					{
						print GAPFILLTEXT "\t-";
						print GAPFILL "<td>-";
					}
					if(exists $matchedSeqTwoByTag->{$_})
					{
						print GAPFILLTEXT "/$matchedSeqTwoByTag->{$_}";
						print GAPFILL "/$matchedSeqTwoByTag->{$_}</td>";
					}
					else
					{
						print GAPFILLTEXT "/-";
						print GAPFILL "/-</td>";
					}
					if(exists $matchedSeqOneByBes->{$_})
					{
						print GAPFILLTEXT "\t$matchedSeqOneByBes->{$_}";
						print GAPFILL "<td>$matchedSeqOneByBes->{$_}";
					}
					else
					{
						print GAPFILLTEXT "\t-";
						print GAPFILL "<td>-";
					}
					if(exists $matchedSeqTwoByBes->{$_})
					{
						print GAPFILLTEXT "/$matchedSeqTwoByBes->{$_}";
						print GAPFILL "/$matchedSeqTwoByBes->{$_}</td>";
					}
					else
					{
						print GAPFILLTEXT "/-";
						print GAPFILL "/-</td>";
					}
					if(exists $matchedSeqOneByGenome->{$_})
					{
						print GAPFILLTEXT "\t$matchedSeqOneByGenome->{$_}";
						print GAPFILL "<td>$matchedSeqOneByGenome->{$_}";
					}
					else
					{
						print GAPFILLTEXT "\t-";
						print GAPFILL "<td>-";
					}
					if(exists $matchedSeqTwoByGenome->{$_})
					{
						print GAPFILLTEXT "\t$matchedSeqTwoByGenome->{$_}";
						print GAPFILL "/$matchedSeqTwoByGenome->{$_}</td>";
					}
					else
					{
						print GAPFILLTEXT "/-";
						print GAPFILL "/-</td>";
					}

					print GAPFILLTEXT "\t$fpcCtgMatch ($fpcCtgMatchLeft-$fpcCtgMatchRight)\t$score\n";
					print GAPFILL "<td>$fpcCtgMatch ($fpcCtgMatchLeft-$fpcCtgMatchRight)</td><td>$score</td></tr>";
				}
				unless ($candidates)
				{
					print GAPFILLTEXT "\tNo Filler Found\n";
					print GAPFILL <<END;
<tr><td></td><td><div class="ui-state-error ui-corner-all">No Filler Found</div></td><td></td><td></td><td></td><td></td></tr>
END
				}
				print GAPFILLTEXT "Ctg$assemblyCtgs[2]\tL-end: $assemblySeqTwo[2]($assemblyCtgTwoSeedSeq) ($seqType{$assemblySequenceTwo[3]},$bacAssignType{$assemblySequenceTwo[7]}: $assemblySequenceTwo[6])\t$wgpTagTwo\tBES\tGenome\t$fpcCtgTwo ($fpcCtgTwoLeft-$fpcCtgTwoRight)\t\n";
				print GAPFILL <<END;
<tr><td style='text-align:right' title='Chromosome $assemblyCtgs[4]'><b>Ctg$assemblyCtgs[2]</b><span class='ui-icon ui-icon-carat-1-nw' style='float: right; margin-right: 0;' title='L-end'></span></td><td>L-end: <a title='$seqType{$assemblySequenceTwo[3]}'>$assemblySeqTwo[2]($assemblyCtgTwoSeedSeq)</a> ($bacAssignType{$assemblySequenceTwo[7]}: <a title='mapped tags'>$assemblySequenceTwo[6]</a>)</td><td>$wgpTagTwo</td><td>BES</td><td>Genome</td><td>$fpcCtgTwo ($fpcCtgTwoLeft-$fpcCtgTwoRight)</td><td></td></tr>
END

				$lastAssemblyCtgChr = $assemblyCtgs[4];
				$lastAssemblyCtgNumber = $assemblyCtgs[2];
				$lastAssemblyCtgSeqs = $assemblyCtgs[8];
				$chrStart = 0;
			}

			print GAPFILL <<END;
</tbody></table>
</div>
<script>
\$('#dialog').dialog("option", "title", "Gap Filler for Assembly $assembly[2]");
\$( "#dialog" ).dialog( "option", "buttons", [ { text: "reRun Gap Filler", click: function() {  closeDialog();openDialog('assemblyGapFillForm.cgi?assemblyId=$assemblyId'); } }, { text: "Print", click: function() {printDiv('gapFillerDiv$$'); } },{ text: "OK", click: function() {closeDialog(); } } ] );
\$('#dialog').dialog("option", "width", 1000);
\$( "#gapFiller$$" ).dataTable({
	"scrollY": "500px",
	"scrollCollapse": true,
	"paging": false,
	"ordering": false,
	"searching": false
});
</script>
END
			close (GAPFILLTEXT);
			`gzip -f $commoncfg->{TMPDIR}/$assemblyId.gapfill.txt`;
			close (GAPFILL);
			print <<END;
		<script>
		parent.loadingHide();
		parent.closeDialog();
		parent.openDialog('$commoncfg->{TMPURL}/$assemblyId.gapfill.html');
		</script>
END
		}	
	}
}
else
{
	print <<END;
<script>
	parent.loadingHide();
	parent.errorPop("Please give an assembly Id!");
</script>	
END
}
