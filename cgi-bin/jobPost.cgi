#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use LWP::Simple qw/getstore/;
use JSON::XS; #JSON::XS is recommended to be installed for handling JSON string of big size 
use Bio::SeqIO;
use DBI;
use lib "lib/";
use lib "lib/pangu";
use pangu;
use user;
use config;
use userCookie;

my $userCookie = new userCookie;
my $userId = (cookie('cid')) ? $userCookie->checkCookie(cookie('cid')) : 0;
exit if (!$userId);

my $user = new user;
my $userDetail = $user->getAllFieldsWithUserId($userId);
my $userName = $userDetail->{"userName"};

my $commoncfg = readConfig("main.conf");
my $config = new config;
my $JOBURL = $config->getFieldValueWithFieldName("JOBURL");
my $JOBREALTIME = $config->getFieldValueWithFieldName("JOBREALTIME");
my $json = JSON::XS->new->allow_nonref;

my $blastn = 'blast+/bin/blastn';
my $jobdir = $commoncfg->{JOBDIR};
my $shortCutoff = param('shortCutoff') || 10000; #short sequences
my $identity = param('identity') || 95; #end overlap sequences
my $minOverlap =  param('minOverlap') || 500; #end overlap sequences
my @vectorId =  param('vectorId'); #vector
my $tagMatchIdentity = 100; #identities for tags
my $besMatchIdentity = 98; #identities for tags
my $tagMatchPercent = 0.80; #total percent of matched tags
my $minCloneTagNumber = 5;
my $format = 'fasta';
my @jobId = param('jobId');
my $goodOnly = param ('goodOnly') || '';
print header;
my $pid = fork();
if ($pid) {
	print <<END;
<script>
parent.closeDialog();
parent.refresh("menu");
parent.informationPop("It's running! This processing might take a while.");
</script>	
END
}
elsif($pid == 0){
 	close (STDOUT);
	#connect to the mysql server
	my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});

	my $vector = "$commoncfg->{TMPDIR}/$$.vector";
	my $vectorLength;
	if(@vectorId)
	{
		open (VECTOR,">>$vector") or die "can't open file: $vector";
		for(@vectorId)
		{
			unless ($_)
			{
				open (DEFAULTVECTOR,"$commoncfg->{VECTOR}") or die "can't open file: $commoncfg->{VECTOR}";
				while(<DEFAULTVECTOR>)
				{
					print VECTOR $_;
				}
				close(DEFAULTVECTOR);
				$vectorLength->{'pAGIBAC1.HindIII'} = $commoncfg->{VECTORLENGTH};
				next;
			}
			my $vectorSeq = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
			$vectorSeq->execute($_);
			my @vectorSeq=$vectorSeq->fetchrow_array();
			$vectorSeq[8] =~ s/^>.*\n//g;
			$vectorSeq[8] =~ s/[^a-zA-Z0-9]//g;
			print VECTOR ">$vectorSeq[0] $vectorSeq[2]\n$vectorSeq[8]\n";
			$vectorLength->{$vectorSeq[0]} = length ($vectorSeq[8]);
		}
		close(VECTOR);
	}
	
	foreach my $input (sort {$a <=> $b} @jobId)
	{
		my @allSequenceId;
		my $outdir = "$jobdir/$input";
		unless (-e $outdir)
		{
			mkdir $outdir;
		}
		if ($JOBREALTIME)
		{
			getstore("$JOBURL/$input/contents/data/polished_assembly.fasta.gz","$outdir/$input.fasta.gz");
		}
		else
		{
			next if(!-e "$commoncfg->{POLISHED}/$input.fasta.gz");
			`cp -fru $commoncfg->{POLISHED}/$input.fasta.gz $outdir`;
		}
		`gunzip -f $outdir/$input.fasta.gz`;
		my $inputSequence = "$outdir/$input.fasta";
		
		##what to do if the sequences of this job have been used for assembly??? --- keep or delete them???
		#clean up
		my $deleteSequence=$dbh->do("DELETE FROM matrix WHERE x = $input AND container LIKE 'sequence'");
		unlink("$outdir/$input.clean") if (-e "$outdir/$input.clean");
		unlink("$outdir/$input.gapped") if (-e "$outdir/$input.gapped");
		unlink("$outdir/$input.short") if (-e "$outdir/$input.short");
		unlink("$outdir/$input.noVector") if (-e "$outdir/$input.noVector");
		unlink("$outdir/$input.vector") if (-e "$outdir/$input.vector");
		unlink("$outdir/$input.partial") if (-e "$outdir/$input.partial");
		unlink("$outdir/$input.aln") if (-e "$outdir/$input.aln");
		unlink("$outdir/$input.log") if (-e "$outdir/$input.log");
		my $updateJobToRun=$dbh->do("UPDATE matrix SET barcode = '-1' WHERE container LIKE 'job' AND name LIKE '$input'");
		my $count;
		$count->{"Assembled"} = 0;
		$count->{"Circularized"} = 0;
		$count->{"Gapped"} = 0;
		$count->{"Insert"} = 0;
		$count->{"NoVector"} = 0;
		$count->{"Partial"} = 0;
		$count->{"SHORT"} = 0;
		$count->{"Vector"} = 0;
		$count->{"Mixer"} = 0;
		open (CLEAN,">$outdir/$input.clean") or die "can't open file: $outdir/$input.clean";
		open (LOG,">$outdir/$input.log") or die "can't open file: $outdir/$input.log";

		#read sequences
		my $in = Bio::SeqIO->new(-file => $inputSequence,
								-format => $format);
		while ( my $seq = $in->next_seq() )
		{
			$count->{"Assembled"}++;
			my $seqLength = $seq->length();
			my $sequence;
			#raw sequences (type-0)
			my $seqDetails;
			$seqDetails->{'id'} = $seq->id;
			$seqDetails->{'description'} = $seq->desc || '';
			$seqDetails->{'sequence'} = $seq->seq;
			$seqDetails->{'sequence'} =~ tr/a-zA-Z/N/c; #replace nonword characters.
			$seqDetails->{'gapList'} = '';
			$seqDetails = $json->encode($seqDetails);
			my $insertSequence=$dbh->prepare("INSERT INTO matrix VALUES ('', 'sequence', ? , 0, ?, ?, 0, 0, ?, ?, NOW())");
			$insertSequence->execute($seq->id(),$input,$seqLength,$seqDetails,$userName);
			my $rawSequenceId = $dbh->{mysql_insertid};
			push @allSequenceId,$rawSequenceId;


			if($seqLength < $shortCutoff)
			{
				#ignore short sequences (type-8)
				$insertSequence=$dbh->prepare("INSERT INTO matrix VALUES ('', 'sequence', '', 8, ?, ?, 0, 0, ?, ?, NOW())");
				$insertSequence->execute($input,$seqLength,$seqDetails,$userName);
				my $shortSequenceId = $dbh->{mysql_insertid};
				push @allSequenceId,$shortSequenceId;
				$sequence = ">" . $seq->id() . "\n" . $seq->seq() ."\n";
				open (SHORT,">>$outdir/$input.short") or die "can't open file: $outdir/$input.short";
				print SHORT $sequence;
				close(SHORT);
				print LOG $seq->id() . "\t$seqLength\tSHORT\n";
				$count->{"SHORT"}++;
				next;
			}
			else
			{
				open (UNITIG,">$commoncfg->{TMPDIR}/UNITIG$$") or die "can't open $commoncfg->{TMPDIR}/UNITIG$$: $!";
				print UNITIG ">" . $seq->id() . "\n" . $seq->seq() ."\n";
				close (UNITIG);
				print LOG $seq->id() . "\t$seqLength\tLONG(>$shortCutoff)\n";
			}

			#break unitigs by vector
			my %breakPoint;
			$breakPoint{1} = "SequenceEnd";
			$breakPoint{$seqLength} = "SequenceEnd";
			open (ALN,">>$outdir/$input.aln") or die "can't open file: $outdir/$input.aln";
			open (CMD,"$blastn -query $commoncfg->{TMPDIR}/UNITIG$$ -subject $vector -dust no -evalue 1e-200 -outfmt 6 |") or die "can't open CMD: $!";
			while(<CMD>)
			{
				print ALN $_;
				/^#/ and next;
				my @vectorHit = split("\t",$_);
				if($vectorHit[8] < $vectorHit[9])
				{
					if($vectorHit[8] == 1)
					{
						if($vectorHit[6]-1 > 1)
						{
							$breakPoint{$vectorHit[6]} = "VectorEnd";
							$breakPoint{$vectorHit[6]-1} = "InsertEnd";
						}
						else
						{
							$breakPoint{1} = "VectorEnd"; #for avoiding vector happens to be at the end of raw sequences;
						}
					}
					if($vectorHit[9] == $vectorLength->{$vectorHit[1]})
					{
						if($vectorHit[7]+1 < $seqLength)
						{
							$breakPoint{$vectorHit[7]} = "VectorEnd";
							$breakPoint{$vectorHit[7]+1} = "InsertEnd";
						}
						else
						{
							$breakPoint{$seqLength} = "VectorEnd";#for avoiding vector happens to be at the end of raw sequences
						}
					}
				}
				else
				{
					if($vectorHit[9] == 1)
					{
						if($vectorHit[7]+1 < $seqLength)
						{
							$breakPoint{$vectorHit[7]} = "VectorEnd";
							$breakPoint{$vectorHit[7]+1} = "InsertEnd";
						}
						else
						{
							$breakPoint{$seqLength} = "VectorEnd";#for avoiding vector happens to be at the end of raw sequences
						}
					}
					if($vectorHit[8] == $vectorLength->{$vectorHit[1]})
					{
						if($vectorHit[6]-1 > 1)
						{
							$breakPoint{$vectorHit[6]} = "VectorEnd";
							$breakPoint{$vectorHit[6]-1} = "InsertEnd";
						}
						else
						{
							$breakPoint{1} = "VectorEnd"; #for avoiding vector happens to be at the end of raw sequences;
						}
					}
				}
			}
			close(CMD);
			close(ALN);		

			#determine sequence piece type
			my $pieceNumber = 0;
			my $pieceLeft;
			my $pieceLeftPosition;
			my $subSeqALeft;
			my $subSeqARight;
			my $subSeqBLeft;
			my $subSeqBRight;
			for(sort {$a <=> $b} keys %breakPoint)
			{
				if($pieceLeft)
				{
					$pieceNumber++;
					$seqLength = $_ - $pieceLeftPosition + 1;
					my $sequenceDetails;
					$sequenceDetails->{'id'} = $seq->id;
					$sequenceDetails->{'description'} = "($pieceLeftPosition-$_)";
					$sequenceDetails->{'sequence'} = $seq->subseq($pieceLeftPosition,$_);
					$sequenceDetails->{'sequence'} =~ tr/a-zA-Z/N/c; #replace nonword characters.
					$sequenceDetails->{'gapList'} = '';
					$sequenceDetails = $json->encode($sequenceDetails);
					$sequence = ">" . $seq->id() . " ($pieceLeftPosition-$_)\n" . $seq->subseq($pieceLeftPosition,$_) ."\n";
					if ($pieceLeft eq "SequenceEnd")
					{
						if($breakPoint{$_} eq "SequenceEnd")
						{
							#no vector(type-3)
							$insertSequence=$dbh->prepare("INSERT INTO matrix VALUES ('', 'sequence', '', 3, ?, ?, 0, 0, ?, ?, NOW())");
							$insertSequence->execute($input,$seqLength,$sequenceDetails,$userName);
							my $novectorSequenceId = $dbh->{mysql_insertid};
							push @allSequenceId,$novectorSequenceId;

							open (NOVECTOR,">>$outdir/$input.noVector") or die "can't open file: $outdir/$input.noVector";
							print NOVECTOR $sequence;
							close(NOVECTOR);
							print LOG $seq->id() . "\t$pieceLeftPosition\t$_\tNoVector\tSequenceEnd-SequenceEnd\n";
							$count->{"NoVector"}++;
						}
						elsif($breakPoint{$_} eq "VectorEnd")
						{
							#vector or mixer(type-6)
							$insertSequence=$dbh->prepare("INSERT INTO matrix VALUES ('', 'sequence', '', 6, ?, ?, 0, 0, ?, ?, NOW())");
							$insertSequence->execute($input,$seqLength,$sequenceDetails,$userName);
							my $vectorSequenceId = $dbh->{mysql_insertid};
							push @allSequenceId,$vectorSequenceId;

							open (VECTOR,">>$outdir/$input.vector") or die "can't open file: $outdir/$input.vector";
							print VECTOR $sequence;
							close(VECTOR);					
							print LOG $seq->id() . "\t$pieceLeftPosition\t$_\tVector\tSequenceEnd-VectorEnd\n";
							$count->{"Vector"}++;
						}
						else # "InsertEnd"
						{
							#left end (sub-sequence one)
							open (SUBSEQONE,">$commoncfg->{TMPDIR}/subSeqA$$") or die "can't open $commoncfg->{TMPDIR}/subSeqA$$: $!";
							print SUBSEQONE ">" . $seq->id() . "-subSeqA\n" . $seq->subseq($pieceLeftPosition,$_) . "\n";
							close (SUBSEQONE);
							$subSeqALeft = $pieceLeftPosition;
							$subSeqARight = $_;
							print LOG $seq->id() . "\t$pieceLeftPosition\t$_\tSubSeqA\tSequenceEnd-InsertEnd\n";
						}
					}
					elsif($pieceLeft eq "VectorEnd")
					{
						if($breakPoint{$_} eq "SequenceEnd")
						{
							#vector or mixer(type-6)
							$insertSequence=$dbh->prepare("INSERT INTO matrix VALUES ('', 'sequence', '', 6, ?, ?, 0, 0, ?, ?, NOW())");
							$insertSequence->execute($input,$seqLength,$sequenceDetails,$userName);
							my $vectorSequenceId = $dbh->{mysql_insertid};
							push @allSequenceId,$vectorSequenceId;

							open (VECTOR,">>$outdir/$input.vector") or die "can't open file: $outdir/$input.vector";
							print VECTOR $sequence;
							close(VECTOR);					
							print LOG $seq->id() . "\t$pieceLeftPosition\t$_\tVector\tVectorEnd-SequenceEnd\n";
							$count->{"Vector"}++;
						}
						elsif($breakPoint{$_} eq "VectorEnd")
						{
							#vector(type-6)
							$insertSequence=$dbh->prepare("INSERT INTO matrix VALUES ('', 'sequence', '', 6, ?, ?, 0, 0, ?, ?, NOW())");
							$insertSequence->execute($input,$seqLength,$sequenceDetails,$userName);
							my $vectorSequenceId = $dbh->{mysql_insertid};
							push @allSequenceId,$vectorSequenceId;

							open (VECTOR,">>$outdir/$input.vector") or die "can't open file: $outdir/$input.vector";
							print VECTOR $sequence;
							close(VECTOR);					
							print LOG $seq->id() . "\t$pieceLeftPosition\t$_\tVector\tVectorEnd-VectorEnd\n";
							$count->{"Vector"}++;
						}
						else # "InsertEnd"
						{
							#mixer (type-7)
							$insertSequence=$dbh->prepare("INSERT INTO matrix VALUES ('', 'sequence', '', 7, ?, ?, 0, 0, ?, ?, NOW())");
							$insertSequence->execute($input,$seqLength,$sequenceDetails,$userName);
							my $mixerSequenceId = $dbh->{mysql_insertid};
							push @allSequenceId,$mixerSequenceId;

							open (PARTIAL,">>$outdir/$input.partial") or die "can't open file: $outdir/$input.partial";
							print PARTIAL $sequence;
							close(PARTIAL);					
							print LOG $seq->id() . "\t$pieceLeftPosition\t$_\tMixer\tVectorEnd-InsertEnd\n";
							$count->{"Mixer"}++;
						}
					}
					else # "InsertEnd"
					{
						if($breakPoint{$_} eq "SequenceEnd")
						{
							#right end (sub-sequence two)
							open (SUBSEQTWO,">$commoncfg->{TMPDIR}/subSeqB$$") or die "can't open $commoncfg->{TMPDIR}/subSeqB$$: $!";
							print SUBSEQTWO ">" . $seq->id() . "-subSeqB\n" . $seq->subseq($pieceLeftPosition,$_) . "\n";
							close (SUBSEQTWO);
							$subSeqBLeft = $pieceLeftPosition;
							$subSeqBRight = $_;
							print LOG $seq->id() . "\t$pieceLeftPosition\t$_\tSubSeqB\tInsertEnd-SequenceEnd\n";
						}
						elsif($breakPoint{$_} eq "VectorEnd")
						{
							#mixer (type-7)
							$insertSequence=$dbh->prepare("INSERT INTO matrix VALUES ('', 'sequence', '', 7, ?, ?, 0, 0, ?, ?, NOW())");
							$insertSequence->execute($input,$seqLength,$sequenceDetails,$userName);
							my $mixerSequenceId = $dbh->{mysql_insertid};
							push @allSequenceId,$mixerSequenceId;

							open (PARTIAL,">>$outdir/$input.partial") or die "can't open file: $outdir/$input.partial";
							print PARTIAL $sequence;
							close(PARTIAL);
							print LOG $seq->id() . "\t$pieceLeftPosition\t$_\tMixer\tInsertEnd-VectorEnd\n";
							$count->{"Mixer"}++;
						}
						else # "InsertEnd"
						{
							$seqLength = $_ - $pieceLeftPosition + 1;
							if($seqLength < $shortCutoff)
							{
								#ignore short sequences (type-8)
								$insertSequence=$dbh->prepare("INSERT INTO matrix VALUES ('', 'sequence', '', 8, ?, ?, 0, 0, ?, ?, NOW())");
								$insertSequence->execute($input,$seqLength,$sequenceDetails,$userName);
								my $shortSequenceId = $dbh->{mysql_insertid};
								push @allSequenceId,$shortSequenceId;

								open (SHORT,">>$outdir/$input.short") or die "can't open file: $outdir/$input.short";
								print SHORT $sequence;
								close(SHORT);
								print LOG $seq->id() . "\t$seqLength\tSHORT-insert\t$pieceLeftPosition-$_\n";
								$count->{"SHORT"}++;
							}
							else
							{
								#good insert sequence (type-1)
								$insertSequence=$dbh->prepare("INSERT INTO matrix VALUES ('', 'sequence', '', 1, ?, ?, 0, 0, ?, ?, NOW())");
								$insertSequence->execute($input,$seqLength,$sequenceDetails,$userName);
								my $goodSequenceId = $dbh->{mysql_insertid};
								push @allSequenceId,$goodSequenceId;

								print CLEAN $sequence;
								print LOG $seq->id() . "\t$pieceLeftPosition\t$_\tInsert\tInsertEnd-InsertEnd\n";
								$count->{"Insert"}++;
							}
						}
					}
					$pieceLeft = '';
				}
				else
				{
					$pieceLeft = $breakPoint{$_};
					$pieceLeftPosition = $_;
				}
			}
			if (-e "$commoncfg->{TMPDIR}/subSeqA$$")
			{
				if(-e "$commoncfg->{TMPDIR}/subSeqB$$")
				{
					# This script use blast2seq function to count length of overlapping region between 2 subseqs
					my $midHit = 0;
					my $queryPosition = 0;
					my $subjectPosition = 0;
					my $overlapLength = 0;
					my $overlapIdentities = 0;
					open (ALN,">>$outdir/$input.aln") or die "can't open file: $outdir/$input.aln";
					open (CMD,"$blastn -query $commoncfg->{TMPDIR}/subSeqA$$ -subject $commoncfg->{TMPDIR}/subSeqB$$ -dust no -evalue 1e-200 -perc_identity $identity |") or die "can't open CMD: $!";
					while(<CMD>)
					{
						print ALN $_;
						last if ($subjectPosition > 0);
						/Strand=Plus\/Minus/ and last;
						if(/Sbjct  (\d*)/)
						{
							$subjectPosition = $1 if ($queryPosition > 0);
						}
						if(/Query  (\d*)/)
						{
							$queryPosition = $1 if ($1 > $midHit);
						}
						if(/Identities = (\d+)\/(\d+)/)
						{
							$midHit = int ($1/2);
							$overlapLength = $1;
							$overlapIdentities = $1/$2;
							last if ($1 < $minOverlap);
						}		
					}
					close(CMD);
					close(ALN);	
					if($queryPosition > 0)
					{
						my $subSeqBEnd = $subSeqBLeft + $subjectPosition - 2;
						my $subSeqAStart = $subSeqALeft + $queryPosition - 1;
						if($subjectPosition == 1) #only subSeqA is used because alignment of the subSeqB start from 1.
						{
							#good Circularized sequence (type-2)
							#this is very special case: alignment looks like below
							#             1234567...
							# subSeqB     ----------
							#             ||||||||||
							# subSeqA -------------------------
							#             ^cut here
						
							$seqLength = $subSeqARight - $subSeqAStart + 1;
							my $sequenceDetails;
							$sequenceDetails->{'id'} = $seq->id;
							$sequenceDetails->{'description'} = "($subSeqAStart-$subSeqARight,Overlap:$overlapLength-$overlapIdentities)";
							$sequenceDetails->{'sequence'} = $seq->subseq($subSeqAStart, $subSeqARight);
							$sequenceDetails->{'sequence'} =~ tr/a-zA-Z/N/c; #replace nonword characters.
							$sequenceDetails->{'gapList'} = "";
							$sequenceDetails = $json->encode($sequenceDetails);
							$insertSequence=$dbh->prepare("INSERT INTO matrix VALUES ('', 'sequence', '', 2, ?, ?, 0, 0, ?, ?, NOW())");
							$insertSequence->execute($input,$seqLength,$sequenceDetails,$userName);
							my $goodSequenceId = $dbh->{mysql_insertid};
							push @allSequenceId,$goodSequenceId;

							$sequence = ">" . $seq->id() ." ($subSeqAStart-$subSeqARight,Overlap:$overlapLength-$overlapIdentities)\n" . $seq->subseq($subSeqAStart, $subSeqARight) ."\n";
							print CLEAN $sequence;
							print LOG $seq->id() . "\t$subSeqAStart\t$subSeqARight\tCircularized\tsubSeqA only(Overlap:$overlapLength,$overlapIdentities)\n";
							$count->{"Circularized"}++;
						}
						else
						{
							#good Circularized sequence (type-2)
							#normal alignment looks like below
							#         1234567...
							# subSeqB --------------
							#             ||||||||||
							# subSeqA     -------------------------
							#                  ^cut here
							$seqLength = $subSeqBEnd - $subSeqBLeft + 1 + $subSeqARight - $subSeqAStart + 1;
							my $sequenceDetails;
							$sequenceDetails->{'id'} = $seq->id;
							$sequenceDetails->{'description'} = "($subSeqBLeft-$subSeqBEnd,$subSeqAStart-$subSeqARight,Overlap:$overlapLength-$overlapIdentities)";
							$sequenceDetails->{'sequence'} = $seq->subseq($subSeqBLeft, $subSeqBEnd) ."\n". $seq->subseq($subSeqAStart, $subSeqARight);
							$sequenceDetails->{'sequence'} =~ tr/a-zA-Z/N/c; #replace nonword characters.
							$sequenceDetails->{'gapList'} = "";
							$sequenceDetails = $json->encode($sequenceDetails);
							$insertSequence=$dbh->prepare("INSERT INTO matrix VALUES ('', 'sequence', '', 2, ?, ?, 0, 0, ?, ?, NOW())");
							$insertSequence->execute($input,$seqLength,$sequenceDetails,$userName);
							my $goodSequenceId = $dbh->{mysql_insertid};
							push @allSequenceId,$goodSequenceId;

							$sequence = ">" . $seq->id() ." ($subSeqBLeft-$subSeqBEnd,$subSeqAStart-$subSeqARight,Overlap:$overlapLength-$overlapIdentities)\n" . $seq->subseq($subSeqBLeft, $subSeqBEnd) ."\n". $seq->subseq($subSeqAStart, $subSeqARight) ."\n";
							print CLEAN $sequence;
							print LOG $seq->id() . "\t$subSeqBLeft-$subSeqBEnd\t$subSeqAStart-$subSeqARight\tCircularized\tsubSeqB-subSeqA(Overlap:$overlapLength,$overlapIdentities)\n";
							$count->{"Circularized"}++;
						}
					}
					else
					{
						if($pieceNumber > 99) #to be polished 
						{
							open (PARTIAL,">>$outdir/$input.partial") or die "can't open file: $outdir/$input.partial";
							#partial sequence;(type-5)
							$seqLength = $subSeqARight - $subSeqALeft + 1;
							my $sequenceDetailsA;
							$sequenceDetailsA->{'id'} = $seq->id;
							$sequenceDetailsA->{'description'} = "($subSeqALeft-$subSeqARight)";
							$sequenceDetailsA->{'sequence'} = $seq->subseq($subSeqALeft,$subSeqARight);
							$sequenceDetailsA->{'sequence'} =~ tr/a-zA-Z/N/c; #replace nonword characters.
							$sequenceDetailsA->{'gapList'} = '';
							$sequenceDetailsA = $json->encode($sequenceDetailsA);
							$insertSequence=$dbh->prepare("INSERT INTO matrix VALUES ('', 'sequence', '', 5, ?, ?, 0, 0, ?, ?, NOW())");
							$insertSequence->execute($input,$seqLength,$sequenceDetailsA,$userName);
							my $partialSequenceIdA = $dbh->{mysql_insertid};
							push @allSequenceId,$partialSequenceIdA;

							$sequence = ">" . $seq->id() . " ($subSeqALeft-$subSeqARight)\n" . $seq->subseq($subSeqALeft,$subSeqARight) ."\n";
							print PARTIAL $sequence;

							#partial sequence;(type-5)
							$seqLength = $subSeqBRight - $subSeqBLeft + 1;
							my $sequenceDetailsB;
							$sequenceDetailsB->{'id'} = $seq->id;
							$sequenceDetailsB->{'description'} = "($subSeqBLeft-$subSeqBRight)";
							$sequenceDetailsB->{'sequence'} = $seq->subseq($subSeqBLeft,$subSeqBRight);
							$sequenceDetailsB->{'sequence'} =~ tr/a-zA-Z/N/c; #replace nonword characters.
							$sequenceDetailsB->{'gapList'} = '';
							$sequenceDetailsB = $json->encode($sequenceDetailsB);
							$insertSequence=$dbh->prepare("INSERT INTO matrix VALUES ('', 'sequence', '', 5, ?, ?, 0, 0, ?, ?, NOW())");
							$insertSequence->execute($input,$seqLength,$sequenceDetailsB,$userName);
							my $partialSequenceIdB = $dbh->{mysql_insertid};
							push @allSequenceId,$partialSequenceIdB;

							$sequence = ">" . $seq->id() . " ($subSeqBLeft-$subSeqBRight)\n" . $seq->subseq($subSeqBLeft,$subSeqBRight) ."\n";
							print PARTIAL $sequence;
							close(PARTIAL);
							print LOG $seq->id() . "\t$subSeqALeft\t$subSeqARight\tPartial\tsubSeqA\n";
							print LOG $seq->id() . "\t$subSeqBLeft\t$subSeqBRight\tPartial\tsubSeqB\n";
							$count->{"Partial"}++;
							$count->{"Partial"}++;
						}
						else
						{
							#gapped sequence(type-4)
							$seqLength = $subSeqBRight - $subSeqBLeft + 1 + 100 + $subSeqARight - $subSeqALeft + 1;
							my $sequenceDetails;
							$sequenceDetails->{'id'} = $seq->id;
							$sequenceDetails->{'description'} = "(subSeqB:$subSeqBLeft-$subSeqBRight,100Ns,subSeqA:$subSeqALeft-$subSeqARight)";
							$sequenceDetails->{'sequence'} = $seq->subseq($subSeqBLeft,$subSeqBRight) . "N" x 100 . $seq->subseq($subSeqALeft,$subSeqARight);
							$sequenceDetails->{'sequence'} =~ tr/a-zA-Z/N/c; #replace nonword characters.
							my $gapStart = $subSeqBRight - $subSeqBLeft + 2;
							my $gapEnd = $gapStart + 99;
							$sequenceDetails->{'gapList'} = "$gapStart-$gapEnd";
							$sequenceDetails = $json->encode($sequenceDetails);

							if($seqLength < $shortCutoff)
							{
								#ignore short sequences (type-8)
								$insertSequence=$dbh->prepare("INSERT INTO matrix VALUES ('', 'sequence', '', 8, ?, ?, 0, 0, ?, ?, NOW())");
								$insertSequence->execute($input,$seqLength,$sequenceDetails,$userName);
								my $shortSequenceId = $dbh->{mysql_insertid};
								push @allSequenceId,$shortSequenceId;

								$sequence = ">" . $seq->id(). " (subSeqB:$subSeqBLeft-$subSeqBRight,100Ns,subSeqA:$subSeqALeft-$subSeqARight)\n" .
											$seq->subseq($subSeqBLeft,$subSeqBRight) . "N" x 100 . $seq->subseq($subSeqALeft,$subSeqARight) ."\n";
								open (SHORT,">>$outdir/$input.short") or die "can't open file: $outdir/$input.short";
								print SHORT $sequence;
								close(SHORT);
								print LOG $seq->id() . "\t$subSeqBLeft-$subSeqBRight\t$subSeqALeft-$subSeqARight\tSHORT-Gapped\tsubSeqB-subSeqA\n";
								$count->{"SHORT"}++;
							}
							else
							{
								$insertSequence=$dbh->prepare("INSERT INTO matrix VALUES ('', 'sequence', '', 4, ?, ?, 0, 0, ?, ?, NOW())");
								$insertSequence->execute($input,$seqLength,$sequenceDetails,$userName);
								my $gapSequenceId = $dbh->{mysql_insertid};
								push @allSequenceId,$gapSequenceId;

								$sequence = ">" . $seq->id(). " (subSeqB:$subSeqBLeft-$subSeqBRight,100Ns,subSeqA:$subSeqALeft-$subSeqARight)\n" .
											$seq->subseq($subSeqBLeft,$subSeqBRight) . "N" x 100 . $seq->subseq($subSeqALeft,$subSeqARight) ."\n";
								open (GAPPED,">>$outdir/$input.gapped") or die "can't open file: $outdir/$input.gapped";
								print GAPPED $sequence;
								close(GAPPED);
								print LOG $seq->id() . "\t$subSeqBLeft-$subSeqBRight\t$subSeqALeft-$subSeqARight\tGapped\tsubSeqB-subSeqA\n";
								$count->{"Gapped"}++;
							}
						}
					}			
				}
				else
				{
					#partial sequence;(type-5)
					$seqLength = $subSeqARight - $subSeqALeft + 1;
					my $sequenceDetails;
					$sequenceDetails->{'id'} = $seq->id;
					$sequenceDetails->{'description'} = "($subSeqALeft-$subSeqARight)";
					$sequenceDetails->{'sequence'} = $seq->subseq($subSeqALeft,$subSeqARight);
					$sequenceDetails->{'sequence'} =~ tr/a-zA-Z/N/c; #replace nonword characters.
					$sequenceDetails->{'gapList'} = '';
					$sequenceDetails = $json->encode($sequenceDetails);
					$insertSequence=$dbh->prepare("INSERT INTO matrix VALUES ('', 'sequence', '', 5, ?, ?, 0, 0, ?, ?, NOW())");
					$insertSequence->execute($input,$seqLength,$sequenceDetails,$userName);
					my $partialSequenceId = $dbh->{mysql_insertid};
					push @allSequenceId,$partialSequenceId;

					$sequence = ">" . $seq->id() . " ($subSeqALeft-$subSeqARight)\n" . $seq->subseq($subSeqALeft,$subSeqARight) ."\n";
					open (PARTIAL,">>$outdir/$input.partial") or die "can't open file: $outdir/$input.partial";
					print PARTIAL $sequence;
					close(PARTIAL);
					print LOG $seq->id() . "\t$subSeqALeft\t$subSeqARight\tPartial\tsubSeqA\n";
					$count->{"Partial"}++;
				}
			}
			else
			{
				if(-e "subSeqB$$")
				{
					#partial sequence;(type-5)
					$seqLength = $subSeqBRight - $subSeqBLeft + 1;
					my $sequenceDetails;
					$sequenceDetails->{'id'} = $seq->id;
					$sequenceDetails->{'description'} = "($subSeqBLeft-$subSeqBRight)";
					$sequenceDetails->{'sequence'} = $seq->subseq($subSeqBLeft,$subSeqBRight);
					$sequenceDetails->{'sequence'} =~ tr/a-zA-Z/N/c; #replace nonword characters.
					$sequenceDetails->{'gapList'} = '';
					$sequenceDetails = $json->encode($sequenceDetails);
					$insertSequence=$dbh->prepare("INSERT INTO matrix VALUES ('', 'sequence', '', 5, ?, ?, 0, 0, ?, ?, NOW())");
					$insertSequence->execute($input,$seqLength,$sequenceDetails,$userName);
					my $partialSequenceId = $dbh->{mysql_insertid};
					push @allSequenceId,$partialSequenceId;

					$sequence = ">" . $seq->id() . " ($subSeqBLeft-$subSeqBRight)\n" . $seq->subseq($subSeqBLeft,$subSeqBRight) ."\n";
					open (PARTIAL,">>$outdir/$input.partial") or die "can't open file: $outdir/$input.partial";
					print PARTIAL $sequence;
					close(PARTIAL);
					print LOG $seq->id() . "\t$subSeqBLeft\t$subSeqBRight\tPartial\tsubSeqB\n";
					$count->{"Partial"}++;
				}
			}
			unlink("$commoncfg->{TMPDIR}/subSeqA$$");
			unlink("$commoncfg->{TMPDIR}/subSeqB$$");
		}
		unlink("$commoncfg->{TMPDIR}/UNITIG$$");
		close(CLEAN);		
		print LOG "=" x 80;
		print LOG "\n";

		my $updateJob=$dbh->prepare("UPDATE matrix SET o = ?, x = ?, y = ?, z = ?, creator = ? WHERE container LIKE 'job' AND name LIKE ?");
		$updateJob->execute($count->{"Assembled"},$count->{"Circularized"} + $count->{"Insert"},$count->{"NoVector"},$count->{"Gapped"},$userName,$input);

		open (ALLLOG,">>$jobdir/jobs.log") or die "can't open file: $jobdir/jobs.log";
		print ALLLOG "$input\tGood:";
		print ALLLOG $count->{"Circularized"} + $count->{"Insert"} . "\t";
		for (sort keys %$count)
		{
			print LOG $_ . ":" . $count->{$_} . "\n";
			print ALLLOG "\t" . $_ . ":" . $count->{$_} ;
		}
		print ALLLOG "\n";
		close(ALLLOG);		

		foreach (@allSequenceId)
		{
			my $getSequences = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
			$getSequences->execute($_);
			while(my @getSequences = $getSequences->fetchrow_array())
			{
				my $sequenceDetails = decode_json $getSequences[8];
				$sequenceDetails->{'sequence'} = '' unless (exists $sequenceDetails->{'sequence'});
				my $seqDir;
				for (my $position = 0; $position < length($getSequences[4]); $position += 2)
				{
					$seqDir .= "/j". substr($getSequences[4],$position,2);
					until (-e "$commoncfg->{DATADIR}/sequences$seqDir")
					{
						mkdir "$commoncfg->{DATADIR}/sequences$seqDir";
					}
				}

				my $seqFile;
				for (my $position = 0; $position < length($getSequences[0]); $position += 2)
				{
					$seqFile .= "/s". substr($getSequences[0],$position,2);
					if (length ($seqFile) % 4 == 0)
					{
						until (-e "$commoncfg->{DATADIR}/sequences$seqDir$seqFile")
						{
							mkdir "$commoncfg->{DATADIR}/sequences$seqDir$seqFile";
						}
					}
				}

				until (-e "$commoncfg->{DATADIR}/sequences$seqDir$seqFile.fa")
				{
					open (SEQ,">$commoncfg->{DATADIR}/sequences$seqDir$seqFile.fa") or die "can't open file: $commoncfg->{DATADIR}/sequences$seqDir$seqFile.fa";
					print SEQ ">$getSequences[0]\n";
					print SEQ &multiLineSeq($sequenceDetails->{'sequence'},80);
					close(SEQ);
				}

				$sequenceDetails->{'sequence'} = "sequences$seqDir$seqFile.fa";
				my $seqDetailsEncoded = $json->encode($sequenceDetails);
				my $updateSequence = $dbh->prepare("UPDATE matrix SET note = ? WHERE id = ?");
				$updateSequence->execute($seqDetailsEncoded,$getSequences[0]);
			}
		}

		#associate BAC and sequence
		my $assignedSequenceNumber = 0;
		my $jobToPool = $dbh->prepare("SELECT parent FROM link WHERE type LIKE 'poolJob' AND child = ?");
		$jobToPool->execute($input);
		while(my @jobToPool = $jobToPool->fetchrow_array())
		{
			#get clone list
			my $cloneTagNumber;
			my $cloneBesNumber;
			my $tagTotalNumber = 0;
			my $besTotalNumber = 0;
			open (TAG,">$commoncfg->{TMPDIR}/$input.$$.tag") or die "can't open file: $commoncfg->{TMPDIR}/$input.$$.tag";
			open (BES,">$commoncfg->{TMPDIR}/$input.$$.bes") or die "can't open file: $commoncfg->{TMPDIR}/$input.$$.bes";
			my $poolToClone = $dbh->prepare("SELECT child FROM link WHERE type LIKE 'poolClone' AND parent = ?");
			$poolToClone->execute($jobToPool[0]);
			while (my @poolToClone = $poolToClone->fetchrow_array())
			{
				#get tag list
				my $cloneToTag = $dbh->prepare("SELECT matrix.* FROM matrix,clones WHERE matrix.container LIKE 'tag' AND (matrix.name LIKE clones.name OR matrix.name LIKE clones.origin) AND clones.name LIKE ? ORDER BY matrix.o");
				$cloneToTag->execute($poolToClone[0]);
				while (my @cloneToTag = $cloneToTag->fetchrow_array())
				{
					$cloneTagNumber->{$cloneToTag[2]}++;
					$tagTotalNumber++;
					print TAG ">$cloneToTag[2].$cloneToTag[3]\n$cloneToTag[8]\n";						
				}			
				#get BES list
				my $cloneToBes = $dbh->prepare("SELECT matrix.* FROM matrix,clones WHERE matrix.container LIKE 'sequence' AND matrix.o = 98 AND (matrix.name LIKE clones.name OR matrix.name LIKE clones.origin) AND clones.name LIKE ? ORDER BY matrix.z");
				$cloneToBes->execute($poolToClone[0]);
				while (my @cloneToBes = $cloneToBes->fetchrow_array())
				{
					$cloneBesNumber->{$cloneToBes[2]}++;
					$besTotalNumber++;
					my $sequenceDetails = decode_json $cloneToBes[8];
					my $sequence = 'ERROR: NO SEQUENCE FOUND! PLEASE CONTACT YOUR ADMINISTRATOR.';
					my $in = Bio::SeqIO->new(-file => "$commoncfg->{DATADIR}/$sequenceDetails->{'sequence'}",
											-format => 'Fasta');
					while ( my $seq = $in->next_seq() )
					{
						$sequence = $seq->seq;
					}
					print BES ">$cloneToBes[2].$cloneToBes[6]\n$sequence\n";						
				}			
			}
			close (BES);
			close (TAG);
			my $jobToSequence;
			if($goodOnly)
			{
				$jobToSequence = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND (o = 1 OR o = 2) AND x = ?"); #get circularized or inserted sequences only
			}
			else
			{
				$jobToSequence = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND (o = 1 OR o = 2 OR o = 3 OR o = 4) AND x = ?"); #get circularized or inserted sequences, non-vector or gapped sequences
			}
			$jobToSequence->execute($input);
			while(my @jobToSequence = $jobToSequence->fetchrow_array() )
			{
				my $bacIdAssigned = 0;
				my $sequenceDetails = decode_json $jobToSequence[8];
				my $sequence = 'ERROR: NO SEQUENCE FOUND! PLEASE CONTACT YOUR ADMINISTRATOR.';
				my $in = Bio::SeqIO->new(-file => "$commoncfg->{DATADIR}/$sequenceDetails->{'sequence'}",
										-format => 'Fasta');
				while ( my $seq = $in->next_seq() )
				{
					$sequence = $seq->seq;
				}
				open (SEQ,">$commoncfg->{TMPDIR}/$jobToSequence[0].seq") or die "can't open file: $commoncfg->{TMPDIR}/$jobToSequence[0].seq";
				print SEQ ">$jobToSequence[0]\n$sequence";
				close (SEQ);
				if($tagTotalNumber > 0)
				{
					my $matchedTagNumber;
					open (CMD,"$blastn -query $commoncfg->{TMPDIR}/$jobToSequence[0].seq -subject $commoncfg->{TMPDIR}/$input.$$.tag -dust no -max_target_seqs 10000 -outfmt 6 |") or die "can't open CMD: $!";
					while(<CMD>)
					{
						my @blastLine = split /\t/, $_;
						next if ($blastLine[2] < $tagMatchIdentity);
						my ($bac,$tagId)  = split /\./,$blastLine[1];
						$matchedTagNumber->{$bac}++;
					}
					close(CMD);
					for (sort { $matchedTagNumber->{$b} <=> $matchedTagNumber->{$a} } keys(%$matchedTagNumber))
					{
						print LOG "$jobToSequence[0]\t$_\t$cloneTagNumber->{$_}\t$matchedTagNumber->{$_}\t";
						if (($matchedTagNumber->{$_}/$cloneTagNumber->{$_} > $tagMatchPercent) && ($cloneTagNumber->{$_} >= $minCloneTagNumber))
						{
							my $updateSequence=$dbh->prepare("UPDATE matrix SET name = ?, z = ?, barcode = 1 WHERE id = ?");
							$updateSequence->execute($_,$matchedTagNumber->{$_},$jobToSequence[0]);
							my $sequencedStatus = 0;
							my $getSequenced = $dbh->prepare("SELECT * FROM clones WHERE name LIKE ? OR origin LIKE ?");
							$getSequenced->execute($_,$_);
							while (my @getSequenced = $getSequenced->fetchrow_array())
							{
								$sequencedStatus = $getSequenced[6];
							}
							if ($jobToSequence[3] == 1 || $jobToSequence[3] == 2)
							{
								my $updateClone=$dbh->prepare("UPDATE clones SET sequenced = ? WHERE name LIKE ? OR origin LIKE ?");
								$updateClone->execute($jobToSequence[3],$_,$_);
								my $updateFpcClone=$dbh->prepare("UPDATE matrix SET x = ? WHERE container LIKE 'fpcClone' AND name LIKE ?");
								$updateFpcClone->execute($jobToSequence[3],$_);
							}
							elsif($jobToSequence[3] == 3)
							{
								if($sequencedStatus != 1 || $sequencedStatus != 2)
								{
									my $updateClone=$dbh->prepare("UPDATE clones SET sequenced = ? WHERE name LIKE ? OR origin LIKE ?");
									$updateClone->execute($jobToSequence[3],$_,$_);
									my $updateFpcClone=$dbh->prepare("UPDATE matrix SET x = ? WHERE container LIKE 'fpcClone' AND name LIKE ?");
									$updateFpcClone->execute($jobToSequence[3],$_);
								}
							}
							elsif($jobToSequence[3] == 4)
							{
								if($sequencedStatus != 1 || $sequencedStatus != 2 || $sequencedStatus != 3)
								{
									my $updateClone=$dbh->prepare("UPDATE clones SET sequenced = ? WHERE name LIKE ? OR origin LIKE ?");
									$updateClone->execute($jobToSequence[3],$_,$_);
									my $updateFpcClone=$dbh->prepare("UPDATE matrix SET x = ? WHERE container LIKE 'fpcClone' AND name LIKE ?");
									$updateFpcClone->execute($jobToSequence[3],$_);
								}
							}
							else
							{
								if($sequencedStatus != 1 || $sequencedStatus != 2 || $sequencedStatus != 3 || $sequencedStatus != 4)
								{
									my $updateClone=$dbh->prepare("UPDATE clones SET sequenced = ? WHERE name LIKE ? OR origin LIKE ?");
									$updateClone->execute($jobToSequence[3],$_,$_);
									my $updateFpcClone=$dbh->prepare("UPDATE matrix SET x = ? WHERE container LIKE 'fpcClone' AND name LIKE ?");
									$updateFpcClone->execute($jobToSequence[3],$_);
								}
							}
							$assignedSequenceNumber++;
							$bacIdAssigned = 1;
							print LOG "Good\n";
							last;
						}
						else
						{
							print LOG "Bad\n";
						}				
					}
				}
				if($besTotalNumber > 0 && $bacIdAssigned < 1)
				{
					my $matchedBesNumber;
					open (CMD,"$blastn -query $commoncfg->{TMPDIR}/$jobToSequence[0].seq -subject $commoncfg->{TMPDIR}/$input.$$.bes -dust no -max_target_seqs 10000 -outfmt 6 |") or die "can't open CMD: $!";
					while(<CMD>)
					{
						my @blastLine = split /\t/, $_;
						next if ($blastLine[2] < $besMatchIdentity);
						next if ($blastLine[6] > 1000 && $blastLine[7] < $jobToSequence[5] - 1000); #determine ends location
						my ($bac,$besDirection)  = split /\./,$blastLine[1];
						$matchedBesNumber->{$bac}++;
					}
					close(CMD);
					for (sort { $matchedBesNumber->{$b} <=> $matchedBesNumber->{$a} } keys(%$matchedBesNumber))
					{
						print LOG "$jobToSequence[0]\t$_\tBES$cloneBesNumber->{$_}\t$matchedBesNumber->{$_}\t";
						if ($matchedBesNumber->{$_} > 0)
						{
							if($matchedBesNumber->{$_} >= $cloneBesNumber->{$_})
							{
								my $updateSequence=$dbh->prepare("UPDATE matrix SET name = ?, z = ?, barcode = 2 WHERE id = ?");
								$updateSequence->execute($_,$matchedBesNumber->{$_},$jobToSequence[0]);
							}
							else
							{
								#this condition needs to be updated
								my $updateSequence=$dbh->prepare("UPDATE matrix SET name = ?, z = ?, barcode = 2 WHERE id = ?");
								$updateSequence->execute($_,$matchedBesNumber->{$_},$jobToSequence[0]);
							}
							my $sequencedStatus = 0;
							my $getSequenced = $dbh->prepare("SELECT * FROM clones WHERE name LIKE ? OR origin LIKE ?");
							$getSequenced->execute($_,$_);
							while (my @getSequenced = $getSequenced->fetchrow_array())
							{
								$sequencedStatus = $getSequenced[6];
							}
							if ($jobToSequence[3] == 1 || $jobToSequence[3] == 2)
							{
								my $updateClone=$dbh->prepare("UPDATE clones SET sequenced = ? WHERE name LIKE ? OR origin LIKE ?");
								$updateClone->execute($jobToSequence[3],$_,$_);
								my $updateFpcClone=$dbh->prepare("UPDATE matrix SET x = ? WHERE container LIKE 'fpcClone' AND name LIKE ?");
								$updateFpcClone->execute($jobToSequence[3],$_);
							}
							elsif($jobToSequence[3] == 3)
							{
								if($sequencedStatus != 1 || $sequencedStatus != 2)
								{
									my $updateClone=$dbh->prepare("UPDATE clones SET sequenced = ? WHERE name LIKE ? OR origin LIKE ?");
									$updateClone->execute($jobToSequence[3],$_,$_);
									my $updateFpcClone=$dbh->prepare("UPDATE matrix SET x = ? WHERE container LIKE 'fpcClone' AND name LIKE ?");
									$updateFpcClone->execute($jobToSequence[3],$_);
								}
							}
							elsif($jobToSequence[3] == 4)
							{
								if($sequencedStatus != 1 || $sequencedStatus != 2 || $sequencedStatus != 3)
								{
									my $updateClone=$dbh->prepare("UPDATE clones SET sequenced = ? WHERE name LIKE ? OR origin LIKE ?");
									$updateClone->execute($jobToSequence[3],$_,$_);
									my $updateFpcClone=$dbh->prepare("UPDATE matrix SET x = ? WHERE container LIKE 'fpcClone' AND name LIKE ?");
									$updateFpcClone->execute($jobToSequence[3],$_);
								}
							}
							else
							{
								if($sequencedStatus != 1 || $sequencedStatus != 2 || $sequencedStatus != 3 || $sequencedStatus != 4)
								{
									my $updateClone=$dbh->prepare("UPDATE clones SET sequenced = ? WHERE name LIKE ? OR origin LIKE ?");
									$updateClone->execute($jobToSequence[3],$_,$_);
									my $updateFpcClone=$dbh->prepare("UPDATE matrix SET x = ? WHERE container LIKE 'fpcClone' AND name LIKE ?");
									$updateFpcClone->execute($jobToSequence[3],$_);
								}
							}
							$assignedSequenceNumber++;
							print LOG "Good\n";
							last;
						}
						else
						{
							print LOG "Bad\n";
						}				
					}
				}
				unlink ("$commoncfg->{TMPDIR}/$jobToSequence[0].seq");
			}
			unlink ("$commoncfg->{TMPDIR}/$input.$$.bes");
			unlink ("$commoncfg->{TMPDIR}/$input.$$.tag");
		}
		my $updateJobBacAssigned=$dbh->do("UPDATE matrix SET barcode = $assignedSequenceNumber WHERE container LIKE 'job' AND name LIKE '$input'");
		close(LOG);	
	}

	if(@vectorId)
	{
		unlink ("$commoncfg->{TMPDIR}/$$.vector");
	}
	exit 0;
}
else{
	die "couldn't fork: $!\n";
} 

__DATA__

#example of job list

{
  "page" : 1,
  "records" : 260,
  "total" : 260,
  "rows" : [ {
    "automated" : true,
    "collectionProtocol" : "MagBead Standard Seq v2",
    "comments" : "RS_Subreads protocol",
    "copy" : false,
    "createdBy" : "dcopetti",
    "custom1" : "User Defined Field 1=",
    "custom2" : "User Defined Field 2=",
    "custom3" : "User Defined Field 3=",
    "custom4" : "User Defined Field 4=",
    "custom5" : "User Defined Field 5=",
    "custom6" : "User Defined Field 6=",
    "description" : null,
    "editable" : false,
    "expanded" : false,
    "groupName" : "16737",
    "groupNames" : [ "all" ],
    "inputCount" : 1,
    "instrumentId" : 1,
    "instrumentName" : "42219",
    "jobId" : 16737,
    "jobStatus" : "Completed",
    "lastHeartbeat" : "2014-10-07T10:49:19-0700",
    "leaf" : true,
    "modifiedBy" : null,
    "name" : "Ampl_extract_all_SubR",
    "nodeId" : 167370,
    "nodeLevel" : 0,
    "parentId" : 0,
    "path" : null,
    "plateId" : "090514_MH63_MTP_24-1~4_24-5~8_24-9~12_AmpExp_MTP_25AB-25CD-25EF-25GH_24-1~4-75_24-1~4-50_24-5~8-75_24-5~8-50_P5_0.25nM_1-10",
    "primaryProtocol" : "BasecallerV1",
    "protocolName" : "RS_Subreads.1",
    "referenceSequenceName" : null,
    "sampleName" : "AmpExp_P5_0.25nM_1-10",
    "sequencingCondition" : null,
    "showHeartbeatErr" : true,
    "version" : "2.2.0",
    "whenCreated" : "2014-10-07T10:44:15-0700",
    "whenEnded" : "2014-10-07T10:49:48-0700",
    "whenModified" : "2014-10-07T10:49:48-0700",
    "whenStarted" : "2014-10-07T10:44:16-0700"
  } ]
}

#single job
{
  "automated" : true,
  "collectionProtocol" : "MagBead Standard Seq v2",
  "comments" : "4k and 5M",
  "copy" : false,
  "createdBy" : "dcopetti",
  "custom1" : "User Defined Field 1=",
  "custom2" : "User Defined Field 2=",
  "custom3" : "User Defined Field 3=",
  "custom4" : "User Defined Field 4=",
  "custom5" : "User Defined Field 5=",
  "custom6" : "User Defined Field 6=",
  "description" : null,
  "editable" : false,
  "expanded" : false,
  "groupName" : "16743",
  "groupNames" : [ "all" ],
  "inputCount" : 1,
  "instrumentId" : 1,
  "instrumentName" : "42219",
  "jobId" : 16743,
  "jobStatus" : "Completed",
  "lastHeartbeat" : "2014-10-11T17:57:24-0700",
  "leaf" : true,
  "modifiedBy" : null,
  "name" : "ZS97_MTP1_5~8_0.2-1-10",
  "nodeId" : 167430,
  "nodeLevel" : 0,
  "parentId" : 0,
  "path" : null,
  "plateId" : "101014_MH63_MTP_41-5~8_41-9~12_42_1~4_42_5~8_ZS97_MTP_1-1~4_1_5~8_1-9~12_2-1~4_P5_0.2nM_1-10",
  "primaryProtocol" : "BasecallerV1",
  "protocolName" : "RS_HGAP_Assembly.3",
  "referenceSequenceName" : null,
  "sampleName" : "ZS97_MTP1_5~8_4k_P5_0.2nM_1-10",
  "sequencingCondition" : null,
  "showHeartbeatErr" : true,
  "version" : "2.2.0",
  "whenCreated" : "2014-10-11T14:54:08-0700",
  "whenEnded" : "2014-10-11T17:58:13-0700",
  "whenModified" : "2014-10-11T17:58:13-0700",
  "whenStarted" : "2014-10-11T14:54:11-0700"
}