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
use userCookie;

my $userCookie = new userCookie;
my $userId = (cookie('cid')) ? $userCookie->checkCookie(cookie('cid')) : 0;
exit if (!$userId);

my $user = new user;
my $userDetail = $user->getAllFieldsWithUserId($userId);
my $userName = $userDetail->{"userName"};

my $commoncfg = readConfig("main.conf");

my $blastn = 'blast+/bin/blastn';
my $jobdir = $commoncfg->{JOBDIR};
my $tagMatchIdentity = 100; #identities for tags
my $besMatchIdentity = 98; #identities for tags
my $tagMatchPercent = 0.80; #total percent of matched tags
my $minCloneTagNumber = 5;
my $format = 'fasta';
my @jobId = param('jobId');
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
	my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});
	#connect to the mysql server
	my $relatedLibrary;
	my $cloneTagNumber;
	my $cloneBesNumber;
	my $tagTotalNumber = 0;
	my $besTotalNumber = 0;

	foreach my $input (sort {$a <=> $b} @jobId)
	{
		my $outdir = "$jobdir/$input";
		unless (-e $outdir)
		{
			mkdir $outdir;
		}
		open (ALLLOG,">>$jobdir/jobs.log") or die "can't open file: $jobdir/jobs.log";
		print ALLLOG "Assign BAC ID for Job (forcing mode): $input\n";
		close(ALLLOG);		
		open (LOG,">>$outdir/$input.log") or die "can't open file: $outdir/$input.log";
		#associate BAC and sequence
		my $updateJobToRunning=$dbh->do("UPDATE matrix SET barcode = '-1' WHERE container LIKE 'job' AND (name LIKE '$input' OR name LIKE '-$input')");
		my $jobToPool = $dbh->prepare("SELECT parent FROM link WHERE type LIKE 'poolJob' AND child = ?");
		$jobToPool->execute($input);
		while(my @jobToPool = $jobToPool->fetchrow_array())
		{
			#get pool to library
			my $pool = $dbh->prepare("SELECT * FROM matrix WHERE id = ?");
			$pool->execute($jobToPool[0]);
			my @pool = $pool->fetchrow_array();
			
			unless (exists $relatedLibrary->{$pool[4]})
			{
				my $library=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
				$library->execute($pool[4]);
				my @library = $library->fetchrow_array();
				open (TAG,">$commoncfg->{TMPDIR}/$pool[4].$$.tag") or die "can't open file: $commoncfg->{TMPDIR}/$pool[4].$$.tag";
				#get tag list
				my $libraryToTag = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'tag' AND (x = ? OR x = ?) ORDER BY o");
				$libraryToTag->execute($library[0],$library[5]);
				while (my @libraryToTag = $libraryToTag->fetchrow_array())
				{
					$cloneTagNumber->{$libraryToTag[2]}++;
					$tagTotalNumber++;
					print TAG ">$libraryToTag[2].$libraryToTag[3]\n$libraryToTag[8]\n";						
				}			
				close (TAG);
				open (BES,">$commoncfg->{TMPDIR}/$pool[4].$$.bes") or die "can't open file: $commoncfg->{TMPDIR}/$pool[4].$$.bes";
				#get tag list
				my $libraryToBes = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 98 AND (x = ? OR x = ?) ORDER BY z");
				$libraryToBes->execute($library[0],$library[5]);
				while (my @libraryToBes = $libraryToBes->fetchrow_array())
				{
					$cloneBesNumber->{$libraryToBes[2]}++;
					$besTotalNumber++;
					my $sequenceDetails = decode_json $libraryToBes[8];
					my $sequence = 'ERROR: NO SEQUENCE FOUND! PLEASE CONTACT YOUR ADMINISTRATOR.';
					my $in = Bio::SeqIO->new(-file => "$commoncfg->{DATADIR}/$sequenceDetails->{'sequence'}",
											-format => 'Fasta');
					while ( my $seq = $in->next_seq() )
					{
						$sequence = $seq->seq;
					}
					print BES ">$libraryToBes[2].$libraryToBes[6]\n$sequence\n";						
				}			
				close (BES);
				$relatedLibrary->{$pool[4]} = 1;
			}

			my $jobToSequence = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND (o = 1 OR o = 2) AND name = '' AND x = ?"); #get circularized or inserted sequences only
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
					open (CMD,"$blastn -query $commoncfg->{TMPDIR}/$jobToSequence[0].seq -subject $commoncfg->{TMPDIR}/$pool[4].$$.tag -dust no -max_target_seqs 10000 -outfmt 6 |") or die "can't open CMD: $!";
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
							my $updateSequence=$dbh->prepare("UPDATE matrix SET name = ?, z = ?, barcode = 4 WHERE id = ?"); #barcode(=4) means tagForced
							$updateSequence->execute($_,$matchedTagNumber->{$_},$jobToSequence[0]);
							my $sequencedStatus = 0;
							my $getSequenced = $dbh->prepare("SELECT * FROM clones WHERE name LIKE ? OR origin LIKE ?");
							$getSequenced->execute($_,$_);
							while (my @getSequenced = $getSequenced->fetchrow_array())
							{
								$sequencedStatus = $getSequenced[6];
							}
							if($sequencedStatus != 1 && $sequencedStatus != 2)
							{
								my $updateClone=$dbh->prepare("UPDATE clones SET sequenced = ? WHERE name LIKE ? OR origin LIKE ?");
								$updateClone->execute($jobToSequence[3],$_,$_);
								my $updateFpcClone=$dbh->prepare("UPDATE matrix SET x = ? WHERE container LIKE 'fpcClone' AND name LIKE ?");
								$updateFpcClone->execute($jobToSequence[3],$_);
							}
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
					open (CMD,"$blastn -query $commoncfg->{TMPDIR}/$jobToSequence[0].seq -subject $commoncfg->{TMPDIR}/$pool[4].$$.bes -dust no -max_target_seqs 10000 -outfmt 6 |") or die "can't open CMD: $!";
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
							if($sequencedStatus != 1 && $sequencedStatus != 2)
							{
								my $updateClone=$dbh->prepare("UPDATE clones SET sequenced = ? WHERE name LIKE ? OR origin LIKE ?");
								$updateClone->execute($jobToSequence[3],$_,$_);
								my $updateFpcClone=$dbh->prepare("UPDATE matrix SET x = ? WHERE container LIKE 'fpcClone' AND name LIKE ?");
								$updateFpcClone->execute($jobToSequence[3],$_);
							}
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
		}
		my $assignedSequenceCount = $dbh->prepare("SELECT COUNT(*) FROM matrix WHERE container LIKE 'sequence' AND name != '' AND name NOT LIKE '-%' AND o > 0 AND x = ?");
		$assignedSequenceCount->execute($input);
		my @assignedSequenceCount = $assignedSequenceCount->fetchrow_array();
		my $updateJob=$dbh->do("UPDATE matrix SET barcode = $assignedSequenceCount[0] WHERE container LIKE 'job' AND (name LIKE '$input' OR name LIKE '-$input')");
		close(LOG);	
	}
	for(keys %$relatedLibrary)
	{
		unlink ("$commoncfg->{TMPDIR}/$_.$$.tag");
		unlink ("$commoncfg->{TMPDIR}/$_.$$.bes");
	}
	exit 0;
}
else{
	die "couldn't fork: $!\n";
} 
