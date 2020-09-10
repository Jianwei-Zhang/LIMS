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
my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});

my $chrNumber = param ('chrNumber');
my @assemblyCtgId = param('assemblyCtgId');
my $position;
print header;

for my $assemblyCtgId(@assemblyCtgId)
{
	my $assemblyCtg=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$assemblyCtg->execute($assemblyCtgId);
	my @assemblyCtg = $assemblyCtg->fetchrow_array();
	my $assembly=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$assembly->execute($assemblyCtg[3]);
	my @assembly = $assembly->fetchrow_array();
	unless ($assembly[7] == 1 || $assembly[7] == 0) # exit if for frozen or running assembly
	{
		print <<END;
<script>
	parent.closeDialog();
	parent.errorPop("This assembly is running or frozen.");
</script>	
END
		exit;
	}

	if ($assembly[5])
	{
		if ($chrNumber == 0)
		{
			my $update=$dbh->do("UPDATE matrix SET x = 0, z = 0 WHERE id = $assemblyCtgId ");
		}
		else
		{
			my $genomeSeq=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 99 AND x = $assembly[5] AND z = ? ");
			$genomeSeq->execute($chrNumber);
			my @genomeSeq = $genomeSeq->fetchrow_array();
			my $maxAlign=0;
			$position=0;
			foreach (split ",", $assemblyCtg[8])
			{
				$_ =~ s/[^a-zA-Z0-9]//g;
				my $assemblySeq=$dbh->prepare ("SELECT * FROM matrix WHERE id = ?");
				$assemblySeq->execute($_);
				while (my @assemblySeq = $assemblySeq->fetchrow_array())
				{
					my $queryDir;
					for (my $position = 0; $position < length($assemblySeq[5]); $position += 2)
					{
						$queryDir .= "/q". substr($assemblySeq[5],$position,2);
					}
					my $subjectDir;
					for (my $position = 0; $position < length($genomeSeq[0]); $position += 2)
					{
						$subjectDir .= "/s". substr($genomeSeq[0],$position,2);
					}
					open (TBL,"$commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$assemblySeq[5]-$genomeSeq[0].tbl") or die "can't open file: $commoncfg->{DATADIR}/alignments/seqToSeq$queryDir$subjectDir/$assemblySeq[5]-$genomeSeq[0].tbl";
					while(<TBL>)
					{
						/^#/ and next;
						my @getAlignment = split("\t",$_);
						if ($assemblySeq[7] > 0)
						{
							my $pos1=$getAlignment[8] < $getAlignment[9] ? $getAlignment[8] : $getAlignment[9];
							my $pos2=$getAlignment[6] < $getAlignment[7] ? $getAlignment[6] : $getAlignment[7];
							$position=$pos1-$pos2 if ($maxAlign < $getAlignment[3]);
							$maxAlign=$maxAlign > $getAlignment[3] ? $maxAlign :$getAlignment[3];
						}
						else
						{
							my $pos1=$getAlignment[8] < $getAlignment[9] ? $getAlignment[8] : $getAlignment[9];
							my $pos2=$getAlignment[6] > $getAlignment[7] ? $getAlignment[6] : $getAlignment[7];
							$position=$pos1-$assemblySeq[6]+$pos2 if ($maxAlign < $getAlignment[3]);
							$maxAlign=$maxAlign > $getAlignment[3] ? $maxAlign :$getAlignment[3];
						}
					}
					close(TBL);
				}
			}
			my $update=$dbh->do("UPDATE matrix SET x = $chrNumber,z = $position WHERE id = $assemblyCtgId ");
		}
	}
	else
	{
		my $update=$dbh->do("UPDATE matrix SET x = $chrNumber WHERE id = $assemblyCtgId ");
	}
}
	print <<END;
<script>
	parent.closeDialog();
	parent.informationPop("Chromosome number assignment is done.");
	parent.refresh("menu");
</script>	
END



