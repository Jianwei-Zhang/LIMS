#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use JSON::XS;
use DBI;
use lib "lib/";
use lib "lib/pangu";
use pangu;
use userCookie;

my $userCookie = new userCookie;
my $userId = (cookie('cid')) ? $userCookie->checkCookie(cookie('cid')) : 0;

print header();
if(!$userId)
{
	print <<eof;
	<h2>Error</h2>
	Please <a onclick="openDialog('login.cgi');">login</a> first!
	<script type="text/javascript">
		openDialog('login.cgi');
		loadingHide();
	</script>
eof
    exit;
}

my $commoncfg = readConfig("main.conf");
my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});

my @itemId = param ('itemId');

my $itemNumber = scalar (@itemId);
if($itemNumber < 1)
{
	print 'Select at least 1 record first!';
	exit;
}

my $readyPhrase = ($itemNumber > 1) ? "$itemNumber selected records are ready for downloading." : "$itemNumber selected record is ready for downloading.";
open (RECORD,">$commoncfg->{TMPDIR}/record.$$.txt") or die "can't open file: $commoncfg->{TMPDIR}/record.$$.txt";
my $recordLineNumber = 0;
foreach (@itemId)
{
	my $item=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$item->execute($_);
	while (my @item = $item->fetchrow_array())
	{
		my $itemDetails = decode_json $item[8];
		unless ($recordLineNumber > 0)
		{
			my $column = 0;
			for (sort {$a <=> $b} keys %$itemDetails)
			{
				$itemDetails->{$_}->{'field'} = '' unless ($itemDetails->{$_}->{'field'});
				if ($column > 0)
				{
					print RECORD "\t$itemDetails->{$_}->{'field'}";
				}
				else
				{
					print RECORD "$itemDetails->{$_}->{'field'}";
				}
				$column++;
			}
			print RECORD "\n";
		}

		my $column = 0;
		for (sort {$a <=> $b} keys %$itemDetails)
		{
			$itemDetails->{$_}->{'value'} = '' unless ($itemDetails->{$_}->{'value'});
			if ($column > 0)
			{
				print RECORD "\t$itemDetails->{$_}->{'value'}";
			}
			else
			{
				print RECORD "$itemDetails->{$_}->{'value'}";
			}
			$column++;
		}
		print RECORD "\n";
	}
	$recordLineNumber++;
}
close (RECORD);
`gzip -f $commoncfg->{TMPDIR}/record.$$.txt`;

undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$\$/$$/g;
$html =~ s/\$readyPhrase/$readyPhrase/g;
$html =~ s/\$commoncfg->{TMPURL}/$commoncfg->{TMPURL}/g;

print $html;

__DATA__
$readyPhrase
<script>
$('#dialog').dialog("option", "title", "Download Records");
$( "#dialog" ).dialog( "option", "buttons", [ { text: "Download", click: function() { location.href='$commoncfg->{TMPURL}/record.$$.txt.gz'; } }, { text: "OK", click: function() {closeDialog(); } } ] );
</script>