#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use JSON; # imports encode_json, decode_json, to_json and from_json.
use DBI;
use lib "lib/";
use lib "lib/pangu";
use pangu;

my $commoncfg = readConfig("main.conf");
my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});

my $keyword = param("term");
my @results;
print header('application/json');
my $same;
my $sth=$dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'library'");
$sth->execute();
while (my @sth = $sth->fetchrow_array())
{
	my $libraryDetails = decode_json $sth[8];
	if ($libraryDetails->{'subspecies'})
	{
		if (exists $same->{$libraryDetails->{'subspecies'}})
		{
			next;
		}
		else
		{
			if ($libraryDetails->{'subspecies'} =~ /$keyword/i)
			{
				my %hash;
				$hash{'id'} = $libraryDetails->{'subspecies'};
				$hash{'label'} = $libraryDetails->{'subspecies'};
				$hash{'value'} = $libraryDetails->{'subspecies'};
				push @results, \%hash;
			}
			$same->{$libraryDetails->{'subspecies'}} = 1;
		}
	}
}
my $json = JSON->new->allow_nonref;
my $json_text = $json->encode(\@results);
print $json_text;
