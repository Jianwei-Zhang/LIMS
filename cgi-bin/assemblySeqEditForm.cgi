#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use DBI;
use lib "lib/";
use lib "lib/pangu";
use pangu;
use userCookie;

my $userCookie = new userCookie;
my $userId = (cookie('cid')) ? $userCookie->checkCookie(cookie('cid')) : 0;
exit if (!$userId);

my $commoncfg = readConfig("main.conf");
my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});

my $assemblySeqId = param ('assemblySeqId') || '';
my $scrollLeft = param ('scrollLeft') || '0';

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
	0=>'NA',
	1=>'f',
	2=>'r'
	);

my %gapType = (
	0=>'None',
	1=>'5\' Gap',
	2=>'3\' Gap',
	3=>'5\' and 3\' Gaps',
	4=>'5\' Telomere',
	5=>'3\' Telomere',
	6=>'5\' and 3\' Telomeres',
	7=>'5\' Telomere and 3\' Gap',
	8=>'5\' Gap and 3\' Telomere'
	);

if ($assemblySeqId)
{
	print header;
	my $assemblySeq=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$assemblySeq->execute($assemblySeqId);
	my @assemblySeq = $assemblySeq->fetchrow_array();

	my $assembly=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$assembly->execute($assemblySeq[3]);
	my @assembly = $assembly->fetchrow_array();

	unless ($assembly[7] == 1 || $assembly[7] == 0) # exit if for frozen or running assembly
	{
		print <<END;
<script>
	closeDialog();
	errorPop("This assembly is running or frozen.");
</script>	
END
		exit;
	}
	my $target=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$target->execute($assembly[4]);
	my @target = $target->fetchrow_array();

	my $assemblySequence=$dbh->prepare("SELECT * FROM matrix WHERE id = ?");
	$assemblySequence->execute($assemblySeq[5]);
	my @assemblySequence = $assemblySequence->fetchrow_array();
	my $seqStart;
	my $seqEnd;
	if($assemblySeq[8])
	{
		($seqStart,$seqEnd) = split ",",$assemblySeq[8];
	}
	else
	{
		$seqStart = 1;
		$seqEnd = $assemblySeq[6];
		my $updateAssemblySeq=$dbh->do("UPDATE matrix SET note = '$seqStart,$seqEnd' WHERE id = $assemblySeq[0]");
	}
	my $sequenceIdSelector .= ($target[1] eq 'library') ? "<tr><td style='text-align:right'><input type='radio' id='editAssemblySeqRadio$assemblySequence[0]' name='seqId' value='$assemblySequence[0]' checked='checked' onchange='setSeqLength($assemblySequence[5])'></td><td><label for='editAssemblySeqRadio$assemblySequence[0]'>$seqType{$assemblySequence[3]} ($bacAssignType{$assemblySequence[7]}:$assemblySequence[6])</label> <a href='download.cgi?seqId=$assemblySequence[0]' title='Download this sequence' target='hiddenFrame'>" .commify($assemblySequence[5]). " bp</a> <sup><a onclick='closeDialog();openDialog(\"jobView.cgi?jobId=$assemblySequence[4]\")'>JobId: $assemblySequence[4]</a></sup></td></tr>"
		: "<tr><td style='text-align:right'><input type='radio' id='editAssemblySeqRadio$assemblySequence[0]' name='seqId' value='$assemblySequence[0]' checked='checked' onchange='setSeqLength($assemblySequence[5])'></td><td><label for='editAssemblySeqRadio$assemblySequence[0]'>$assemblySequence[2] (<a href='download.cgi?seqId=$assemblySequence[0]' title='Download this sequence' target='hiddenFrame'>" .commify($assemblySequence[5]). " bp</a>) in Genome <a onclick='closeDialog();openDialog(\"genomeView.cgi?genomeId=$assemblySequence[4]\")'>$target[2]</a></label></td></tr>";


	my $gapType;
	for (sort keys %gapType)
	{
		$gapType .= ($_ == $assemblySeq[4]) ? "<option value='$_' selected>$gapType{$_}</option>" : "<option value='$_'>$gapType{$_}</option>";
	}

	undef $/;# enable slurp mode
	my $html = <DATA>;
	$html =~ s/\$scrollLeft/$scrollLeft/g;
	$html =~ s/\$assemblySeqId/$assemblySeqId/g;
	$html =~ s/\$assemblySeqName/$assemblySeq[2]/g;
	$html =~ s/\$seqLength/$assemblySeq[6]/g;
	$html =~ s/\$gapType/$gapType/g;
	$html =~ s/\$seqStart/$seqStart/g;
	$html =~ s/\$seqEnd/$seqEnd/g;
	$html =~ s/\$sequenceIdSelector/$sequenceIdSelector/g;
	$html =~ s/\$\$/$$/g;
	print $html;
}
else
{
	print header(-type=>'text/html',-status=>'402 Invalid operation');
	exit;
}


__DATA__
	<form id="editAssemblySeq$assemblySeqId" name="editAssemblySeq$assemblySeqId" action="assemblySeqEdit.cgi" enctype="multipart/form-data" method="post" target="hiddenFrame">
	<input name="assemblySeqId" id="editAssemblySeqId$$" type="hidden" value="$assemblySeqId" />
	<input name="seqLength" id="editAssemblySeqLength$$" type="hidden" value="$seqLength" />
	<input name="scrollLeft" id="editAssemblySeqScrollLeft$$" type="hidden" value="$scrollLeft" />
	<table>
	<tr><td style='text-align:right;white-space: nowrap;'><b>Sequence</b></td>
	<td></td>
	</tr>
	$sequenceIdSelector
	<tr><td style='text-align:right;white-space: nowrap;'><b>Gap Type:</b></td>
	<td><select class='ui-widget-content ui-corner-all' name='gapType' id='gapType$$'>$gapType</select></td>
	</tr>
	<tr><td style='text-align:right;white-space: nowrap;'><b>Sequence Region:</b></td>
	<td><input type="text" name="seqStart" id="seqStart$$" value="$seqStart" onblur="changeSeqRegion();"> to <input type="text" name="seqEnd" id="seqEnd$$" value="$seqEnd" onblur="changeSeqRegion();"></td>
	</tr>
	<tr><td style='text-align:right'></td>
	<td><div id="slider-range$$"></div></td>
	</tr>
	<tr><td style='text-align:right'></td>
	<td>The selected region will be exported.</td>
	</tr>
	</table>
	</form>
<script>
$( "#slider-range$$" ).slider({
	  range: true,
	  min: 1,
	  max: $seqLength,
	  values: [ $seqStart, $seqEnd ],
	  slide: function( event, ui ) {
		$( "#seqStart$$" ).val( ui.values[ 0 ]);
		$( "#seqEnd$$" ).val( ui.values[ 1 ] );
	  }
});
$('#dialog').dialog("option", "title", "Set Sequence of $assemblySeqName");
$( "#dialog" ).dialog( "option", "buttons", [ { text: "Set", click: function() { submitForm('editAssemblySeq$assemblySeqId'); } }, { text: "Cancel", click: function() {closeDialog();} } ] );
function changeSeqRegion()
{
	var seqStart = $( "#seqStart$$" ).val();
	var seqEnd = $( "#seqEnd$$" ).val();
	var max = $( "#slider-range$$" ).slider( "option", "max" );
	var min = $( "#slider-range$$" ).slider( "option", "min" );
	if(seqStart < min)
	{
		seqStart = min;
		$( "#seqStart$$" ).val(min);
	}	
	if(seqStart > max)
	{
		seqStart = max;
		$( "#seqStart$$" ).val(max);
	}	
	if(seqEnd < min)
	{
		seqEnd = min;
		$( "#seqEnd$$" ).val(min);
	}	
	if(seqEnd > max)
	{
		seqEnd = max;
		$( "#seqEnd$$" ).val(max);
	}	
	if (seqStart > seqEnd)
	{
		seqStart = seqEnd;
		$( "#seqStart$$" ).val(seqEnd);
	}
	$( "#slider-range$$" ).slider( "values", [seqStart, seqEnd]);
}
function setSeqLength(seqLength)
{
	$( "#slider-range$$" ).slider( "option", "max", seqLength );
	$( "#slider-range$$" ).slider( "values", [1, seqLength]);
	$( "#editAssemblySeqLength$$" ).val(seqLength);
	$( "#seqStart$$" ).val(1);
	$( "#seqEnd$$" ).val(seqLength);
}
</script>