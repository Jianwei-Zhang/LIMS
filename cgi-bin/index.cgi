#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser ); 
use URI::Escape::XS;
#use Date::Calc qw(:all);
use JSON::XS;
use DBI;
use Bio::SeqIO;
use lib "lib/";
use lib "lib/pangu";
use pangu;
use config;

my $config = new config;
my $copyright = $config->getFieldValueWithFieldName("COPYRIGHT");
$copyright = escapeHTML($copyright);
$copyright =~ s/\n/<br>/g;
my $siteName = $config->getFieldValueWithFieldName("SITENAME");
$siteName = escapeHTML($siteName);
my $slogan = $config->getFieldValueWithFieldName("SLOGAN");
$slogan = escapeHTML($slogan);
my $keywords = $config->getFieldValueWithFieldName("KEYWORDS");
$keywords = escapeHTML($keywords);

# my $year1 = "2011";
# my $month1 = "8";
# my $day1 = "1";
# my ($year2,$month2,$day2) = Today();
# my $Dd = Delta_Days($year1,$month1,$day1,$year2,$month2,$day2);
my $active = cookie('menu') || 0;;
my $commoncfg = readConfig("main.conf");
undef $/;# enable slurp mode
my $html = <DATA>;
$html =~ s/\$active/$active/;
$html =~ s/\$commoncfg->{HTDOCS}/$commoncfg->{HTDOCS}/g;
$html =~ s/\$copyright/$copyright/g;
$html =~ s/\$siteName/$siteName/g;
$html =~ s/\$slogan/$slogan/g;
$html =~ s/\$keywords/$keywords/g;

## Data files will be saved to $commoncfg->{TMPDIR}
unless (-e "$commoncfg->{DATADIR}")
{
	mkdir "$commoncfg->{DATADIR}";
}
## sequences will be saved to $commoncfg->{DATADIR}/sequences
unless (-e "$commoncfg->{DATADIR}/sequences")
{
	mkdir "$commoncfg->{DATADIR}/sequences";
}
## alignments will be saved to $commoncfg->{DATADIR}/alignments
unless (-e "$commoncfg->{DATADIR}/alignments")
{
	mkdir "$commoncfg->{DATADIR}/alignments";
}
unless (-e "$commoncfg->{DATADIR}/alignments/seqToSeq")
{
	mkdir "$commoncfg->{DATADIR}/alignments/seqToSeq";
}
unless (-e "$commoncfg->{DATADIR}/alignments/seqToSet")
{
	mkdir "$commoncfg->{DATADIR}/alignments/seqToSet";
}
unless (-e "$commoncfg->{DATADIR}/alignments/setToSet")
{
	mkdir "$commoncfg->{DATADIR}/alignments/setToSet";
}
## Temporary files will be saved to $commoncfg->{TMPDIR}
unless (-e "$commoncfg->{TMPDIR}")
{
	mkdir "$commoncfg->{TMPDIR}";
}

my $update = 0;

my $pid = fork();
if ($pid)
{
	print header;
	print $html;
}
elsif($pid == 0)
{
 	close (STDOUT);
 	exit unless ($update);
	my $dbh=DBI->connect("DBI:mysql:$commoncfg->{DATABASE}:$commoncfg->{DBHOST}",$commoncfg->{USERNAME},$commoncfg->{PASSWORD});
# 	#update old assembly data
# 	my $assembly = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assembly'");
# 	$assembly->execute();
# 	while(my @assembly = $assembly->fetchrow_array())
# 	{
# 		my $updateToAssemblySeq = $dbh->do("UPDATE matrix SET container = 'assemblySeq' WHERE container LIKE 'assemblyClone' AND o = $assembly[0]");
# 		my $assemblyCtgs = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblyCtg' AND o = ?");
# 		$assemblyCtgs->execute($assembly[0]);
# 		while(my @assemblyCtgs = $assemblyCtgs->fetchrow_array())
# 		{
# 			for (split ",", $assemblyCtgs[8])
# 			{
# 				next unless ($_);
# 				$_ =~ s/[^a-zA-Z0-9]//g;
# 				my $seqName = $_;
# 				my $assemblySeq = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblySeq' AND o = ? AND name LIKE ?");
# 				$assemblySeq->execute($assembly[0],$seqName);
# 				while (my @assemblySeq = $assemblySeq->fetchrow_array())
# 				{
# 					$assemblyCtgs[8] =~ s/$seqName/($assemblySeq[0])/g;
# 				}
# 			}
# 			my $updateToAssemblyCtg = $dbh->prepare("UPDATE matrix SET note = ? WHERE id = ?");
# 			$updateToAssemblyCtg->execute($assemblyCtgs[8],$assemblyCtgs[0]);
# 		}
# 	}

# 	#update old sequence data
# 	my $getSequences = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence'");
# 	$getSequences->execute();
# 	while (my @getSequences =  $getSequences->fetchrow_array())
# 	{
# 		if ($getSequences[8] =~ /^>/)
# 		{
# 			my $seqFile = "$commoncfg->{TMPDIR}/$getSequences[0].seq";
# 			open (FILE, ">$seqFile");
# 			print FILE $getSequences[8];
# 			close FILE;
# 			my $in = Bio::SeqIO->new(-file => $seqFile,
# 				-format => 'Fasta');
# 			while ( my $seq = $in->next_seq() )
# 			{
# 				my $seqDetails;
# 				$seqDetails->{'id'} = $seq->id;
# 				$seqDetails->{'description'} = $seq->desc || '';
# 				$seqDetails->{'sequence'} = $seq->seq;
# 				$seqDetails->{'gapList'} = '';
# 				my $seqend=0;
# 				foreach (split(/([N|n]{20,})/,$seq->seq)) #at least 20 Ns to be a gap
# 				{
# 					my $seqstart=$seqend+1;
# 					$seqend=$seqend + length($_);
# 					if($_ =~ /^[N|n]+$/)
# 					{
# 						$seqDetails->{'gapList'} .= ($seqDetails->{'gapList'} ne '') ? ",$seqstart-$seqend" : "$seqstart-$seqend" ;
# 					}
# 				}
# 				my $json = JSON::XS->new->allow_nonref;
# 				$seqDetails = $json->encode($seqDetails);
# 				my $updateToSequence = $dbh->prepare("UPDATE matrix SET note = ? WHERE id = ?");
# 				$updateToSequence->execute($seqDetails,$getSequences[0]);
# 			}
# 			unlink($seqFile);			
# 		}
# 	}

	#update old Assembly data
# 	my $getAssembly = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assembly'");
# 	$getAssembly->execute();
# 	while (my @getAssembly =  $getAssembly->fetchrow_array())
# 	{
# 		if ($getAssembly[8] !~ /^{/)
# 		{
# 			my $assemblyDetails;
# 			$assemblyDetails->{'description'} = $getAssembly[8];
# 			my $json = JSON::XS->new->allow_nonref;
# 			$assemblyDetails = $json->encode($assemblyDetails);
# 			my $updateToAssembly = $dbh->prepare("UPDATE matrix SET note = ? WHERE id = ?");
# 			$updateToAssembly->execute($assemblyDetails,$getAssembly[0]);
# 		}
# 	}
# 	my $sequenceNameToId;
# 	my $genomeId = '3039246';
# 	my $getSequences = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'sequence' AND o = 99 AND x = ?");
# 	$getSequences->execute($genomeId);
# 	while (my @getSequences =  $getSequences->fetchrow_array())
# 	{
# 		$sequenceNameToId->{$getSequences[2]} = $getSequences[0];
# 	}
# 	my $assemblyId = '3037999';
# 	my $getAssemblySequences = $dbh->prepare("SELECT * FROM matrix WHERE container LIKE 'assemblySeq' AND o = ?");
# 	$getAssemblySequences->execute($assemblyId);
# 	while (my @getAssemblySequences =  $getAssemblySequences->fetchrow_array())
# 	{
# 		my $updateAssemblySeq = $dbh->prepare("UPDATE matrix SET y = ? WHERE id = ?");
# 		$updateAssemblySeq->execute($sequenceNameToId->{$getAssemblySequences[2]},$getAssemblySequences[0]);
# 	}
	exit 0;
}
else
{
	die "couldn't fork: $!\n";
} 


__DATA__
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	 "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" lang="en-US" xml:lang="en-US">
<head>
<title>$siteName - $slogan</title>
<link rev="made" href="mailto:zhangjw@email.arizona.edu" />
<meta name="keywords" content="$keywords" />
<meta name="copyright" content="coypright $copyright" />
<link type="image/png" rel="icon" href="$commoncfg->{HTDOCS}/favicon.ico" />
<link type="image/x-icon" rel="shortcut icon" href="$commoncfg->{HTDOCS}/favicon.ico" />
<link rel="stylesheet" type="text/css" href="$commoncfg->{HTDOCS}/css/smoothness/jquery-ui-1.10.4.custom.min.css" />
<link rel="stylesheet" type="text/css" href="$commoncfg->{HTDOCS}/css/dataTables.jqueryui.css" />
<link rel="stylesheet" type="text/css" href="$commoncfg->{HTDOCS}/css/lims.css" />
<script src="$commoncfg->{HTDOCS}/js/jquery-1.10.2.js" type="text/javascript"></script>
<script src="$commoncfg->{HTDOCS}/js/jquery-ui-1.10.4.min.js" type="text/javascript"></script>
<script src="$commoncfg->{HTDOCS}/js/d3.v3.min.js" type="text/javascript"></script>
<script src="$commoncfg->{HTDOCS}/js/newick.js" type="text/javascript"></script>
<script src="$commoncfg->{HTDOCS}/js/d3.phylogram.js" type="text/javascript"></script>
<script src="$commoncfg->{HTDOCS}/js/jquery.dataTables.min.js" type="text/javascript"></script>
<script src="$commoncfg->{HTDOCS}/js/dataTables.jqueryui.js" type="text/javascript"></script>
<script src="$commoncfg->{HTDOCS}/js/jquery.ui-contextmenu.min.js" type="text/javascript"></script>
</head>
<body id="top">
<table width="100%" border="0" cellspacing="0" cellpadding="0">
	<tr>
		<td width="15%"><a href="$commoncfg->{HTDOCS}"><img src="pangulogo.cgi"/></a></td>
		<td style="vertical-align: bottom;align:left;"><div id="slogan" style="float:right;padding: .7em;"><i>$slogan</i></div><h2>$siteName</h2></td>
		<td width="200" style="vertical-align: top;align:left;" class="header">
		<div id="login" class="ui-widget-content ui-corner-all" style="padding: .5em .5em; margin-top: .1em; margin-bottom: .1em;">
			</div>
		</td>
	</tr>
	<tr>
		<td colspan="3">
			<div id="menu">
				<ul>
					<li><a href="home.cgi">Home</a></li>
					<li><a href="resource.cgi">Resource</a></li>
					<li><a href="storage.cgi">Storage</a></li>
					<li><a href="assemblyHome.cgi">Assembly</a></li>
					<li><a href="general.cgi">General</a></li>
					<li><a href="setting.cgi">Settings</a></li>
					<li><a href="feedback.cgi">Contact</a></li>
				</ul>
			</div>
		</td>
	</tr>
	<tr align="center">
		<td colspan="3">
	    &copy;$copyright <a href="http://www.genome.arizona.edu" target="_blank">Arizona Genomics Institute</a>.
		</td>
	</tr>
</table>
<div id="back-top">
	<a href="#top"><span></span>Back to Top</a>
</div>
<div id='loadingPage' class='ui-state-error ui-corner-all' onclick='loadingHide()'><img src='$commoncfg->{HTDOCS}/css/images/loading.gif'> Page is loading...</div>
<div id='savingPage' class='ui-state-highlight ui-corner-all'><img src='$commoncfg->{HTDOCS}/css/images/loading.gif'> Saving...</div>
<div id="information" class="ui-state-highlight ui-corner-all" style="display:none;padding:1em;z-index:10000;"></div>
<div id="error" class="ui-state-error ui-corner-all" style="display:none;padding:1em;z-index:10000;" onclick="$(this).hide();"></div>
<div id='dialog' style="display:none"></div>
<div id='viewer' style="display:none"></div>
<div id='hiddenDiv' style="display:none">If you can see this line, your browser is not working properly.</div>

<iframe name='hiddenFrame' id="hiddenFrame" style='display:none'>If you can see this line, your browser is not working properly.</iframe>

<script type="text/javascript">
	var myTimeout;
	$.ajaxSetup({cache: false});
	$('#login').fadeOut("slow").load('status.cgi').fadeIn("slow");
	$("#menu").tabs({
		// loading spinner
		beforeLoad: function(event, ui) {
			ui.panel.html('<img src="$commoncfg->{HTDOCS}/css/images/loading.gif" style="vertical-align:middle;"> Loading...');
		},
		active: $active
	});
	$("#menu .ui-tabs-panel").css('padding','0px');
	$( ".toolTip" ).tooltip();
	$( "#dialog" ).dialog().dialog("close");
	$( "#viewer" ).dialog().dialog("close");
	$('#hiddenDiv').hide();
	loadingHide();
	savingHide();

	// hide #back-top first
	$("#back-top").hide();
	// fade in #back-top
	$(window).scroll(function () {
		if ($(this).scrollTop() > 100) {
			$('#back-top').fadeIn();
		} else {
			$('#back-top').fadeOut();
		}
	});
	// scroll body to 0px on click
	$('#back-top a').click(function () {
		$('body,html').animate({
			scrollTop: 0
		}, 'slow');
		return false;
	});
	function refresh(tabId) {
		var current_index = $("#" + tabId).tabs("option","active");
		$("#" + tabId).tabs('load',current_index);
	}
	function checkAll(checkBoxName) {
		$("input[name="+ checkBoxName +"]").prop('checked', true);
		return false;
	}
	function uncheckAll(checkBoxName) {
		$("input[name="+ checkBoxName +"]").prop('checked', false);
		return false;
	}
	function checkOrUncheck(checkBoxId){
		if ($('#'+checkBoxId).prop('checked')){
			$('#'+checkBoxId).prop('checked', false);
		}
		else{
			$('#'+checkBoxId).prop('checked', true);
		}
	}
	function checkClass(className) {
		$("input."+ className).prop('checked', true);
		return false;
	}
	function uncheckClass(className) {
		$("input."+ className).prop('checked', false);
		return false;
	}
	function loadingShow()
	{
		$('#loadingPage').show();
		return false;
	}
	function loadingHide()
	{
		$('#loadingPage').hide();
		return false;
	}
	function savingShow()
	{
		$('#savingPage').show();
		return false;
	}
	function savingHide()
	{
		$('#savingPage').hide();
		return false;
	}

	function informationPop(message)
	{
		$('#information')
		.html(message)
		.css({
			position:'fixed',
			left: ($(window).width() - $('#information').outerWidth())/2,
			top: ($(window).height() - $('#information').outerHeight())/5
			})
		.fadeIn("fast")
		.delay(5000)
		.fadeOut("slow");
	}

	function errorPop(message)
	{
		$('#error')
		.html(message)
		.css({
			position:'fixed',
			left: ($(window).width() - $('#error').outerWidth())/2,
			top: ($(window).height() - $('#error').outerHeight())/5
			})
		.fadeIn("fast")
		.delay(5000)
		.fadeOut("slow");
	}

	function closePops() {
		$('#error').fadeOut();
		$('#information').fadeOut();
	}

	function loaddiv(divname,address)
	{
		$('#'+divname).empty().append("<img src='$commoncfg->{HTDOCS}/css/images/loading.gif'>Loading...").load(address);
		return false;
	}
	function submitForm(formId)
	{
		$( "#" + formId ).submit(); 
		return false;
	}
	function openDialog(dialogAddress)
	{
		$('#dialog').dialog({
			modal: true,
			open: function ()
			{
				$(this).empty();
		        $(this).append("<img src='$commoncfg->{HTDOCS}/css/images/loading.gif'>Loading...");
				$(this).load(dialogAddress);
			},
			close: function ()
			{
				clearTimeout(myTimeout);
			},
			width: '650',
			position: { my: "center", at: "center", of: window },
			title: "^_^",
			buttons:[],
			draggable: true,
			resizable: false
		});
		return false;
	}
	function openViewer(viewerAddress)
	{
		$('#viewer').dialog({
			open: function ()
			{
				$(this).empty();
		        $(this).append("<img src='$commoncfg->{HTDOCS}/css/images/loading.gif'>Loading...");
				$(this).load(viewerAddress);
			},
			close: function ()
			{
				clearTimeout(myTimeout);
			},
			width: '600',
			position: { my: "center", at: "center", of: window },
			title: "^_^",
			buttons:[],
			draggable: true,
			resizable: true
		});
		return false;
	}
	function openDialogForm(formURL,formId)
	{
		var str = $( "#" + formId ).serialize();
		$.post(formURL,str,
			function(data){
				$('#dialog').dialog({
				modal: true,
				open: function ()
				{
					$(this).empty();
					$(this).append("<img src='$commoncfg->{HTDOCS}/css/images/loading.gif'>Loading...");
					$(this).html(data);
				},
				close: function ()
				{
					clearTimeout(myTimeout);
				},
				width: '650',
				position: { my: "center", at: "center", of: window },
				title: "^_^",
				buttons:[],
				draggable: true,
				resizable: false
			});
		});
		return false;
	}

	function closeDialog()
	{
		$('#dialog').dialog("close");
		return false;
	}
	function closeViewer()
	{
		$('#viewer').dialog("close");
		return false;
	}
	function highlight(id){
		$("#"+id).addClass("ui-state-hover");
	}
	function unhighlight(id){
		$("#"+id).removeClass("ui-state-hover");
	}
	function editIconShow(id)
	{
		$("#"+id).append("<span style='position: absolute;left: 0px;top: -10px;display:inline-block;' class='ui-icon ui-icon-pencil'></span>");
	}
	function editIconHide(id)
	{
		$("#"+id).children("span").remove();
	}

	function mouseoverIcon(iconName)
	{
		$('#'+iconName).addClass('ui-state-hover').css('cursor', 'pointer');
	}

	function mouseoutIcon(iconName)
	{
		$('#'+iconName).removeClass('ui-state-hover');
	}

	function deleteItem(itemId,option)
	{
		var r=confirm("Are you sure to delete this?");
		if (r==true)
		{
			loaddiv("hiddenDiv" , "itemDelete.cgi?items=" + itemId + "&option=" + option);
		}
		return false;
	}
	function revokeAssignment(seqId)
	{
		var r=confirm("Are you sure to revoke this assigned BAC name?");
		if (r==true)
		{
			loaddiv("hiddenDiv" , "jobRevokeBac.cgi?seqId=" + seqId );
		}
		return false;
	}

	function printDiv(divId) {
		var HTMLdata= $("#"+divId).html();
		var mywindow = window.open('', 'new div', 'height=400,width=600');
		mywindow.document.write('<html><head><title>Printing...</title>');
		mywindow.document.write('<link rel="stylesheet" type="text/css" href="$commoncfg->{HTDOCS}/css/smoothness/jquery-ui-1.10.4.custom.min.css" /><link rel="stylesheet" type="text/css" href="$commoncfg->{HTDOCS}/css/dataTables.jqueryui.css" /><link rel="stylesheet" type="text/css" href="$commoncfg->{HTDOCS}/css/lims.css" />');
		mywindow.document.write('</head><body >');
		mywindow.document.write(HTMLdata);
		mywindow.document.write('</body></html>');
		mywindow.print();
		//mywindow.close();
	    //return true;
	}

	function buttonInit()
	{
		$("button").button();
		$("input:submit").button();
		$("input:reset").button();
		$("input:button").button();
	}

	function datepickerInit()
	{
		$( ".datepicker" ).datepicker({
			changeMonth: true,
			changeYear: true,
			yearRange: "1900:-18"
		});
	}

	function wordCount(noun){
		$('.word_count').each(function() {
		var input = '#' + this.id;
		var count = input + '_count';
		$(count).show();
		word_count(input, count, noun);
		$(this).keyup(function() { word_count(input, count, noun) });
		});
	}
	function word_count(field, count, noun) {
		var number = 0;
		var matches = $(field).val().match(/\b/g);
		if(matches) {
			number = matches.length/2;
		}
		$(count).text( number + ' ' + noun + (number > 1 ? 's' : ''));
	}

	//check browser if IE or not
	if( navigator.userAgent.toLowerCase().match(/msie/))
	{
		alert( 'Your IE browser is poorly supported by this site! Some pages might not be fully functional.' );
	}
</script>
</body>
</html>
