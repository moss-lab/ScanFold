#!/usr/bin/perl

#This program will take a a FASTA file, divide it into windows, randomize those
#windows, and calculate Z-Scores across the sequence. This file can be used as
#the input for ScanFold-Fold (https://github.com/moss-lab/ScanFold).

# Usage:
#
# $ perl ScanFold-Scan inputfile stepsize windowsize randomizations temperature > outname
#
# Here the inputfile is your sequence in fasta format
# the stepsize and windowsize determine the extent of the length of the sequence
# fragment analyzed and how many windows to calculate a good step and window to
# start would be 40 and 120 
#
# randomizations determines how many random sequences are used in the z-score
# calc. a good start would be 10
#
# outname is the name of the outputfile you want to create. Each column will have the fist and last nt of the window, the z-score, and the minimum free energy (MFE) of the native
# a file called RNAfold.out will contain the predicted MFE structures of each window

#use threads;
#use strict;
#use threads::shared;
use warnings;
use Bio::SeqIO;



my $fastafile = $ARGV[0];
my $StepSize = $ARGV[1];
my $WindowSize = $ARGV[2];
my $Randomizations = $ARGV[3];
my $Temperature = $ARGV[4];

#Opens the file or quits and prints and error if the file cannot open
open(FASTAFILE, "$fastafile") || die "Can't open fasta file\n";
my @fasta = <FASTAFILE>;
close FASTAFILE;

my $seqio = Bio::SeqIO->new(-file => $fastafile, '-format' => 'Fasta');
while(my $seq = $seqio->next_seq) {
  #my $string = $seq->seq;
  # do stuff with $string


my $Name = $seq->display_id;
#chomp $Name;
#print $Name;
my $SEQ = $seq->seq();

#Divide the sequence into windows scramble them and calc Z-Zcore
my $Length = length $SEQ;


print "i\tj\tTemperature\tNative_dG\tZ-score\tP-value\tED\tfMFE\tSequence\tMFE_fold\tCentroid\t#A's\t#G's\t#C's\t#U's\t$Name\n";

for (my $i = 0; $i <= ($Length - $WindowSize); $i += $StepSize) {
    my $Frag = substr($SEQ, $i, $WindowSize);

#Uncomment these two lines to take calculate for reverse strand
    #Take Reverse Complement
    # my $Frag = ReverseComplement ($Frag);

		my $StartNt = $i + 1;

        my @Out = `echo $Frag | RNAfold -T $Temperature --noPS -p`;
        #open (OUT, ">> $fastafile.out") || die "can't open $fastafile.out\n";
        #print OUT ">$StartNt \n";
        #print OUT @Out;
        #close OUT;

        my $Fold_seq = $Out[0];
        chomp $Fold_seq;
        my $FoldData = $Out[1];
          my @FoldData = split(/\s+/, $FoldData);
          my $Fold = $FoldData[0];
          my $CentroidData = $Out[3];
          my @Centroid = split(/\s+/, $CentroidData);
          my $CentroidFold = $Centroid[0];

        #chomp $Fold;
        my $Data = $Out[4];
        my @Data = split(/\s+/, $Data);
        my $FreqMFE = $Data[7];
        #chomp $FreqMFE;
        my @FreqMFE_SF = split(/\;/, $FreqMFE);
        my $FreqMFE_SF = $FreqMFE_SF[0];
        #my $FreqMFE_SF = $FreqMFE_SSF[0];
        #chomp $FreqMFE_SF;
        my $EnsembleDiversity = $Data[10];

        #Put the native sequence in the first position always!
        my @seqarray = ();
        push (@seqarray, $Frag);
        my @ScrambledSeqs = Scramble ($Frag);
        push (@seqarray, @ScrambledSeqs);

    my $NTFreqs = &NucFreqs($Frag);
	my @EnergyArray = &Energy(@seqarray);
	my $Zscore = &ZScore(@EnergyArray);
    my $NativeDG = $EnergyArray[0];
	chomp $NativeDG;
	my $PValue = &PValue(@EnergyArray);

	$WinStart = $i + 1;
	$WinEnd = $i + $WindowSize;

   	print "$WinStart\t$WinEnd\t$Temperature\t$NativeDG\t$Zscore\t$PValue\t$EnsembleDiversity\t$FreqMFE_SF\t$Frag\t$Fold\t$CentroidFold\t$NTFreqs\n";


}

print "\n\n";


######Sub-routine to scramble RNAs################################
sub Scramble {

        my $InSeq = $_[0];
        my @Out = ();

        for (my $i = 0; $i < $Randomizations; $i++) {
                my $OutSeq = "";
                my @InSeq = split ("", $InSeq);

                while (@InSeq > 0) {
                        my $Rand = rand(@InSeq);
                        my $RandBase = splice(@InSeq, $Rand, 1);
                        $OutSeq .= $RandBase;
                        }
                push(@Out, $OutSeq);
        }
        return @Out;
}

######Sub-routine to calculate MFEs using RNAfold#################
sub Energy {

    my @engarr = @_;
    my $k = 0;
    my @returnarray = ();

    foreach my $Sequence (@engarr) {
	#print "SEQ: $Sequence\n";
    #open (TEMP, "> Temp.seq") || die "can't open tempfile\n";
    #my $input = "$Sequence";
    #print TEMP "$input";
    #close TEMP;

    my @Out = `echo $Sequence | RNAfold -T $Temperature --noPS`;
    #print @Out;
    #print "\n";
    my $STR_EN =  $Out[1];
    my @STR_EN = split(/\s+\(/, $STR_EN);
    my $EN = $STR_EN[1];
    #print "$EN\n";
    $EN =~ s/\(//g;
    $EN =~ s/\)//g;

    push (@returnarray, $EN);

}
return @returnarray;
@returnarray = ();
}

######Sub-routine to calculate Nucleotide frequencies#######
sub NucFreqs {

        my $InSeq = $_[0];
        $InSeq =~ s/T/U/g;

        $Gs = 0;
        $Cs = 0;
        $As = 0;
        $Us = 0;

        while ($InSeq =~ /G/g) {$Gs += 1;}
        while ($InSeq =~ /C/g) {$Cs += 1;}
        while ($InSeq =~ /A/g) {$As += 1;}
        while ($InSeq =~ /U/g) {$Us += 1;}

	my $length = $As + $Cs + $Gs + $Us;
        my $OutFreqs = "$As\t$Gs\t$Cs\t$Us";
        return  $OutFreqs;
}

######Sub-routine to calculate Z-scores########
sub ZScore {

	my @arr = @_;
	#print "This is my arr @arr\n";
	my $sum = 0;
	my $count = @arr;
	foreach my $l(@arr){


		$sum += $l;
                #print "DG: $l\n";

}
	my $Average = $sum/$count;
        my $average = ($sum-$arr[0])/($count-1);
        #print "sum: $sum\n";
  #print "average = $average\n";
  #print "Average = $Average\n";
  $sum = "";
  my $MFE = $arr[0];
	#print "AVG: $Average\n";
  #print "MFE = $MFE\n";
	my $Sigma = 0;
	my $Sum = 0;

	foreach my $m(@arr){
    #print "m is $m\n";
		#print "My i is $i\n";
		$Sigma = ($m - $average)**2;
		#print "My Sigma is $Sigma\n";
		$Sum += $Sigma;
}
	my $SD = sqrt($Sum/$count);
        #print "SUM: $Sum\n";
        #print "SD: $SD\n";
	$Sum = "";

  my $return = (($SD ne 0) ? (($arr[0]-$Average)/$SD):"Undefined");
        my $ZScore = "$return\t$average";
        chomp $ZScore;
				my @ZScore_Data = split(/\s+/, $ZScore);
        #$ZScore =~ s/\\n//g;
        #print "My return is $return\n";

		#my $Output = substr($ZScore, 0, 5);
		my $Output = substr($ZScore_Data[0], 0, 5);
	return $Output;
}

######Sub-routine to calculate P-value (fraction of scrambled dG < native dG########
sub PValue {

	my @arr = @_;
	my $BelowNative = 0;
	my $TotalCount = @arr;

	my $Native = $arr[0];
	foreach my $l (@arr) {
	    if ($l < $Native) {$BelowNative += 1;}
	}

	my $Fraction = ($BelowNative / $TotalCount);
	my $Output = substr($Fraction, 0, 4);

	return $Output;

}

##########Sub-routine to generate reverse complement of window sequence###########

sub ReverseComplement {

    my $Input = $_[0];

    # reverse the sequence
    my $RevComp = reverse ($Input);

    # complement the reversed sequence
    $RevComp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
    return $RevComp;

}
}
