use warnings;
use strict;
BEGIN{
    use Getopt::Long;
    die "Getopt::Long not updated. Please fix" unless $Getopt::Long::VERSION >= 2.51;
}
use Getopt::Long;

my($genomeFile, $chrTag, $aux);
GetOptions("-g=s"=> \$genomeFile, "-c=s"=> \$chrTag);
die "Incorrect arguments. Try again: -g genome.fasta -c chromosomeTag" unless ($genomeFile && $chrTag);
open(IN, $genomeFile) || die "Error opening file/file non existent. Try again";
$aux=0;
open(OUT, ">${chrTag}.fasta");
while(<IN>){
    if($_=~/^>.?\Q$chrTag\E.+/){
        $aux=1;
        print OUT ("$_\n");
    }elsif ($aux==1){
        if($_=~/^>/){
            print OUT ("\n");
            last;
        }else{
            print OUT ("$_");
        }
    }else{
        next;
    }
}
close IN;
close OUT;