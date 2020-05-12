use strict;
use warnings;
BEGIN{
    use Getopt::Long;
    die "Getopt::Long not updated, please fix" unless $Getopt::Long::VERSION >= 2.51;
}
use Getopt::Long;
my ($file, $out);
GetOptions("-f=s"=> \$file, "-o=s"=> \$out);
die "Error in command line arguments. Try again: -f blast result txt file -o output file name" 
unless ($file);
open (OUT, ">${out}BlastMatches.txt");
print OUT "Query ID\t Blast best match ID\n";
open (IN, $file) || die "Error opening file/file non existent. Try again";
my ($aux, @line, $id);
$aux=0;
while (<IN>){
    if($_=~ /^Query\s#\d+:/ && $aux==0){
        @line=split /\s/, $_;
        $id=$line[2];
        $aux=1#Secuencia actual
    }elsif($_=~/^No\ssignficant\ssimilarity.+/ && $aux==1){
        $aux=0;
        next;#No la incluyas si no ha encontrado ninguna proteina similar
    }elsif($_=~/^Description/ && $aux==1){
        $aux=2;#Hay resultados para la secuencia
    }elsif($aux==2){
        @line=split /\s+/,$_;
        print OUT "$id\t$line[-1]\n";
        $aux=0;#Extrae exclusivamente el mejor resultado
    }
}
close IN;
close OUT; 
