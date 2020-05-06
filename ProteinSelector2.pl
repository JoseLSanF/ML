use strict;
use warnings;
BEGIN{
    use Getopt::Long;
    die "Getopt::Long is not updated. Please fix" unless $Getopt::Long::VERSION >= 2.51;
}
use Getopt::Long;

my ($tabFile, @GoTerms, $gffFile);
GetOptions("-t=s"=> \$tabFile, "-g=s@"=> \@GoTerms, "-f=s"=> \$gffFile);
die "Error in command line arguments. Try again: -t TabFile with proteins -f GffFile with genes and
-g GO:x1 -g GO:x2 -g GO:x3...etc" unless($tabFile && @GoTerms && $gffFile);
my %GoTerms=map{$_=>1} @GoTerms; #Convierte array en hash para buscar más rápido
open(IN, $tabFile) || die "Error opening file/File non existent. Try again";
my(%Protein, @line, $size, @GoTermsTargets, $term);
while(<IN>){ #Selecciona las proteínas que contenga al menos uno de los términos query
    chomp $_;
    @line=split /\t/, $_;
    $size=@line;
    if ($size>=14){#Comprueba que acaso haya términos GO, que haya 14 columnas
       if($line[13]=~/\d+|GO/){#Varios GO
           @GoTermsTargets=split /\|/, $line[13];
           foreach $term (@GoTermsTargets){
               if(exists $GoTerms{$term}){#Comprueba que el termino sea uno de los buscados
                   if(exists $Protein{$line[0]}){#Comprueba que la proteina ya haya sido incluida
                        $Protein{$line[0]}=join('|', ($Protein{$line[0]}, $term));
                   }else{
                        $Protein{$line[0]}=$term;
                   }
               }
           }
       }else{
           if(exists $GoTerms{$line[13]}){
               if(exists $Protein{$line[0]}){
                   $Protein{$line[0]}=join('|', ($Protein{$line[0]}, $line[13]));
               }else{
                   $Protein{$line[0]}=$line[13];
               } 
           }
       }
    }else{
        next;
    }
}
close IN;
open(GFF, $gffFile) || die "Error opening file/File non existent. Try again";
my (%GenPosition, $aux, $id);
$aux=0;
while(<GFF>){#Extrae las posiciones de los genes de las proteinas seleccionadas
    if($_=~/^#\sstart\sgene\s(.+)/ && $aux==0){
        if(exists $Protein{$1}){
            $aux=1;
            $id=$1;
        }
    }elsif($_=~/gene/ && $aux==1){
            @line=split /\t/, $_;
            $GenPosition{$id}=join '-', ($line[3], $line[4]);
            $aux=0;
    }
}
close GFF;
open (OUT, ">InterestingProteins2.txt");
foreach my $k (keys %Protein){
    print OUT "$k\t$Protein{$k}\t$GenPosition{$k}\n";
}
close OUT;

