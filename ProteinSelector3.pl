use strict;
use warnings;
BEGIN{
    use Getopt::Long;
    die "Getopt::Long is not updated. Please fix" unless $Getopt::Long::VERSION >= 2.51;
    #use Data::Dumper;
    #die "Data::Dumper is no updated. Please fix" unless $Data::Dumper::VERSION >= 2.158;
}
use Getopt::Long;
#use Data::Dumper;
my ($tsvFile, $GoFile, $out);
GetOptions("-t=s"=> \$tsvFile, "-g=s"=> \$GoFile, "-o=s"=> \$out);
die "Error in command line arguments. Try again: -t BioMart tsv file with proteins 
and -g txt file with GO terms one per line below its protein beggining with #" unless($tsvFile && $GoFile);
open (GO, $GoFile) || die "Error opening file/File non existent. Try again";
my (%GoTerms, $id, @refProt);
print "Reading and indexing GOfile...\n";
while (<GO>){
    chomp $_;
    if($_=~/^#(\w.+)/){
        $id=$1;
        push @refProt, $id;
    }elsif ($_=~/^GO:\d+/){
        $GoTerms{$_}=$id; #Asocia cada GO a su proteina
    }
}
#print Dumper \%GoTerms;
close GO;
open(IN, $tsvFile) || die "Error opening file/File non existent. Try again";
my(%Protein, %Similarity, @line, $term);
print "Reading Protein tsv file...\n";
while(<IN>){ #Selecciona las proteínas que contenga al menos uno de los términos query
    chomp $_;
    @line=split /\t/, $_;
    if ($line[2]=~/^(GO:\d+)/){# Saltate las lineas sin ningun GO
        if (exists $GoTerms{$1}){# Comprueba si es un GO de la lista proporcionada
            $term=$1;
            if(exists $Protein{$line[-1]}){# Concatena si la proteina ya ha sido seleccionada
                if($Protein{$line[-1]}=~/\Q$term\E/){#Concatena solo si el termino GO es nuevo
                    next;
                }else{             
                    $Protein{$line[-1]}= join '_', ($Protein{$line[-1]}, $term);
                    $Similarity{$line[-1]}{$GoTerms{$term}}++;#Actualiza 
                }         
            }else{
                $Protein{$line[-1]}=$term;
                foreach (@refProt){#Inicializa el hash de similitud a 0 con cada nueva entrada
                    $Similarity{$line[-1]}{$_}=0;
                }
                $Similarity{$line[-1]}{$GoTerms{$term}}++;#Guarda primer incremento
            }
        }
    }
}
#print Dumper \%Similarity;
close IN;
if($out){
    open OUT, ">${out}.txt";
}else{
    open OUT, ">SelectedProteins.txt";
}
my $i;
print "Writing results...\n";
print OUT "ProteinID\tGOs\t#GOs\t";#Escribe cabecera del fichero
for ($i=0; $i<=$#refProt; $i++){
    if($i<$#refProt){
    print OUT "#$refProt[$i]\t";
    }else{
        print OUT "#$refProt[$i]\n"
    }
}
my $GOsize;
foreach my $k (keys %Protein){
    if($Protein{$k}=~/_/){#Extrae numero de terminos GO de cada proteina encontrada
        @line=split /_/, $Protein{$k};
        $GOsize=@line; 
    }else{
        $GOsize=1;
    }
    $Protein{$k}=~s/_/|/g;
    print OUT "$k\t$Protein{$k}\t$GOsize\t";#Escribe linea de fichero
    $i=0;
    for ($i=0; $i<=$#refProt; $i++){
        if ($i<$#refProt){
            print OUT "$Similarity{$k}{$refProt[$i]}\t";
        }else{
            print OUT "$Similarity{$k}{$refProt[$i]}\n";
        }
    }
}
print "Done\n";
close OUT;