use strict;
use warnings;
BEGIN{
    use Getopt::Long;
    die "Getopt::Long module is not up to date. Please fix" unless $Getopt::Long::VERSION >=2.51;
}
use Getopt::Long;
my ($gapdDir, $id, @line, $file);
GetOptions("-d=s"=> \$gapdDir);
die "Error in comand line arguments. Try again: -d dir. with the gapdfiles" unless ($gapdDir);
opendir(DIR, $gapdDir) || die "Error finding directory. Try again";
open (OUT, ">InterestingGOs.txt");
while($file=readdir DIR){
    if($file eq '.' || $file eq '..'){
        next;
    }else{
        open(IN, "${gapdDir}$file") || die "Error opening file/file non existent";
        while(<IN>){
            if($_=~/^!https:.+&geneProductId=(.+)/){
                print OUT "#$1\n";
            }elsif($_=~/^UniProtKB\t/){
                @line=split /\s+/, $_;
                print OUT "$line[3]\n";
            }
        }
        close IN;
    }
}
closedir DIR;
close OUT;
