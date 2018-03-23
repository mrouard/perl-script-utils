#!/usr/bin/perl

use strict;
use Carp qw (cluck confess croak);
use warnings;
use File::Basename;
use Bio::DB::Fasta;
use Bio::Seq;
use Bio::SeqIO;
use Statistics::Basic qw(:all);
use Data::Dumper;



####### input  ########
my $dir_input = shift;
my $fasta_file = shift;
#######################


if (!-e $dir_input)
{
    confess("directory does not exist");
}

print "search trees\n";

my $get_files_cmd = 'find ' .  $dir_input . " -name \'*.nwk.out\'";
my @output = qx($get_files_cmd);

if (scalar @output == 0 )
{
    confess("no newick file found");
}
#print "DEBGUG $output[0]\n"; 

my %topologies;
my @tree_type1;
my @tree_type2;


print "Process trees\n";

#load on files
foreach my $newick_file (@output)
{
    chomp $newick_file;
    
     my($filename, $dirs, $suffix) = fileparse($newick_file , '.nwk.out');
   
    #open file
    open (FH, $newick_file) or die "cannot open $newick_file : $!";

    while (<FH>)
    {
        my $tree = $_;
        $tree =~ tr/\r\n//d;
        my $tree_taxon = $tree;
        $tree_taxon =~ s/Mabur[\w_]+:/MABUR:/g;
	 	    $tree_taxon =~ s/Maban[\w_]+:/MABAN:/g;
	 	    $tree_taxon =~ s/Musba[\w_]+:/MUSBA:/g;
	 	    $tree_taxon =~ s/Mazeb[\w_]+:/MAZEB:/g;
	 	    $tree_taxon =~ s/Ma\d{2}_[gpt][\w_\.]+:/MUSAC:/g;
        while ( $tree_taxon =~ m/([0-9Ee+\-\.]+)([:,;\)])/ )
        {
            $tree_taxon =~ s/([0-9Ee+\-\.:]+)([:,;\)])/$2/;
        }
            
        #The genetic distance between banksii and burmannica in gene trees where malaccensis and banksii are sister 
        if ($tree_taxon =~ m/\(MUSAC,MABAN\)|\(MABAN,MUSAC\)/)
        {
            if ($tree_taxon =~ m/MABUR/)
            {
                #print  "DEBGUG" . $tree . "\n";
                $newick_file =~ m/(\w+_cds)/;
                push @tree_type1, $dirs.$1 . '.mafft' ;
                #print "case 1 : $newick_file \n";
            }
        }
        #The genetic distance between banksii and burmannica in gene trees where burmannica and banksii are sister
        elsif ($tree_taxon =~ m/\(MABUR,MABAN\)|\(MABAN,MABAN\)/)
        {
                #print  "2\t" . $tree . "\n";
                $newick_file =~ m/(\w+_cds)/;
                push @tree_type2, $dirs.$1 . '.mafft' ;
                #print "case 2 : $newick_file \n";
        }
        $topologies{$tree_taxon}++;
     }
    close FH;
} 
print "total number of processed trees :" . scalar @output . "\n";
print "Number of trees with topology [dAC|A,B]  ". scalar  @tree_type1 ."\ttopology [dAC|B,C] " . scalar @tree_type2  ."\n";


### calculate pairwise genetic distance with dismat emboss
print "compute genetic distances\n";

my @removed;
my @list_d1;
foreach  (@tree_type1)
{
    my $d =  compute_genetic_distance($_, 'SEQ0000002', 'SEQ0000003');
    
    if ( ($d =~ m/\d\.?\d+/) && ($d < 10) )
    {
        push @list_d1, $d;
    }
    else { push @removed, $d;}
}

my @list_d2;
foreach  (@tree_type2)
{
    my $d =  compute_genetic_distance($_, 'SEQ0000002', 'SEQ0000003');
    #print "B $d  \n";
    
    if ( ($d =~ m/\d\.?\d+/) && ($d < 10) )
    {
        push @list_d2, $d;
    }
    else { push @removed, $d;}
}

print "Number of distance values list1= ". scalar  @list_d1 ."\tlist2= " . scalar @list_d2  ."\n";
print "Number of removed values (non numerical or > 10) :" . scalar  @removed ."\n";

### calculate average and substract
print "compute D2\n";
my $avg_list1 = mean(@list_d1);
my $avg_list2 = mean(@list_d2);
my $d2 = $avg_list1 - $avg_list2;
print "\n D2 value  =  $d2  ($avg_list1 - $avg_list2) \n";

### calculate stats for pvalue
print "compute support\n";
my ($mean, $sdev) = d2_random(\@list_d1, \@list_d2, 1000);
my $z = ( $d2 - $mean)/$sdev;
print "\n mean = $mean \t stdev = $sdev \t zscore = $z \n";




###############################
####### function ##############
###############################

sub compute_genetic_distance
{
    my ($aln, $seq1, $seq2) = @_;
    my $output_file = $aln . '.dist';
    my $distmat_cmd = '/usr/local/bioinfo/EMBOSS/6.6.0/bin/distmat -sequence ' . $aln  . ' -nucmethod 4 ' . ' -outfile ' . $output_file;
    
    #does not overwrite
    if (!-e $output_file)	
	  {
	    my $aln_dist_exec  = qx($distmat_cmd);
	  }	

	  open (FH, "$output_file") or die "cannot find file : $!";
	  my $line;
	  my (@data, %lookup);
	  while ($line = <FH>)
	  {
       chomp($line);
	     if ($line =~ m/^(?:\s+(?:\d+)){3,}$/g)
	     {
          @data = (['columns', ("$line" =~ m/(\d+)/g)]);
          while (($line = <FH>) && ($line =~ m/\d\s*\t\s*\S/))
          {
            push(@data, [("\t$line" =~ m/ *\t *([^\s]*)/g)]);
            my $current_line_index = $#data;
            my $current_line_data = $data[$current_line_index];
            my $label_column_index = scalar(@$current_line_data) - 1;
            $lookup{$data[$current_line_index][$label_column_index]} = $current_line_index;
          }
	   }
	}
	close FH;
  return $data[$lookup{$seq1}][$lookup{$seq2}];
}




sub d2_random {
    my ($list1, $list2, $permut_number) = @_;
    
    my @list1 = @$list1;
    my @list2 = @$list2;
    my @d2_values;
    
    # for each permutation
    for (my $i = 1; $i <= $permut_number; $i++)
    {
        my @rand_list;
        # need to reinitialise table values
        my @merge_list = (@list1, @list2);
        
        while (@merge_list)
        {
            push @rand_list, splice(@merge_list, rand(scalar(@merge_list)), 1);
        }

        my @rand_list2 = splice (@rand_list, scalar(@list1), scalar(@list2));

        my $avg_list1 = mean(@rand_list);
        my $avg_list2 = mean(@rand_list2);
        my $d2 = $avg_list1 - $avg_list2;
        push @d2_values,  $d2;
    }
    #calculate mean and stdev
    my $d2_mean     = mean(@d2_values);
    my $d2_stdev    = stddev(@d2_values);
  
    return ($d2_mean, $d2_stdev);
}

