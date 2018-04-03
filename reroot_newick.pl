#!/usr/bin/perl

use strict;
use Carp qw (cluck confess croak);
use warnings;
use File::Basename;

my $dir_input = shift;
my $dir_output = shift;


if (!-e $dir_input)
{
        confess("directory does not exist");
}

##find core_output/ -name '*.nwk.out' -exec cp {} ./trees \;
my $get_files_cmd = 'find ' .  $dir_input . " -name \'*.nwk.out\'";
my @output = qx($get_files_cmd);

#print "DEBGUG $output[0]\n";
my $skipped_trees = 0;

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
	 		#print "DEBGUG $tree \n";
             
      my @musba = $tree =~ m/(Musba[\w_]+):/g;
      
      if (scalar(@musba) > 1) 
      {
          print "WARNING: Skipped tree : $filename \n $tree \n";
          $skipped_trees++;
          next;
      }
      
	 		my $outgroup = $musba[0];
         
	 		# reroot file
			my $reroot_cmd = '/homedir/rouard/bin/newick-utils-1.6/bin/nw_reroot ' . $newick_file . " $outgroup";
			#print "DEBGUG: $reroot_cmd \n";
			
		  my $rooted_tree = qx($reroot_cmd);
		  #print "DEBGUG: $rooted_tree \n";
        
		  # create file
	 		my $rooted_tree_file = $dir_output . $filename . '.rooted.nwk';
      #print "DEBGUG: $rooted_tree_file \n";
      
	 		open (OUT, ">$rooted_tree_file") or die "cannot create file : $!";
	 		print OUT $rooted_tree . "\n";
	 		close OUT;
	 }
  #close file
  close FH;
  
}


print "DEBGUG Number of skipped trees : $skipped_trees\n";

