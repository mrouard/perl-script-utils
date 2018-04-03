#!/usr/bin/perl

=pod

=head1 NAME

run_pipeline.pl - Execute GreenPhyl pipeline

=head1 SYNOPSIS

    run_pipeline.pl

=head1 REQUIRES

Perl5

=head1 DESCRIPTION

This script runs the GreenPhyl phylogeny pipeline. It takes as input a
multiple sequence FASTA file and computes a phylogeny analysis.

=cut

use strict;
use Carp qw (cluck confess croak);
use warnings;
use Getopt::Long;
use Pod::Usage;
use Error qw(:try);
use Fatal qw (open close);
use Cwd qw(abs_path getcwd);
use threads;
use threads::shared;

use Bio::SeqIO;
use Bio::AlignIO;
use Bio::TreeIO;
use Bio::LocatableSeq;
use Bio::SimpleAlign;

use GreenPhylConfig;



# Script global constants
##########################

=pod

=head1 CONSTANTS

B<$GREENPHYL_NAME>: (string)

name (displayed) of the pipeline.

B<$GREENPHYL_GROUP_ID>: (integer)

Unix group ID of GreenPhyl group. The valud can be obtained using the command:
% perl -e 'print [getgrnam('greenphyl')]->[2] . "\n";'

B<$DEFAULT_ADMIN_EMAIL>: (string)

e-mail address used as sender for notification e-mails.

B<$FIRST_COLUMN_WIDTH>: (integer)

number of characters in the first column in screen reports.

B<$INDENTATION>: (integer)

number of spaces to use for indentation.

B<$INPUT_FILE_EXTENSION>: (string)

file extension used for input files.

B<$OUTPUT_FILE_EXTENSION>: (string)

file extension used by output files.

B<$DICTIONARY_FILE_EXTENSION>: (string)

file extension used by dictionaries.

B<$MAFFT_FILE_EXTENSION>: (string)

file extension for MAFFT alignment file.

B<$HMM_FILE_EXTENSION>: (string)

file extension for HMM file.

B<$CLUSTAL_FILE_EXTENSION>: (string)

file extension for CLUSTAL files.

B<$MASKING_FILE_EXTENSION>: (string)

file extension for masked alignment files.

B<$FILTERING_FILE_EXTENSION>: (string)

file extension for filtered alignment files.

B<$BSP_MASKING_FILE_EXTENSION>: (string)

file extension for bsp masked alignment files.

B<$PHYLIP_FILE_EXTENSION>: (string)

file extension for phylip alignment files.

B<$NEWICK_FILE_EXTENSION>: (string)

file extension for newick tree files.

B<$NEWICK_EXTENDED_FILE_EXTENSION>: (string)

file extension for extended newick tree files.

B<$PHYLIP_SCRIPT_FILENAME>: (string)

script file name used to issue command to PHYLIP programs.

B<$PHYLIP_OUTTREE_FILENAME>: (string)

name of the output file for tree (PHYLIP package).

B<$PHYML_TREE_FILE_SUFFIX>: (string)

suffix of the file name used by PhyML for phylogeny tree.
Note: name is imposed by PhyML program.

B<$PHYML_BOOTSTRAP_FILE_SUFFIX>: (string)

suffix of the file name used by PhyML for bootstrap trees.
Note: name is imposed by PhyML program.

B<$PHYML_DISTANCE_MATRIX_FILE_SUFFIX>: (string)

suffix of the file name used by PhyML for distance matrix.
Note: name is imposed by PhyML program.

B<$PHYLOGENY_FILE_SUFFIX>: (string)

name of phylogeny tree.

B<$SDI_OUTPUT_FILE_EXT>: (string)

SDI output file name.

B<$BOOTSTRAP_TREE_FILENAME_SEED>: (string)

seed for bootstrap filenames.

B<$RIO_OUTPUT_FILE_SEED>: (string)

RIO output file name seed used to generate other names.

B<$RIO_DATA_OUTPUT_FILE_EXTENSION>: (string)

RIO output data file extension.

B<$RIO_TREE_OUTPUT_FILE_EXTENSION>: (string)

RIO output tree file extension.

B<$STEP_ALIGNMENT>: (string)

alignment step short name.

B<$STEP_HMM>: (string)

HMM step short name.

B<$STEP_MASKING>: (string)

masking step short name.

B<$STEP_PHYLOGENY>: (string)

phylogeny step short name.

B<$STEP_TREE_ROOTING>: (string)

tree rooting step short name.

B<$STEP_RESAMPLE_INFERENCE_OF_ORTHOLOGY>: (string)

inference of orthologs step short name.

=cut

our $GREENPHYL_NAME                       = 'GreenPhyl';
our $GREENPHYL_GROUP_ID                   = 33462;
our $DEFAULT_ADMIN_EMAIL                  = 'guignon@cirad.fr';
our $FIRST_COLUMN_WIDTH                   = 30;
our $INDENTATION                          = 4;
our $INPUT_FILE_EXTENSION                 = '.src';
our $OUTPUT_FILE_EXTENSION                = '.out';
our $DICTIONARY_FILE_EXTENSION            = '.dic';
our $FASTA_FILE_EXTENSION                 = '.fasta';
our $MAFFT_FILE_EXTENSION                 = '.mafft';
our $HMM_FILE_EXTENSION                   = '.hmm';
our $CLUSTAL_FILE_EXTENSION               = '.clw';
our $MASKING_FILE_EXTENSION               = '.mask';
our $FILTERING_FILE_EXTENSION             = '.fltr';
our $BSP_MASKING_FILE_EXTENSION           = '.bsp';
our $MASKING_HTML_FILE_EXTENSION          = "_masking.html";
our $FILTERED_HTML_FILE_EXTENSION         = "_filtered.html";
our $FILTERED_SEQ_FILE_EXTENSION          = "_filtered.seq";
our $PHYLIP_FILE_EXTENSION                = '.phy';
our $NEWICK_FILE_EXTENSION                = '.nwk';
our $NEWICK_EXTENDED_FILE_EXTENSION       = '.nhx';
our $PHYLOXML_FILE_EXTENSION              = '.xml';
our $PHYLIP_SCRIPT_FILENAME               = 'phylip.script';
our $PHYLIP_OUTTREE_FILENAME              = 'outtree';
our $PHYML_TREE_FILE_SUFFIX               = '_phyml_tree'; # do not change: see doc
our $PHYML_BOOTSTRAP_FILE_SUFFIX          = '_phyml_boot_trees'; # do not change: see doc
our $PHYML_DISTANCE_MATRIX_FILE_SUFFIX    = '_phyml_dist.txt'; # do not change: see doc
our $PHYLOGENY_FILE_SUFFIX                = '_phylogeny_tree' . $NEWICK_FILE_EXTENSION;
our $ROOTED_TREE_FILE_SUFFIX              = '_rap_gene_tree' . $NEWICK_FILE_EXTENSION;
our $RECONCILIED_TREE_FILE_SUFFIX         = '_rap_reconcilied_gene_tree' . $NEWICK_FILE_EXTENSION;
our $STATS_TREE_FILE_SUFFIX               = '_rap_stats_tree.txt';
our $XML_TREE_FILE_SUFFIX                 = '_rap_tree.xml';

our $SDI_OUTPUT_FILE                      = 'sdir_outfile.xml';
our $SDI_OUTPUT_FILE_EXT                  = '.sdi.xml';
our $BOOTSTRAP_TREE_FILENAME_SEED         = 'bootstrap_';
our $RIO_OUTPUT_FILE_SEED                 = 'rio_';
our $RIO_DATA_OUTPUT_FILE_EXTENSION       = '.txt';
our $RIO_TREE_OUTPUT_FILE_EXTENSION       = '.xml';

our $STEP_ALIGNMENT                       = 'alignment';
our $STEP_HMM                             = 'hmm';
our $STEP_MASKING                         = 'masking';
our $STEP_PHYLOGENY                       = 'phylogeny';
our $STEP_TREE_ROOTING                    = 'rooting';
our $STEP_RESAMPLE_INFERENCE_OF_ORTHOLOGY = 'orthology';

our $SPECIES_CODE_REGEXP                  = '_[A-Z0-9]{3,5}';




# Script global variables
##########################

=pod

=head1 VARIABLES

B<$g_debug>: (boolean)

debug mode flag. Set to 1 when "-debug" command line parameter is used.
Default: 0 (false)

B<$g_output_dir>: (string)

path of the output directory where computation should be operated.

B<$g_sequences_count>: (integer)

number of sequence in the source FASTA file.

B<$g_main_filename>: (string)

main file name without any extension.

B<$g_resume>: (string)

if set, resume an analyse from the specified step. Allowed steps are:
'alignment', 'masking', 'phylogeny', 'rooting', 'rio'.

B<$g_autoresume>: (string)

if set, only start a step if its outputs are missing. Default: no auto-resume.

B<$g_end>: (string)

if set, the anlyse will stop after the specified step.

=cut

our $g_debug = 1;

# working directory path
our $g_output_dir = $GreenPhylConfig::OUTPUT_PATH;

our $g_sequences_count = 0;

our $g_main_filename = '';

our $g_resume = $STEP_ALIGNMENT;

our $g_autoresume = 0;

our $g_end = $STEP_RESAMPLE_INFERENCE_OF_ORTHOLOGY;

our $gs_thread_messages :shared;

our $g_stdout = \*STDOUT;
our $g_stderr = \*STDERR;





# Script global functions
##########################

=pod

=head1 FUNCTIONS

=head2 printWarning

B<Description>: display a warning message on STDERR.

B<ArgsCount>: 1

=over 4

=item $message: (string) (R)

warning message to display.

=back

B<Example>:

    printWarning("the file $filename has been replaced!\nYou have been warned!");

=cut

sub printWarning
{
    my ($message) = @_;
    # check arguments
    if (1 != @_)
    {
        confess "usage: printWarning(message);";
    }

    # remove trailing invisible characters
    $message =~ s/\s+$//s;
    print {$g_stderr} "\n";
    # check for multi-lines warning
    if ($message =~ m/[\r\n]/)
    {
        # multi-lines
        my @message_lines = split(/(?:\n|\r\n|\r)/, $message);
        print {$g_stderr} "WARNING: " . shift(@message_lines) . "\n         ";
        print {$g_stderr} join("\n         ", @message_lines);
        print {$g_stderr} "\n";
    }
    else
    {
        # single line
        print {$g_stderr} "WARNING: " . $message . "\n";
    }
    return;
}




=pod

=head2 printError

B<Description>: display an error message on STDERR.

B<ArgsCount>: 1

=over 4

=item $message: (string) (R)

error message to display.

=back

B<Example>:

    printError("Computation failed!");

=cut

sub printError
{
    my ($message) = @_;
    # check arguments
    if (1 != @_)
    {
        confess "usage: printError(message);";
    }

    # remove trailing invisible characters
    $message =~ s/\s+$//s;
    print {$g_stderr} "\n";
    # check for multi-lines warning
    if ($message =~ m/[\r\n]/)
    {
        # multi-lines
        my @message_lines = split(/(?:\n|\r\n|\r)/, $message);
        print {$g_stderr} "ERROR: " . shift(@message_lines) . "\n         ";
        print {$g_stderr} join("\n       ", @message_lines);
        print {$g_stderr} "\n\n";
    }
    else
    {
        # single line
        print {$g_stderr} "ERROR: " . $message . "\n\n";
    }
    return;
}




=pod

=head2 printDebug

B<Description>: display a debug message on STDERR if debug mode is ON.

B<ArgsCount>: 1

=over 4

=item $message: (string) (R)

debug message to display.

=back

B<Example>:

    printDebug("Starting the process using paramaters:!\na=25\nb=50\nc=3");

=cut

sub printDebug
{
    my ($message) = @_;
    # check arguments
    if (1 != @_)
    {
        confess "usage: printDebug(message);";
    }

    # check if debugging is disabled
    if (!$g_debug)
    {
        #debugging disabled, nothing to do!
        return;
    }

    # remove trailing invisible characters
    $message =~ s/\s+$//s;
    print {$g_stderr} "\n";
    # check for multi-lines warning
    if ($message =~ m/[\r\n]/)
    {
        # multi-lines
        my @message_lines = split(/(?:\n|\r\n|\r)/, $message);
        print {$g_stderr} "DEBUG: " . shift(@message_lines) . "\n       ";
        print {$g_stderr} join("\n       ", @message_lines);
        print {$g_stderr} "\n";
    }
    else
    {
        # single line
        print {$g_stderr} "DEBUG: " . $message . "\n";
    }
    return;
}




=pod

=head2 printStageStart

B<Description>: display the name of a stage when started.

B<ArgsCount>: 1

=over 4

=item $stage_name: (string) (R)

Name of the stage without any EOL character.

=back

B<Example>:

    printStageStart('Alignment');

=cut

sub printStageStart
{
    my ($stage_name) = @_;
    # check arguments
    if (1 != @_)
    {
        confess "usage: printStageStart(stage_name);";
    }

    # compute spaces to add
    my $reminding_spaces_count = $FIRST_COLUMN_WIDTH - $INDENTATION - length($stage_name);
    if (0 > $reminding_spaces_count)
    {
        $reminding_spaces_count = 0;
    }
    print {$g_stdout} ' ' x $INDENTATION . $stage_name . '.' x $reminding_spaces_count;
    # in case of debug mode, add a line-feed
    if ($g_debug)
    {
        print {$g_stderr} "\n"; #
    }
}




=pod

=head2 printStageEnd

B<Description>: display the status of a stage when terminating.

B<ArgsCount>: 1

=over 4

=item $status: (string) (R)

Status of the stage execution.

=back

B<Example>:

    printStageEnd('OK');

=cut

sub printStageEnd
{
    my ($status) = @_;
    # check arguments
    if (1 != @_)
    {
        confess "usage: printStageEnd(status);";
    }

    # indent if debug mode
    if ($g_debug)
    {
        print {$g_stdout} ' ' x $FIRST_COLUMN_WIDTH;
    }

    # remove trailing invisible characters
    $status =~ s/\s+$//s;
    print {$g_stdout} $status . "\n";
}




=pod

=head2 executeMAFFT

B<Description>: execute MAFFT alignment program on the given FASTA input
file and returns the alignment.

B<ArgsCount>: 1

=over 4

=item $input_filename: (string) (R)

the path to the input FASTA file (not aligned).

=back

B<Return>: (string)

the alignment filename.

B<Exception>:

B<Example>:

    my $alignment_filename = executeMAFFT("sequences.fasta");

=cut

sub executeMAFFT
{
    my ($input_filename) = @_;
    # parameters check
    if (1 != @_)
    {
        confess "usage: executeMAFFT(input_fasta_name);";
    }
    # check if a filename was provided
    if (not $input_filename)
    {
        confess "Input file name is missing!";
    }
    # make sure the infile exists and is readable
    if (not -r $input_filename)
    {
        confess "Unable to read file '$input_filename'!";
    }

    # display stage start
    printStageStart("Alignment (MAFFT)");

    # set up the output file name
    my $alignment_filename = $g_main_filename . $MAFFT_FILE_EXTENSION;

    # make sure the output file does not exist
    if (-e $alignment_filename)
    {
        printWarning("MAFFT: existing output file \"$alignment_filename\" will be replaced!");
        unlink $alignment_filename;
    }

    # prepare environment
    $ENV{'MAFFT_BINARIES'} = $GreenPhylConfig::MAFFT_BINARY_PATH;
    # adjust command line parameters
    my $cmd = $GreenPhylConfig::MAFFT_CMD;
    if ((2 < $g_sequences_count) && ($g_sequences_count <= 200))
    {
        printDebug("MAFFT parameters einsi maxiterate 1000");
        $cmd = "$GreenPhylConfig::MAFFT_CMD --ep 0 --maxiterate 1000 --genafpair $input_filename >$alignment_filename 2>$g_main_filename$GreenPhylConfig::MAFFT_LOG_FILENAME";
    }
    elsif ((200 < $g_sequences_count) && ($g_sequences_count <= 500))
    {
        printDebug("MAFFT parameters einsi maxiterate 1");
        $cmd = "$GreenPhylConfig::MAFFT_CMD --maxiterate 1 $input_filename >$alignment_filename 2>$g_main_filename$GreenPhylConfig::MAFFT_LOG_FILENAME";
    }
    elsif ((500 < $g_sequences_count) && ($g_sequences_count <= 10000))
    {
        printDebug("MAFFT FFTNS");
        $cmd = "$GreenPhylConfig::MAFFT_FFTNS_CMD --ep 0 $input_filename >$alignment_filename 2>$g_main_filename$GreenPhylConfig::MAFFT_LOG_FILENAME";
    }

    # remove previous error logs
    unlink "$g_output_dir/mafft_err.log";

    # run MAFFT
    printDebug("COMMAND: $cmd");
    if (0 != system($cmd))
    {
        if (not $!)
        {
            confess "MAFFT: see logs ($g_main_filename$GreenPhylConfig::MAFFT_LOG_FILENAME) for errors.";
        }
        confess "MAFFT: $!";
    }

    # check if we got an output alignment
    if (-s $alignment_filename)
    {
        printStageEnd("OK");
    }
    else
    {
        printStageEnd("Failed!");
        if (-e $alignment_filename)
        {
            confess "Failed to compute alignment: alignment file '$alignment_filename' was empty!";
        }
        else
        {
            confess "Failed to compute alignment: alignment file '$alignment_filename' was missing!";
        }
    }
    return $alignment_filename;
}




=pod

=head2 getMAFFTOutput

B<Description>: returns file names generated by MAFFT alignment program.

B<ArgsCount>: 1

=over 4

=item $input_filename: (string) (R)

the path to the input FASTA file (not aligned).

=back

B<Return>: (string)

the alignment filename.

B<Exception>:

B<Example>:

    my $alignment_filename = getMAFFTOutput("sequences.fasta");

=cut

sub getMAFFTOutput
{
    my ($input_filename) = @_;
    # parameters check
    if (1 != @_)
    {
        confess "usage: getMAFFTOutput(input_fasta_name);";
    }
    # check if a filename was provided
    if (not $input_filename)
    {
        confess "Input file name is missing!";
    }
    return $g_main_filename . $MAFFT_FILE_EXTENSION;
}


=pod

=head2 executeHMMBUILD

B<Description>: execute HMMBUILD program on MAFTT output
file and returns the Hmm.

B<ArgsCount>: 2

=over 4

=item $input_filename: (string) (R)

the path to the input MAFFT file (aligned).

=back

B<Return>: (string)

the hmm filename.

B<Exception>:

B<Example>:

    my $hmm_filename = executeHMMBUILD("output.mafft");

=cut

sub executeHMMBUILD
{
    my ($input_filename) = @_;
    # parameters check
    if (1 != @_)
    {
        confess "usage: executeHMMBUILD(output.mafft);";
    }
    # check if a filename was provided
    if (not $input_filename)
    {
        confess "Input file name is missing!";
    }
    # make sure the infile exists and is readable
    if (not -r $input_filename )
    {
        confess "Unable to read file '$input_filename'!";
    }

    # display stage start
    #printStageStart("HMM (HMMBUILD)");

    # set up the output file name
    my $hmm_filename = getHMMBUILDOutput($input_filename);

    # make sure the output file does not exist
    if (-e $hmm_filename)
    {
        printWarning("HMMBUILD: existing output file \"$hmm_filename\" will be replaced!");
        unlink $hmm_filename;
    }

    # adjust command line parameters
    # my $cmd = "$GreenPhylConfig::HMMBUILD_CMD --informat afa $hmm_filename $input_filename 1>$g_main_filename$GreenPhylConfig::HMMBUILD_LOG_FILENAME 2>>$g_main_filename$GreenPhylConfig::HMMBUILD_LOG_FILENAME";
    my $cmd = "$GreenPhylConfig::HMMBUILD_CMD $hmm_filename $input_filename 1>$g_main_filename$GreenPhylConfig::HMMBUILD_LOG_FILENAME 2>>$g_main_filename$GreenPhylConfig::HMMBUILD_LOG_FILENAME";

    # run HMMBUILD
    printDebug("COMMAND: $cmd");
    if (0 != system($cmd))
    {
        if (not $!)
        {
            confess "HMMBUILD: see logs ($g_main_filename$GreenPhylConfig::HMMBUILD_LOG_FILENAME) for errors.";
        }
        confess "HMMBUILD: $!";
    }

    # check if we got an hmm output
    printStageStart("HMM (HMMBUILD)");
    if (-s $hmm_filename)
    {
        printStageEnd("OK");
    }
    else
    {
        printStageEnd("Failed!");
        if (-e $hmm_filename)
        {
            confess "Failed to compute hmm: hmm file '$hmm_filename' was empty!";
        }
        else
        {
            confess "Failed to compute hmm: hmm file '$hmm_filename' was missing!";
        }
    }
    return $hmm_filename;
}

=pod

=head2 getHMMBUILDOutput

B<Description>: returns file names generated by HMMBUILD program.

B<ArgsCount>: 1

=over 4

=item $input_filename: (string) (R)

the path to the alignment file (mafft output).

=back

B<Return>: (string)

the hmm filename.

B<Exception>:

B<Example>:

    my $hmm_filename = getHMMBUILDOutput("alignment.mafft");

=cut

sub getHMMBUILDOutput
{
    my ($input_filename) = @_;
    # parameters check
    if (1 != @_)
    {
        confess "usage: getHMMBUILDOutput(input_alignment_name);";
    }
    # check if a filename was provided
    if (not $input_filename)
    {
        confess "Input file name is missing!";
    }
    return $g_main_filename . $HMM_FILE_EXTENSION;
}

=pod

=head2 executeRefiner

B<Description>: [function description]. #+++

B<ArgsCount>: [count of arguments] #+++

=over 4

=item variable_name: ([variable nature]) ([requirements]) #+++ see below

#--- requirement can be:
#--- (R)=required,
#--- (O)=optional
#--- (U)=optional and must be undef if omitted
[variable description]. #+++

=item variable_name2: ([variable nature]) ([requirements]) #+++

[variable description]. #+++

=back

B<Return>: ([return type]) #+++

[return description]. #+++

B<Exception>:

=over 4

[=item * exception_type:

description, case when it occurs.] #+++ see below

=back

#--- Example:
#---
#---=over 4
#---
#---=item * Range error:
#---
#---thrown when "argument x" is outside the supported range, for instance "-1".
#---
#---=back

B<Example>:

    my $refined_alignment_filename = executeRefiner("alignment.fasta");

=cut

sub executeRefiner
{
    my ($input_filename) = @_;
    # parameters check
    if (1 != @_)
    {
        confess "usage: executeRefiner(input_alignment_name);";
    }
    # check if a filename was provided
    if (not $input_filename)
    {
        confess "Input file name is missing!";
    }
    # make sure the infile exists and is readable
    if (not -r $input_filename )
    {
        confess "Unable to read file '$input_filename'!";
    }

    printStageStart("Refiner (GBlocks)");

# my $gblocks_outfile = $input_filename . "-gb";
# $cmd = join(" ",$GreenPhylConfig::GBLOCKS_CMD, $outfile, "-t=p ", "1>$g_output_dir/mask4.log", "2>$g_output_dir/mask4_err.log"); #+val
#+val: less stringent:  $cmd = join(" ",$GreenPhylConfig::GBLOCKS_CMD, $outfile, "-t=p -b4=5 -b5=h", "1>$g_output_dir/mask4.log", "2>$g_output_dir/mask4_err.log"); #+val
#+val: even less stringent:  $cmd = join(" ",$GreenPhylConfig::GBLOCKS_CMD, $outfile, "-t=p -b4=2 -b5=a", "1>$g_output_dir/mask4.log", "2>$g_output_dir/mask4_err.log"); #+val
#+val: WARNING: sometimes, Gblocks does not find any relevant blocks and returns an empty alignment!
#+val: WARNING: sometimes, Gblocks selects a very small part of the original alignment and such a refined alignment may not be relevant!
#+val: therefore, the amount of selected position should be checked (above a minimum % of positions) before running into the next step and the refined alignment my have to be discarded
# my $outfile = $gblocks_outfile;
    my $output_filename = $input_filename;
    printStageEnd("Not implemented!");

    return $output_filename;
}




=pod

=head2 executeMaskingAndFiltering

B<Description>: remove some columns of the alignment.

B<ArgsCount>: 1

=over 4

=item $input_filename: (string) (R)

input alignment file name (FASTA format).

=back

B<Return>: (list of string)

the name of the masked alignment file, the FASTA version of that file and the
name of the filtered sequence file.
...

B<Example>:

    my ($output_filename, $filtered_fasta_alignment_filename, $filtered_seq_filename, $output_mask_filename1, $output_mask_filename2, $output_mask_filename3, $bsp_filename1, $bsp_filename2, $masking_html_filename1, $masking_html_filename2, $filtered_html_filename) = executeMaskingAndFiltering("alignment.fasta");

=cut

sub executeMaskingAndFiltering
{
    my ($input_filename) = @_;
    # parameters check
    if (1 != @_)
    {
        confess "usage: executeMaskingAndFiltering(input_alignment_name);";
    }
    # check if a filename was provided
    if (not $input_filename)
    {
        confess "Input file name is missing!";
    }
    # make sure the infile exists and is readable
    if (not -r $input_filename )
    {
        confess "Unable to read file '$input_filename'!";
    }

    printStageStart("Mask. & Filt. (Trimal)");

    # update what is the input file and what will be the output file
    my ($output_filename,
        $filtered_fasta_alignment_filename,
        $filtered_seq_filename,
        $output_mask_filename1,
        $output_mask_filename2,
        $output_mask_filename3,
        $bsp_filename1,
        $bsp_filename2,
        $masking_html_filename1,
        $masking_html_filename2,
        $filtered_html_filename,
       ) = getMaskingAndFilteringOutputs($input_filename);
    my $filtered_sequence_ratio = 0;

    if ((-e $output_mask_filename1)
        || (-e $output_mask_filename2)
        || (-e $output_mask_filename3)
        || (-e $filtered_fasta_alignment_filename)
        || (-e $output_filename)
        || (-e $bsp_filename1)
        || (-e $bsp_filename2)
        || (-e $masking_html_filename1)
        || (-e $masking_html_filename2)
        || (-e $filtered_html_filename)
        || (-e $filtered_seq_filename))
    {
        printWarning("Trimal: existing Trimal output files will be replaced!");
        unlink(
            $output_mask_filename1,
            $output_mask_filename2,
            $output_mask_filename3,
            $filtered_fasta_alignment_filename,
            $output_filename,
            $bsp_filename1,
            $bsp_filename2,
            $masking_html_filename1,
            $masking_html_filename2,
            $filtered_html_filename,
            $filtered_seq_filename,
        );
    }


    # Do the masking
    # -first step: remove gaps
    my $cmd = "$GreenPhylConfig::TRIMAL_CMD -in $input_filename -out $output_mask_filename1 -colnumbering -htmlout \"$masking_html_filename1\" -gt 0.5 >$bsp_filename1 2>$g_main_filename$GreenPhylConfig::TRIMAL_LOG_FILENAME";
    printDebug("COMMAND: $cmd");
    if (0 != system($cmd))
    {
        if (not $!)
        {
            confess "Trimal: see logs ($g_main_filename$GreenPhylConfig::TRIMAL_LOG_FILENAME) for errors.";
        }
        confess "Trimal 1: $!";
    }

    # -second step: do the masking
    #$cmd = "$GreenPhylConfig::TRIMAL_CMD -in $output_mask_filename1 -out $output_mask_filename2 -colnumbering -htmlout \"$masking_html_filename2\" -cons 50 -st 0.001 -w 1 >$bsp_filename2 2>>$g_main_filename$GreenPhylConfig::TRIMAL_LOG_FILENAME";
    $cmd = "$GreenPhylConfig::TRIMAL_CMD -in $output_mask_filename1 -out $output_mask_filename2 -colnumbering -htmlout \"$masking_html_filename2\" -cons 50 -w 1 >$bsp_filename2 2>>$g_main_filename$GreenPhylConfig::TRIMAL_LOG_FILENAME";
    printDebug("COMMAND: $cmd");
    if (0 != system($cmd))
    {
        if (not $!)
        {
            confess "Trimal: see logs ($g_main_filename$GreenPhylConfig::TRIMAL_LOG_FILENAME) for errors.";
        }
        confess "Trimal 2: $!";
    }
    
    # check if we got an output alignment
    if (!-s $output_mask_filename2)
    {
        printStageEnd("Failed!");
        if (-e $output_mask_filename2)
        {
            confess "Failed to mask alignment: masked alignment file was empty!";
        }
        else
        {
            confess "Failed to mask alignment: masked alignment file was missing!";
        }
    }

    # Do the sequence filtering
    $cmd = "$GreenPhylConfig::TRIMAL_CMD -in $output_mask_filename2 -out $output_mask_filename3 -resoverlap 0.1 -seqoverlap 10 -htmlout \"$filtered_html_filename\" >>$g_main_filename$GreenPhylConfig::TRIMAL_LOG_FILENAME 2>>$g_main_filename$GreenPhylConfig::TRIMAL_LOG_FILENAME";
    printDebug("COMMAND: $cmd");

    if (0 != system($cmd))
    {
        if (not $!)
        {
            confess "Trimal: see logs ($g_main_filename$GreenPhylConfig::TRIMAL_LOG_FILENAME) for errors.";
        }
        confess "Trimal 3: $!";
    }

    # check if we got an output alignment
    if (-s $output_mask_filename3)
    {
        printStageEnd("OK");
    }
    else
    {
        printStageEnd("Failed!");
        if (-e $output_mask_filename3)
        {
            confess "Failed to filter alignment: filtered alignment file was empty!";
        }
        else
        {
            confess "Failed to filter alignment: filtered alignment file was missing!";
        }
    }

    # list filtered sequences
    # open source alignment file
    my $input_fasta_align = Bio::AlignIO->new(
        '-file'   => $input_filename,
        '-format' => 'fasta',
        '-alphabet' => 'protein',
    );
    my %sequences_set; # keys are sequence name, values are juste true values (1)
    # list sequences in a hash
    my $alignment = $input_fasta_align->next_aln();
    foreach my $bio_sequence ($alignment->each_seq)
    {
        $sequences_set{$bio_sequence->display_id} = 1;
    }
    my $initial_sequence_count = scalar(keys(%sequences_set));

    # open filtered alignment file
    my $output_phylip_align = Bio::AlignIO->new(
        '-file'   => $output_mask_filename3,
        '-format' => 'fasta',
        '-alphabet' => 'protein',
    );

    # list sequence removed (delete hash keys of kept sequences that must
    # contain a minimal number of residue corresponding to 50% of the alignment
    # length)
    if ($alignment = $output_phylip_align->next_aln())
    {
        my $half_alignment_length = 0;
        foreach my $bio_sequence ($alignment->each_seq)
        {
            # get alignment length
            $half_alignment_length ||= $bio_sequence->length()/2;
            if (exists($sequences_set{$bio_sequence->display_id}))
            {
                # only keep sequences that contain more than half of the alignment length in residues
                if ((1 < $bio_sequence->length())
                    && ($half_alignment_length > ($bio_sequence->seq() =~ tr/-//)))
                {
                    # remove key from hash
                    delete($sequences_set{$bio_sequence->display_id});
                }
                else
                {
                    # if we get here, it's because trimal did not its job: sequence should have already been removed
                    $alignment->remove_seq($bio_sequence);
                    # print "DEBUG: sequence " . $bio_sequence->display_id . " does not contain enough residue to be kept\n";
                }
            }
            else
            {
                cluck "Warning: sequence name '" . $bio_sequence->display_id . "' in filtered alignment ($output_mask_filename3) not found in original alignment ($input_filename)!\n";
            }
        }
    }
    else
    {
        # all sequence removed!
        cluck "Warning: all sequences of original alignment ($input_filename) have been filtered!\n";
    }


    # check if more than 10% of the sequences have been removed
    my $removed_sequence_count = scalar(keys(%sequences_set));
    $filtered_sequence_ratio = $removed_sequence_count / $initial_sequence_count;
    if ($filtered_sequence_ratio > 0.1)
    {
        cluck sprintf("more than a 10%% of the sequences (%.1f%%) have been dropped!\n", 100*$removed_sequence_count/$initial_sequence_count);
    }

    # create PHYLIP version of the alignment
    $output_phylip_align = Bio::AlignIO->new(
        '-file'   => ">$output_filename",
        '-format' => 'phylip',
        '-alphabet' => 'protein',
    );
    $output_phylip_align->write_aln($alignment);

    # create FASTA version of the alignment
    my $output_fasta_align = Bio::AlignIO->new(
        '-file'   => ">$filtered_fasta_alignment_filename",
        '-format' => 'fasta',
        '-alphabet' => 'protein',
    );
    $output_fasta_align->write_aln($alignment);


    # save list a file
    my $filtered_seq_fh;
    if (!open($filtered_seq_fh, ">$filtered_seq_filename"))
    {
        confess "Failed to save list of filtered sequences in '$filtered_seq_filename'!\n";
    }
    print {$filtered_seq_fh} join("\n", keys(%sequences_set));
    close($filtered_seq_fh);
    
    return ($output_filename, $filtered_fasta_alignment_filename, $filtered_seq_filename, $output_mask_filename1, $output_mask_filename2, $output_mask_filename3, $bsp_filename1, $bsp_filename2, $masking_html_filename1, $masking_html_filename2, $filtered_html_filename, $filtered_sequence_ratio);
}




=pod

=head2 getMaskingAndFilteringOutputs

B<Description>: returns the file names generated by the masking step.

B<ArgsCount>: 1

=over 4

=item $input_filename: (string) (R)

input alignment file name (FASTA format).

=back

B<Return>: (list of string)

the name of the masked alingnment file, the FASTA version of that file and
the name of the filtered sequence file.
...

B<Example>:

    my ($output_filename, $filtered_fasta_alignment_filename, $filtered_seq_filename, $output_mask_filename1, $output_mask_filename2, $output_mask_filename3, $bsp_filename1, $bsp_filename2, $masking_html_filename1, $masking_html_filename2, $filtered_html_filename) = getMaskingAndFilteringOutputs("alignment.fasta");

=cut

sub getMaskingAndFilteringOutputs
{
    my ($input_filename) = @_;
    # parameters check
    if (1 != @_)
    {
        confess "usage: getMaskingAndFilteringOutputs(input_alignment_name);";
    }
    # check if a filename was provided
    if (not $input_filename)
    {
        confess "Input file name is missing!";
    }

    my $output_filename                   = $g_main_filename . $PHYLIP_FILE_EXTENSION;
    my $filtered_fasta_alignment_filename = $g_main_filename . $FILTERING_FILE_EXTENSION . $FASTA_FILE_EXTENSION;
    my $filtered_seq_filename             = $g_main_filename . $FILTERED_SEQ_FILE_EXTENSION;
    my $output_mask_filename1             = $g_main_filename . '.1' . $MASKING_FILE_EXTENSION;
    my $output_mask_filename2             = $g_main_filename . '.2' . $MASKING_FILE_EXTENSION;
    my $output_mask_filename3             = $g_main_filename . '.3' . $MASKING_FILE_EXTENSION;
    my $bsp_filename1                     = $g_main_filename . '.1' . $BSP_MASKING_FILE_EXTENSION;
    my $bsp_filename2                     = $g_main_filename . '.2' . $BSP_MASKING_FILE_EXTENSION;
    my $masking_html_filename1            = $g_main_filename . '.1' . $MASKING_HTML_FILE_EXTENSION;
    my $masking_html_filename2            = $g_main_filename . '.2' . $MASKING_HTML_FILE_EXTENSION;
    my $filtered_html_filename            = $g_main_filename . $FILTERED_HTML_FILE_EXTENSION;

    return ($output_filename, $filtered_fasta_alignment_filename, $filtered_seq_filename, $output_mask_filename1, $output_mask_filename2, $output_mask_filename3, $bsp_filename1, $bsp_filename2, $masking_html_filename1, $masking_html_filename2, $filtered_html_filename);
}


=head2 checkQualityAlignment

B<Description>: print warnings if masked alignement is too short compared to the original one or
                if percentage of identity is too low

B<ArgsCount>: 2

=over 4

=item $input_filename: (string) (R)

input alignment file name (FASTA format).

=item $input_masked_filename: (string) (R)

input alignment file name (FASTA format).

=back

B<Return>: (string)

1 (print warnings if necessary)

=cut

sub checkQualityAlignment
{
    my ($input_filename, $input_masked_filename ) = @_;
    # parameters check
    if (2 != scalar @_)
    {
        confess "usage: checkQualityAlignment(input_alignment_name, input_masked_filename);";
    }
    # check if a filename was provided
    if (not $input_filename)
    {
        confess "Input file name is missing!";
    }

    # display stage start
    printStageStart("Alignments checking");

    open(F, ">>$g_main_filename$GreenPhylConfig::TRIMAL_LOG_FILENAME");

    # use Bio::AlignIO to read in the alignment
    my $str = Bio::AlignIO->new('-file' => "$input_filename", -format=>'fasta');
    my $aln = $str->next_aln();

    # use Bio::AlignIO to read in the masked alignment
    my $str2 = Bio::AlignIO->new('-file' => "$input_masked_filename", -format=>'fasta');
    my $aln2 = $str2->next_aln();

    my $warning = 0;
    if ($aln2->percentage_identity() <= 40)
    {
        printWarning('Percentage identity after masking is low:' .  $aln2->percentage_identity() .  "\n");
        print F 'Percentage identity after masking is low:' .  $aln2->percentage_identity() .  "\n";

        $warning = 1;
    }

    my $diff_length = ( ($aln2->length()) / ($aln->length()) );

    if ($diff_length <= 0.3)
    {
        printWarning('ratio length masked - non masked is low:' .  $diff_length .  "\n");
        print F 'ratio length masked - non masked is low:' .  $diff_length .  "\n";

        $warning = 1;
    }

    printStageEnd("OK") unless ($warning);
    print F 'Quality Alignment: OK' unless ($warning);
    close F;

    return 1;
}




=pod

=head2 executePhyML

B<Description>: execute PhyML phylogeny program on the given FASTA input
file and returns the phylogenetic trees.

B<ArgsCount>: 2

=over 4

=item $input_filename: (string) (R)

input alignment file name.

=item $root_tree: (boolean) (O)

tells if the returned tree should be rooted using midpoint rooting method.
Default: False (midpoint rooting)

=back

B<Return>: (list of 3 strings)

the first element is the phylogenetic tree, the second is the set of
bootstrap trees, both in Newick format and the last one is the distance
matrix.

B<Example>:

    my ($phylo_tree_filename, $bootstrap_trees_filename, $distance_matrix_filename) = executePhyML("alignment.fasta");

=cut

sub executePhyML
{

    my ($input_filename, $root_tree) = @_;
    # parameters check
    if ((1 > @_) || (2 < @_))
    {
        confess "usage: executePhyML(input_fasta_name, rooting);";
    }
    # check if a filename was provided
    if (not $input_filename)
    {
        confess "Input file name is missing!";
    }
    # make sure the infile exists and is readable
    if (not -r $input_filename )
    {
        confess "Unable to read file '$input_filename'!";
    }

    # display stage start
    printStageStart("Phylogeny (PhyML)");

    my ($output_tree_filename, $bootstrap_trees_filename, $distance_matrix_filename, $decoded_tree_filename) = getPhyMLOutputs($input_filename);
    my $phylogeny_tree_filename  = $input_filename . $PHYML_TREE_FILE_SUFFIX;

    # alignment_file    data_type    nb_data_sets    nb_bootstrap    subst_model    prop_invariable_sites    nb_subst_categories    gamma_parameter   starting_tree    optimise_topology__branch_lengths_and_rate_parameters  Tree topology search operation option
    # -i $in_align      -d aa        -n 1            -b 100          -m LG          -v estimate              -c 4                   -a estimate       default (BIONJ)  -o tlr -s BEST
    #$cmd = "$GreenPhylConfig::PHYML_CMD -i $input_filename -d aa -n 1 -b $GreenPhylConfig::BOOTSTRAP_COUNT -m LG -v e -c 4 -a e -o n -s BEST >$g_main_filename$GreenPhylConfig::PHYML_LOG_FILENAME 2>$g_main_filename$GreenPhylConfig::PHYML_ERROR_LOG_FILENAME";
    my $cmd;
    if ($GreenPhylConfig::SEQUENCE_TYPE eq 'nucleic')
    {
        $cmd = "$GreenPhylConfig::PHYML_CMD -i $input_filename -d nt -n 1 -s SPR -b $GreenPhylConfig::BOOTSTRAP_COUNT -m HKY85  -c 4 -a e -o tlr >$g_main_filename$GreenPhylConfig::PHYML_LOG_FILENAME 2>$g_main_filename$GreenPhylConfig::PHYML_ERROR_LOG_FILENAME"  ;
    }
    else
    {
        $cmd = "$GreenPhylConfig::PHYML_CMD -i $input_filename -d aa -n 1 -s SPR -b $GreenPhylConfig::BOOTSTRAP_COUNT -m LG -v e -c 4 -a e -o tlr >$g_main_filename$GreenPhylConfig::PHYML_LOG_FILENAME 2>$g_main_filename$GreenPhylConfig::PHYML_ERROR_LOG_FILENAME"  
    }
    
    printDebug("COMMAND: $cmd");
    if ((-e $phylogeny_tree_filename)
        || (-e $bootstrap_trees_filename)
        || (-e $output_tree_filename)
        || (-e $distance_matrix_filename))
    {
        printWarning("PhyML: existing output files will be replaced!");
    }

    if (0 != system($cmd))
    {
        if (not $!)
        {
            confess "PhyML: see logs ($g_main_filename$GreenPhylConfig::PHYML_LOG_FILENAME and $g_main_filename$GreenPhylConfig::PHYML_ERROR_LOG_FILENAME) for errors.";
        }
        confess "PhyML: $!";
    }

    # check for rooting
    if ($root_tree)
    {
        # remove script file
        unlink "$g_output_dir/$PHYLIP_SCRIPT_FILENAME";
        # remove outtree file
        unlink "$g_output_dir/$PHYLIP_OUTTREE_FILENAME";
        # create script
        my $phylip_script_handle;
        open($phylip_script_handle, ">$g_output_dir/$PHYLIP_SCRIPT_FILENAME")
            or confess "Failed to create PHYLIP script file!";
        ($phylogeny_tree_filename) = ($phylogeny_tree_filename =~ m/([^\/]+)$/);
        print $phylip_script_handle "Y\n$phylogeny_tree_filename\nM\nW\nR\nQ\nQ\n";
        close($phylip_script_handle);

        # run retree
        $cmd = "cd $g_output_dir; $GreenPhylConfig::RETREE_CMD <$PHYLIP_SCRIPT_FILENAME >$g_main_filename$GreenPhylConfig::PHYLIP_LOG_FILENAME 2>&1; cd -";
        printDebug("Midpoint-rooting using Retree on \"$phylogeny_tree_filename\"...\nCOMMAND: $cmd");
        if (0 != system($cmd))
        {
            if (not $!)
            {
                confess "Retree: see logs ($g_main_filename$GreenPhylConfig::PHYLIP_LOG_FILENAME) for errors.";
            }
            confess "Retree: $!";
        }
        if (not -e "$g_output_dir/$PHYLIP_OUTTREE_FILENAME")
        {
            confess "Retree: failed to create rooted tree ($g_output_dir/$PHYLIP_OUTTREE_FILENAME)!";
        }
        # rename file
        if (!rename("$g_output_dir/$PHYLIP_OUTTREE_FILENAME", $output_tree_filename))
        {
            confess "ERROR: Failed to rename output rooted phylogeny tree file '$g_output_dir/$PHYLIP_OUTTREE_FILENAME' into '$output_tree_filename'!\n";
        }
    }
    else
    {
        # rename file
        if (!rename($phylogeny_tree_filename, $output_tree_filename))
        {
            confess "ERROR: Failed to rename output phylogeny tree file '$phylogeny_tree_filename' into '$output_tree_filename'!\n";
        }
    }


    # check if we got output trees
    if ((-s $output_tree_filename)
        && ((0 >= $GreenPhylConfig::BOOTSTRAP_COUNT)
            || ((-s $bootstrap_trees_filename) && (-s $distance_matrix_filename))
           )
       )
    {
        printStageEnd("OK");
    }
    else
    {
        printStageEnd("Failed!");
        if (-s $output_tree_filename)
        {
            if (0 < $GreenPhylConfig::BOOTSTRAP_COUNT)
            {
                # we got a phylogeny tree but a problem with bootstrap trees or distance matrix
                if (-s $bootstrap_trees_filename)
                {
                    # problem with the distance matrix file
                    if (-e $distance_matrix_filename)
                    {
                        # warn if we use a version of PhyML that does not output matrices
                        warn "Failed to compute phylogeny: distance matrix file '$distance_matrix_filename' is empty!";
                    }
                    else
                    {
                        # warn if we use a version of PhyML that does not output matrices
                        warn "Failed to compute phylogeny: distance matrix file '$distance_matrix_filename' is missing!";
                    }
                }
                elsif (-e $bootstrap_trees_filename)
                {
                    confess "Failed to compute phylogeny: bootstrap tree file '$bootstrap_trees_filename' is empty!";
                }
                else
                {
                    confess "Failed to compute phylogeny: bootstrap tree file '$bootstrap_trees_filename' is missing!";
                }
            }
        }
        elsif (-e $output_tree_filename)
        {
            # empty phylogeny tree
            confess "Failed to compute phylogeny: phylogeny tree file '$output_tree_filename' is empty!";
        }
        else
        {
            confess "Failed to compute phylogeny: phylogeny tree file '$output_tree_filename' is missing!";
        }
    }
    return ($output_tree_filename, $bootstrap_trees_filename, $distance_matrix_filename);
}




=pod

=head2 getPhyMLOutputs

B<Description>: returns the file names generated by PhyML phylogeny program.

B<ArgsCount>: 1

=over 4

=item $input_filename: (string) (R)

input alignment file name.


=back

B<Return>: (list of 4 strings)

the first element is the phylogenetic tree, the second is the set of
bootstrap trees, both in Newick format, then the distance matrix and finally
the phylogenetic tree with decoded names.

B<Example>:

    my ($phylo_tree_filename, $bootstrap_trees_filename, $distance_matrix_filename, $decoded_tree_filename) = getPhyMLOutputs("alignment.fasta");

=cut

sub getPhyMLOutputs
{

    my ($input_filename) = @_;
    # parameters check
    if (1 != @_)
    {
        confess "usage: getPhyMLOutputs(input_fasta_name);";
    }
    # check if a filename was provided
    if (not $input_filename)
    {
        confess "Input file name is missing!";
    }

    my $phylogeny_tree_filename  = $g_main_filename . $PHYLOGENY_FILE_SUFFIX;
    my $bootstrap_trees_filename = $input_filename . $PHYML_BOOTSTRAP_FILE_SUFFIX;
    my $distance_matrix_filename = $input_filename . $PHYML_DISTANCE_MATRIX_FILE_SUFFIX;
    my $decoded_tree_filename = $phylogeny_tree_filename . $OUTPUT_FILE_EXTENSION;

    return ($phylogeny_tree_filename, $bootstrap_trees_filename, $distance_matrix_filename, $decoded_tree_filename);
}




=pod

=head2 convertNewickToPhyloXML

B<Description>: convert the given Newick tree (binary and unrooted) into PhyloXML tree format.

B<ArgsCount>: 1-2

=over 4

=item $input_filename: (string) (R)

input tree file name (Newick format, binary and unrooted).

=item $output_filename: (string) (O)

the desired output tree file name (in PhyloXML format).

=back

B<Return>: (string)

the name of the PhyloXML tree file.

B<Example>:

    my $phyloxml_tree_filename = convertNewickToPhyloXML("tree.nwk");

=cut

sub convertNewickToPhyloXML
{
    my ($input_filename, $output_filename) = @_;

    # parameters check
    if ((1 > @_) || (2 < @_))
    {
        confess "usage: phyloxml_tree = convertNewickToPhyloXML(input_newick_tree)";
    }
    # check if a filename was provided
    if (not $input_filename)
    {
        confess "Input file name is missing!";
    }
    # make sure the infile exists and is readable
    if (not -r $input_filename)
    {
        confess "Unable to read file '$input_filename'!";
    }

    if (not $output_filename)
    {
        $output_filename = $input_filename . $PHYLOXML_FILE_EXTENSION;
    }
    my $rooted_tree_filename = $input_filename;
    $rooted_tree_filename =~s/$NEWICK_FILE_EXTENSION$//i;
    $rooted_tree_filename .= $NEWICK_EXTENDED_FILE_EXTENSION;

    printDebug("Newick to PhyloXML\n");
    if (-e $output_filename)
    {
        printWarning("NWK2XML: existing output file \"$output_filename\" will be replaced!");
        unlink $output_filename;
    }

    # remove old file
    if (-e $rooted_tree_filename)
    {
        printWarning("NWK2XML: existing output file \"$rooted_tree_filename\" will be replaced!");
        unlink $rooted_tree_filename;
    }

    # convert Newick tree into Phylo XML tree

    # -reformat standard newick to New Hampshire eXtended
    printDebug("Conversion from Newick to New Hampshire eXtended...\n");
    my $newick_handle;
    open($newick_handle, "<$input_filename")
        or confess "Failed to open input Newick file ($input_filename)!";
    my $newick_data = join('', <$newick_handle>);
    close($newick_handle);
    printDebug("Newick:\n    $newick_data\n");

    my $nhx_handle;
    open($nhx_handle, ">$rooted_tree_filename")
        or confess "Failed to create NHX file ($rooted_tree_filename)!";
    $newick_data =~ s/([\w.]+)_([A-Z]+):([\d.]+)/$1:$3\[\&\&NHX:S=$2\]/gs;
    printDebug("NHX:\n    $newick_data\n");
    print $nhx_handle $newick_data;
    close($nhx_handle);

    # -run forester "phyloxml_converter" converter
    my $cmd = "$GreenPhylConfig::PHYLOXML_CONV_CMD -f=gn -i $rooted_tree_filename $output_filename >$g_main_filename$GreenPhylConfig::PHYLOXML_LOG_FILENAME 2>$g_main_filename$GreenPhylConfig::PHYLOXML_ERROR_LOG_FILENAME";
    printDebug("Tree conversion using PhyloXML converter:\nCOMMAND: $cmd");
    if (0 != system($cmd))
    {
        confess "PhyloXML: see logs ($g_main_filename$GreenPhylConfig::PHYLOXML_LOG_FILENAME and $g_main_filename$GreenPhylConfig::PHYLOXML_ERROR_LOG_FILENAME) for errors.";
    }
    elsif (not -r $output_filename)
    {
        # make sure the converted file exists and is readable
        confess "Unable to read file '$output_filename'! Tree conversion failed!";
    }

    # convert <scientific_name> tags into <code> tags as NHX format doesn't offer a field for taxonomy code, so we used scientific_name instead
    printDebug("Replace <scientific_name> tags by <code> tags...\n");
    my $outxml_handle;
    open($outxml_handle, "<$output_filename")
        or confess "Failed to open output PhyloXML tree file ($output_filename) for reading!";
    my $xml_data = join('', <$outxml_handle>);
    close($outxml_handle);

    open($outxml_handle, ">$output_filename")
        or confess "Failed to open output PhyloXML tree file ($rooted_tree_filename) for writing!";
    $xml_data =~ s/scientific_name/code/gi;
    print $outxml_handle $xml_data;
    close($outxml_handle);

    return $output_filename;
}



=pod

=head2 executeSDI

B<Description>: execute the SDI program to infer gene duplications on a gene tree. Both gene and species tree have to be binary.

B<ArgsCount>: 1-2

=over 4

=item $input_filename: (string) (R)

input tree file name (Newick format).

=item $root_sdi: (boolean) (O)

tells if the returned tree should be rooted using SDI R method.
Default: False (no SDI R rooting)
rooted by minimizing the mapping cost L (and also the sum of duplications)
and by minimizing tree height (-mh + -ml )

=back

B<Return>: (string)

the name of the rooted tree file.

B<Example>:

    my $rooted_tree_filename = executeSDI("unrooted_tree.nwk");

=cut

sub executeSDI
{
    my ($input_filename, $root_sdi) = @_;

    # parameters check
    if ((1 > @_) || (2 < @_))
    {
        confess "usage: executeSDI(input_newick_tree)";
    }
    # check if a filename was provided
    if (not $input_filename)
    {
        confess "Input file name is missing!";
    }
    # make sure the infile exists and is readable
    if (not -r $input_filename)
    {
        confess "Unable to read file '$input_filename'!";
    }

    # if more than 1 lines, the files contains multiple trees (bootstrapped trees)
    my $count = qx/wc -l < $input_filename/;
    confess "wc failed: $?" if $?;
    chomp($count);
    my $phyloxml_suffix = ($count < 2 ) ? "_tree" : "_trees";

    my $xml_tree_filename    = $g_main_filename . $phyloxml_suffix . $PHYLOXML_FILE_EXTENSION; # _tree.xml or _trees.xml
    #my $rooted_tree_filename = $g_output_dir . '/' . $SDI_OUTPUT_FILE_EXT;
    my $rooted_tree_filename = $g_main_filename . $phyloxml_suffix . '.sdi' . $PHYLOXML_FILE_EXTENSION;

    printStageStart("SDI");

    # remove old files
    if (-e $xml_tree_filename)
    {
        printWarning("SDI: existing output file \"$xml_tree_filename\" will be replaced!");
        unlink $xml_tree_filename;
    }
    if (-e $rooted_tree_filename)
    {
        printWarning("SDI: existing output file \"$rooted_tree_filename\" will be replaced!");
        unlink $rooted_tree_filename;
    }

    # convert Newick tree into Phylo XML tree
    $xml_tree_filename = convertNewickToPhyloXML($input_filename, $xml_tree_filename);

    # get absolute path to change current working directory because output file
    # cannot be specified in SDI_R
    sub rewriteJarPathToAbsolute
    {
        my $cmd = shift;
        # find jar path
        my ($jar_path) = ($cmd =~ m/\s["']?(\S+\.jar)["']?\s/i);
        # get absolute path to jar
        my $jar_abs_path = abs_path($jar_path);
        # replace current jar path by its absolute path in the command line
        $cmd =~ s/$jar_path/$jar_abs_path/;
        return $cmd;
    }

    my $current_directory      = abs_path(getcwd());
    $xml_tree_filename         = abs_path($xml_tree_filename);
    my $tree_of_life_file_path = abs_path($GreenPhylConfig::TREE_OF_LIFE);
    my $sdi_r_cmd              = rewriteJarPathToAbsolute($GreenPhylConfig::SDI_R_CMD);
    my $sdi_cmd                = rewriteJarPathToAbsolute($GreenPhylConfig::SDI_CMD);

    chdir($g_output_dir);

    if (-e $SDI_OUTPUT_FILE)
    {
        printWarning("SDI: existing output file \"$SDI_OUTPUT_FILE\" will be replaced!");
        unlink $SDI_OUTPUT_FILE;
    }

    # prepare SDI command line (note: appends logs to the ones from the tree conversion)
    my $cmd;

    if ($root_sdi)
    {
        $cmd = "$sdi_r_cmd -ml -mh $xml_tree_filename $tree_of_life_file_path >>$g_main_filename$GreenPhylConfig::SDI_LOG_FILENAME 2>>$g_main_filename$GreenPhylConfig::SDI_ERROR_LOG_FILENAME";
        $rooted_tree_filename = $xml_tree_filename;
        $rooted_tree_filename =~ s/$PHYLOXML_FILE_EXTENSION/$GreenPhylConfig::SDI_EXTENSION$PHYLOXML_FILE_EXTENSION/i;
    }
    else
    {
        $cmd = "$sdi_cmd -b $xml_tree_filename $tree_of_life_file_path $rooted_tree_filename >>$g_main_filename$GreenPhylConfig::SDI_LOG_FILENAME 2>>$g_main_filename$GreenPhylConfig::SDI_ERROR_LOG_FILENAME";
    }
    printDebug("COMMAND: $cmd");

    if (0 != system($cmd))
    {
        confess "SDI: see logs ($g_main_filename$GreenPhylConfig::SDI_LOG_FILENAME and $g_main_filename$GreenPhylConfig::SDI_ERROR_LOG_FILENAME) for errors.";
    }

    # go back to working directory (because output file cannot be specified with SDI_R)
    chdir($current_directory);

    # check if we got an output rooted tree
    if (-s $rooted_tree_filename)
    {

        printStageEnd("OK");
    }
    else
    {
        printStageEnd("Failed!");
        if (-e $rooted_tree_filename)
        {
            confess "Failed to root tree: rooted tree file '$rooted_tree_filename' was empty!";
        }
        else
        {
            confess "Failed to root tree: rooted tree file '$rooted_tree_filename' was missing!";
        }
    }
    return $rooted_tree_filename;
}



=pod

=head2 getSDIOutput

B<Description>: returns the file names generated by SDI.

B<ArgsCount>: 1

=over 4

=item $input_filename: (string) (R)

input tree file name (Newick format).

=back

B<Return>: (string)

the name of the rooted tree file.

B<Example>:

    my $rooted_tree_filename = getSDIOutput("unrooted_tree.nwk");

=cut

sub getSDIOutput
{
    my ($input_filename) = @_;

    # parameters check
    if (1 != @_)
    {
        confess "usage: getSDIOutput(input_newick_tree)";
    }
    # check if a filename was provided
    if (not $input_filename)
    {
        confess "Input file name is missing!";
    }
    my @name = split (/_/, $input_filename);
    return $name[0]. '_tree' . $SDI_OUTPUT_FILE_EXT;
}


sub getSDIBootstrapOutput
{
    my ($input_filename) = @_;

    # parameters check
    if (1 != @_)
    {
        confess "usage: getSDIOutput(input_newick_tree)";
    }
    # check if a filename was provided
    if (not $input_filename)
    {
        confess "Input file name is missing!";
    }

    my @name = split (/_/, $input_filename);
    return  $name[0]. '_trees' . $SDI_OUTPUT_FILE_EXT;
}


=pod

=head2 executeRIO

B<Description>: execute the RIO program to infer orthologs.

B<ArgsCount>: 4

=over 4

=item $dictionary_filename: (string) (R)

dictionary file name.

=item $rooted_tree_filename: (string) (R)

input rooted tree file name (PhyloXML format).

=item $bootstrap_trees_filename: (string) (R)

bootstrap trees file name (PhyloXML format).

=item $distance_matrix_filename: (string) (R)

distance matrix file name (PHYLIP format).

=back

B<Return>: (list of array of 2 strings)

each array of the list contains 2 string, the first one being the file name of RIO outputs and the
second one being RIO tree output.

B<Example>:

    my @files = executeRIO("dictionary.txt", "rooted_tree.nwk", "bootstrap_trees.nwk", "distance_matrix.txt");

=cut

sub executeRIO
{
    my ($dictionary_filename,
        $rooted_tree_filename,
        $bootstrap_trees_filename,
        $distance_matrix_filename) = @_;

    # parameters check
    if (4 != @_)
    {
        confess "usage: executeRIO(dictionary, rooted_newick_tree, bootstrap_trees, distance_matrix)";
    }
    # check if a filename was provided
    if ((not $dictionary_filename)
        || (not $rooted_tree_filename)
        || (not $bootstrap_trees_filename)
        || (not $distance_matrix_filename))
    {
        confess "Some input files names are missing!";
    }
    # make sure the infile exists and is readable
    if (not -r $dictionary_filename)
    {
        confess "Unable to read file '$dictionary_filename'!";
    }
    if (not -r $rooted_tree_filename)
    {
        confess "Unable to read file '$rooted_tree_filename'!";
    }
    if (not -r $bootstrap_trees_filename)
    {
        confess "Unable to read file '$bootstrap_trees_filename'!";
    }
    if (not -r $distance_matrix_filename)
    {
        confess "Unable to read file '$distance_matrix_filename'!";
    }

    printStageStart("RIO");

    my @output_files = ();
    my @queries = ();

    # Get queries
    my $dictionary_fh;
    open($dictionary_fh, "<$dictionary_filename")
        or confess "Failed to open dictionary file ($dictionary_filename)!";
    while (my $ligne = <$dictionary_fh>)
    {
        my ($query) = ($ligne =~ /\t(.+)*$SPECIES_CODE_REGEXP$/o);
        if ($query)
        {
            push @queries, $query;
        }
    }
    close($dictionary_fh);
    printDebug("Found " . @queries . " query sequence.");

    # create single distance matrix data
    my $single_distance_matrix_filename = $distance_matrix_filename . ".dist";
    my $distance_matrix_fh;
    open($distance_matrix_fh, "<$distance_matrix_filename")
        or confess "Failed to open distance matrix file ($distance_matrix_filename)!";
    my $distance_matrix_data = join('', <$distance_matrix_fh>);
    close($distance_matrix_fh);
    open($distance_matrix_fh, ">$single_distance_matrix_filename")
        or confess "Failed to open distance matrix file ($single_distance_matrix_filename)!";
    ($distance_matrix_data) = ($distance_matrix_data =~ m/^(.*?)\n\n/sm);
    # remove species codes
    $distance_matrix_data =~ s/^([\w.]+)$SPECIES_CODE_REGEXP /$1 /gsmo;
    print $distance_matrix_fh $distance_matrix_data;
    close($distance_matrix_fh);

    #if resume
    unlink "$g_main_filename$GreenPhylConfig::RIO_LOG_FILENAME" if (-e "$g_main_filename$GreenPhylConfig::RIO_LOG_FILENAME");
    unlink "$g_main_filename$GreenPhylConfig::RIO_ERROR_LOG_FILENAME" if (-e "$g_main_filename$GreenPhylConfig::RIO_ERROR_LOG_FILENAME");

    # Treatment
    my $dorio_output_filename_seed = $g_output_dir . '/' . $RIO_OUTPUT_FILE_SEED;
    foreach my $query (@queries)
    {
        # get a valid name for files
        my $query_name          = $query;
        #$query_name             =~ s/\W/_/g; # is it necessary
        my $query_output_data   = $dorio_output_filename_seed . $query_name . $RIO_DATA_OUTPUT_FILE_EXTENSION;

        printDebug("Working on $query\n data --> $query_output_data");

        # check if files exist
        if (-e $query_output_data)
        {
            printWarning("RIO: existing output file \"$query_output_data\" will be replaced!");
            unlink $query_output_data or die "cannot remove file $query_output_data: $!";
        }

        # add files to list
        push @output_files, [$query_output_data];

        my $cmd = "$GreenPhylConfig::DORIO_CMD M=$bootstrap_trees_filename N=$query S=$GreenPhylConfig::TREE_OF_LIFE O=$query_output_data T=$rooted_tree_filename D=$single_distance_matrix_filename P=13 L=30 p 1>>$g_main_filename$GreenPhylConfig::RIO_LOG_FILENAME 2>>$g_main_filename$GreenPhylConfig::RIO_ERROR_LOG_FILENAME";

        printDebug("COMMAND: $cmd");

        if (0 != system($cmd))
        {
            confess "RIO: $! see logs ($g_main_filename$GreenPhylConfig::RIO_LOG_FILENAME and $g_main_filename$GreenPhylConfig::RIO_ERROR_LOG_FILENAME) for errors.";
        }
    }

    foreach my $file_pair (@output_files)
    {
        if (not -s $file_pair->[0])
        {
            warn "Failed to infer orthologs for: file '$file_pair->[0]' is missing!";
        }
    }
    printStageEnd("OK");

    return @output_files;
}



=pod

=head2 executeRapGreen

B<Description>: execute the RapGreen program to root the tree and infer orthologs using a defined species tree.

B<ArgsCount>: 1

=over 1

=item $unrooted_tree_filename: (string) (R)

input rooted tree file name (PhyloXML format).

=back

B<Return>: (list of array of 2 strings)

each array of the list contains 2 strings, the first one being


B<Example>:

    my @files = executeRapGreen("unrooted_tree.nwk");

=cut

sub executeRapGreen
{

     my ($unrooted_tree_filename) = @_;

    # parameters check
    if (1 != @_)
    {
        confess "usage: executeRapGreen(unrooted_newick_tree)";
    }
    # check if a filename was provided
    if (not $unrooted_tree_filename)
    {
        confess "input file name is missing!";
    }
    # make sure the infile exists and is readable
    if (not -r $unrooted_tree_filename)
    {
        confess "Unable to read file '$unrooted_tree_filename'!";
    }

    printStageStart("RAP");
    
    # encode IDs
    my $rap_dictionary_filename = $g_main_filename . '_rap' . $DICTIONARY_FILE_EXTENSION;
    my $enocded_unrooted_tree_filename = &convertIDs({'conversion' => 'encode', 'format' => 'newick', 'keep_species_code' => 1}, $unrooted_tree_filename, $rap_dictionary_filename);

    my $enc_rooted_tree_filename      = $g_main_filename . '.enc' . $ROOTED_TREE_FILE_SUFFIX;
    my $rooted_tree_filename          = $g_main_filename . $ROOTED_TREE_FILE_SUFFIX;
    my $enc_reconcilied_tree_filename = $g_main_filename . '.enc' . $RECONCILIED_TREE_FILE_SUFFIX;
    my $reconcilied_tree_filename     = $g_main_filename . $RECONCILIED_TREE_FILE_SUFFIX;
    my $enc_stats_tree_filename       = $g_main_filename . '.enc' . $STATS_TREE_FILE_SUFFIX;
    my $stats_tree_filename           = $g_main_filename . $STATS_TREE_FILE_SUFFIX;
    my $enc_xml_tree_filename         = $g_main_filename . '.enc' . $XML_TREE_FILE_SUFFIX;
    my $xml_tree_filename             = $g_main_filename . $XML_TREE_FILE_SUFFIX;

    #if resume
    unlink "$g_main_filename$GreenPhylConfig::RAP_LOG_FILENAME" if (-e "$g_main_filename$GreenPhylConfig::RAP_LOG_FILENAME");
    unlink "$g_main_filename$GreenPhylConfig::RAP_ERROR_LOG_FILENAME" if (-e "$g_main_filename$GreenPhylConfig::RAP_ERROR_LOG_FILENAME");

    # remove old file
    if (-e $rooted_tree_filename)
    {
        printWarning("RAP: existing output file \"$rooted_tree_filename\" will be replaced!");
        unlink $rooted_tree_filename;
        unlink $enc_rooted_tree_filename;
    }

    my $cmd = "$GreenPhylConfig::RAP_CMD -g $enocded_unrooted_tree_filename -s $GreenPhylConfig::TREE_OF_LIFE -og $enc_rooted_tree_filename -or $enc_reconcilied_tree_filename -gt " . (0 < $GreenPhylConfig::BOOTSTRAP_COUNT ? '80' : '0.95') . " -st 10.0 -pt 0.00 -stats $enc_stats_tree_filename -phyloxml $enc_xml_tree_filename 1>$g_main_filename$GreenPhylConfig::RAP_LOG_FILENAME 2>$g_main_filename$GreenPhylConfig::RAP_ERROR_LOG_FILENAME";
    printDebug("COMMAND: $cmd");
    
    if (0 != system($cmd))
    {
        confess "RAP: $! see logs ($g_main_filename$GreenPhylConfig::RAP_LOG_FILENAME and $g_main_filename$GreenPhylConfig::RAP_ERROR_LOG_FILENAME) for errors.";
    }


    if ((-s $enc_rooted_tree_filename) && (-s $enc_reconcilied_tree_filename))
    {
        printStageEnd("OK");
    }
    else
    {
        printStageEnd("Failed!");
    }
     
    &convertIDs({'conversion' => 'decode'}, $enc_rooted_tree_filename,      $rap_dictionary_filename, $rooted_tree_filename);
    &convertIDs({'conversion' => 'decode'}, $enc_reconcilied_tree_filename, $rap_dictionary_filename, $reconcilied_tree_filename);
    &convertIDs({'conversion' => 'decode'}, $enc_stats_tree_filename,       $rap_dictionary_filename, $stats_tree_filename);
    &convertIDs({'conversion' => 'decode'}, $enc_xml_tree_filename,         $rap_dictionary_filename, $xml_tree_filename);

    return ($rooted_tree_filename, $reconcilied_tree_filename);
}

=pod

=head2 getRapGreenOutput

B<Description>: returns the file names generated by RapGreen.

B<ArgsCount>: 1

=over 4

=item $input_filename: (string) (R)

input tree file name (Newick format).

=back

B<Return>: (string)

the name of the rooted tree file.

B<Example>:

    my $rooted_tree_filename = getRapGreenOutput("unrooted_tree.nwk");

=cut

sub getRapGreenOutput
{
    my ($input_filename) = @_;

    # parameters check
    if (1 != @_)
    {
        confess "usage: getSDIOutput(input_newick_tree)";
    }
    # check if a filename was provided
    if (not $input_filename)
    {
        confess "Input file name is missing!";
    }
    my $rooted_tree_filename      = "$g_main_filename" . $ROOTED_TREE_FILE_SUFFIX;
    my $reconcilied_tree_filename = "$g_main_filename" . $RECONCILIED_TREE_FILE_SUFFIX;
    return ($rooted_tree_filename, $reconcilied_tree_filename);
}


=pod

=head2 convertIDs

B<Description>: encodes/decodes sequences names.

B<ArgsCount>: 3-6

=over 4

=item $parameters: (hash ref) (R)

Hash of parameters. Possible keys-values:
'conversion'        => 'encode': encode names
                    => 'decode': decode names
'format'            => input file format. Can be one of:
                       'FASTA', 'NEWICK'
'keep_species_code' => true: keep species code in names

=item $input_filename: (string) (R)

Input filename to encode/decode.

=item $dictionary_filename: (string) (R)

dictionary filename to use.

=item $output_filename: (string) (U)

output filename to use. If omitted, $input_filename is used to generate
the output filename by appending $OUTPUT_FILE_EXTENSION to $input_filename.

=item $name_match_re: (string) (U)

Used with the 'decode' methode. $name_match_re is a regular expression
that can be used to find and match the dictionnary name. Do not use matching
parenthesis or the boundary characters ^ and $.
Default: 'SEQ\d{7}'

=item $name_replace_sub: (sub) (O)

Used with the 'decode' methode. $name_replace_sub is a reference to a function
that can process the matched name and return an altered 'on-purpose' version
of that name. The matched name is the first argument of the function, the
second is a hash reference of the dictionnary to use (keys are SEQxxxx names
and values are original sequence names).

=back

B<Return>: (string)

the encoded/decoded filename.

B<Example>:

    my $encoded_filename = &convertIDs({'conversion' => 'encode', 'format' => 'fasta'}, $input_filename, $dictionary_filename);
    ...
    my $decoded_filename = &convertIDs({'conversion' => 'decode'}, $encoded_filename, $dictionary_filename);

=cut

sub convertIDs
{
    my ($paramaters, $input_filename, $dictionary_filename, $output_filename, $name_match_re, $name_replace_sub) = @_;

    # parameters check
    if ((3 > @_) || (6 < @_))
    {
        confess "usage: convertIDs(conversion_mode, input_filename, dictionary_filename[, output_filename])";
    }
    # check if a filename was provided
    if (not $input_filename)
    {
        confess "Need an input file!";
    }
    # make sure the infile exists and is readable
    if (not -r $input_filename)
    {
        confess "Input file is not readable!";
    }

    # check if an output filename has been specified
    if ((not defined $output_filename) || ('' eq $output_filename))
    {
        # append extension to input filename
        $output_filename = $input_filename . $OUTPUT_FILE_EXTENSION;
    }

    my $debug_message = <<"___1187_DEBUG_MESSAGE___";
convertIDs:
conversion     : $paramaters->{'conversion'}
input file     : $input_filename
output file    : $output_filename
dictionary file: $dictionary_filename
___1187_DEBUG_MESSAGE___

    printDebug($debug_message);

    if ($paramaters->{'conversion'} =~ m/^encode$/i)
    {
        # encode...
        if ($paramaters->{'format'} =~ m/^fasta$/i)
        {
            my $cpt        = 0;
            my $seq_stream = Bio::SeqIO->new(-file => $input_filename,     -format => 'Fasta');
            my $seq_out    = Bio::SeqIO->new(-file => ">$output_filename", -format => 'Fasta');
            # create dictionary...
            my $dictionary_handle;
            open($dictionary_handle, ">$dictionary_filename")
                or confess "Failed to create dictionary file ($dictionary_filename)!";
            while (my $seq = $seq_stream->next_seq())
            {
                ++$cpt;
                my $seq_id = $seq->id();
                my $new_id = '';
                if ($paramaters->{'keep_species_code'})
                {
                    $new_id = $seq_id;
                    # remove sequence name part and just keep species info
                    $new_id =~ s/^.*($SPECIES_CODE_REGEXP)\s*$/$1/o;
                }
                $new_id = "SEQ" . "0" x (7 - length($cpt)) . $cpt . $new_id;
                $seq->id($new_id);
                $seq_out->write_seq($seq);
                # if species code has been kept in sequence name, do not store it in dictionary
                if ($paramaters->{'keep_species_code'})
                {
                    $new_id =~ s/$SPECIES_CODE_REGEXP$//o;
                    $seq_id =~ s/$SPECIES_CODE_REGEXP$//o;
                }
                print $dictionary_handle "$new_id\t$seq_id\n";
            }
            close($dictionary_handle);
        }
        elsif ($paramaters->{'format'} =~ m/^newick$/i)
        {
            my $cpt        = 0;
            my $tree_stream = Bio::TreeIO->new(-file => $input_filename, -format => 'newick');
            my $tree_out    = Bio::TreeIO->new(-file => ">$output_filename", -format => 'newick');
            # create dictionary...
            my $dictionary_handle;
            open($dictionary_handle, ">$dictionary_filename")
                or confess "Failed to create dictionary file ($dictionary_filename)!";
            my $tree = $tree_stream->next_tree();
            if (!$tree)
            {
                confess "Newick file '$input_filename' doesn't see to contain a valid tree!\n";
            }
            foreach my $taxon ($tree->get_leaf_nodes())
            {
                ++$cpt;
                my $taxon_id = $taxon->id();
                my $new_id = '';
                if ($paramaters->{'keep_species_code'})
                {
                    $new_id = $taxon_id;
                    # remove sequence name part and just keep species info
                    $new_id =~ s/^.*($SPECIES_CODE_REGEXP)$/$1/o;
                }
                $new_id = "SEQ" . "0" x (7 - length($cpt)) . $cpt . $new_id;
                $taxon->id($new_id);
                if ($paramaters->{'keep_species_code'})
                {
                    $new_id   =~ s/$SPECIES_CODE_REGEXP$//o;
                    $taxon_id =~ s/$SPECIES_CODE_REGEXP$//o;
                }
                print $dictionary_handle "$new_id\t$taxon_id\n";
            }
            $tree_out->write_tree($tree);
            close($dictionary_handle);
        }
        elsif ($paramaters->{'format'})
        {
            confess "Unsupported format '$paramaters->{'format'}'!\n";
        }
        else
        {
            confess "No input format supplied!\n";
        }
    }
    elsif ($paramaters->{'conversion'} =~ m/^decode$/i)
    {
        # decode...
        if (!-e $dictionary_filename)
        {
            confess "Need a dictionary file!";
        }

        my %dictionary; # dictionary hash
        # parse dictionary file
        my $dictionary_handle;
        open($dictionary_handle, $dictionary_filename)
            or confess "Failed to open dictionary file ($dictionary_filename)!";
        while(<$dictionary_handle>)
        {
            chomp;
            my @line = split("\t");
            $dictionary{$line[0]} = $line[1];
        }
        close($dictionary_handle);

        # decode file
        my ($input_handle, $output_handle);
        open($input_handle, $input_filename)
            or confess "Failed to open input file for name decoding!";
        open($output_handle, ">$output_filename")
            or confess "Failed to create output file for name decoding!";
        # check how to match dictionnary names
        if ((!defined($name_match_re)) || ('' eq $name_match_re))
        {
            # use default
            $name_match_re = 'SEQ\d{7}';
        }
        # check if original name should be processed
        if (!$name_replace_sub)
        {
            # use default
            $name_replace_sub = sub {my ($match, $dictionary) = @_; return $dictionary->{$match};};
        }

        while (<$input_handle>)
        {
            s/($name_match_re)/&$name_replace_sub($1, \%dictionary)/ge;
            print {$output_handle} $_;
        }
        close($output_handle);
        close($input_handle);
    }
    else
    {
        confess "ERROR: no conversion mode selected!\n";
    }

	if (not -s $output_filename)
    {
        confess "Failed to " . $paramaters->{'conversion'} . " IDs";
	}
    return $output_filename;
}


=pod

=head2 validate_mail

B<Description>: validate a given e-mail address.

B<ArgsCount>: 1

=over 4

=item $target_email: (string) (R)

the target e-mail address to check.

=back

B<Return>: (boolean)

returns true (1) if the mail address is valid, false (0) otherwise.

B<Caller>:

general

B<Example>:

    my $target_mail = 'user@server.com';
    validate_mail($target_mail) or confess 'Not a valid e-mail address!\n';

=cut

sub validate_mail
{
    my ($target_email) = @_;

    # check parameters
    if (1 != @_)
    {
        confess "usage: validate_mail(target_email) or confess 'invalid e-mail address!';";
    }

    # RFC used:
    # -local mailbox: http://tools.ietf.org/html/rfc2822
    # -domain: http://tools.ietf.org/html/rfc1034
    my ($local_mailbox, $domain) = ($target_email =~ m/^(.*)@([a-zA-Z][a-zA-Z0-9\-]*(?:\.[a-zA-Z0-9\-]{1,63})*[a-zA-Z0-9])$/);
    if ((not $local_mailbox)
        || (not $domain)
        || (64 < length($local_mailbox))
        || (255 < length($domain)))
    {
        # invalid local address or domain
        return 0;
    }
    if (($local_mailbox !~ m/^[a-zA-Z0-9!~&'=_#\-\^\$\|\*\+\?\{\}\%\/\`]+(?:\.[a-zA-Z0-9!~&'=_#\-\^\$\|\*\+\?\{\}\%\/\`]+)*$/)
        && ($local_mailbox !~ m/^"(?:[^\\]*(?:\\"[^"]*\\")*(?:\\[^"])*)*"$/))
    {
        # invalid local address
        return 0;
    }
    return 1;
}



=pod

=head2 send_mail

B<Description>: Sends an e-mail.

B<ArgsCount>: hash

=over 4

=item email: (string) (R)

Target e-mail address.

=item subject: (string) (O)

Mail subject.

=item body: (string) (O)

Mail body.

=back

B<Return>: (nothing)

B<Example>:

    send_mail('email' => 'v.guignon@cgiar.org', 'subject' => 'Done', 'body' => 'done!');

=cut

sub send_mail
{
    my (%parameters) = @_;
    
    my $target_email = $parameters{'email'};
    my $mail_subject = $parameters{'subject'};
    my $mail_body    = $parameters{'body'};
    
    if (!$target_email)
    {
        confess "ERROR: send_mail: no target e-mail address provided!\n";
    }
    elsif (!validate_mail($target_email))
    {
        confess "ERROR: Invalid e-mail address!\n";
    }

    if (!$mail_subject)
    {
        $mail_subject = "[GreenPhyl] Pipeline $$";
    }
    $mail_subject =~ s/[\n\r]+//g;

    if (!$mail_body)
    {
        $mail_body = "-GreenPhyl-";
    }
    $mail_body .= "\n";

    my $sendmail_h;
    if (open($sendmail_h, "|/usr/sbin/sendmail -t"))
    {
        print {$sendmail_h} "To: $target_email\n";
        print {$sendmail_h} "From: $DEFAULT_ADMIN_EMAIL\n";
        print {$sendmail_h} "Subject: $mail_subject\n\n";
        print {$sendmail_h} $mail_body;
        close($sendmail_h);
    }
    else
    {
        warn "WARNING: Failed to send notrification e-mail!\n$!\n";
    }
}




# Script options
#################

=pod

=head1 OPTIONS

    run_pipeline.pl -help | -f <family_id> -i <fasta_in_file>
        [-r alignment | masking | phylogeny | rooting | rio]
        [-e alignment | masking | phylogeny | rooting | rio]
        [-autoresume] [-skip_bad_masking] [-output_dir <directory>] [-rap|-norap]
        [-email <email>]

=head2 Parameters

=over 4

=item B<-help>:

display this help message

=item B<-f family_id> (string):

the family identifier (name)

=item B<-show_config>:

display current configuration and exists.

=item B<-i fasta_in_file> (string):

the FASTA (multiple set of sequences) input file name.

=item B<-r alignment | masking | phylogeny | rooting | rio> (string):

resume the pipeline from the specified step.

=item B<-e alignment | masking | phylogeny | rooting | rio> (string):

end the pipeline after the specified step.

=item B<-auto-resume> (flag):

Skip steps when all requested outputs are available.

=item B<-skip_bad_masking> (flag):

When a masking process failed (generated gap-only sequences, no site,...) the
original alignment will be used instead of the curated one.
Default: stop on alignment error.

=item B<-output_dir> (string):

Path to the directory to use for output. By default, the pipeline creates a
subdirectory using the family name in $GreenPhylConfig::OUTPUT_PATH.
If '-output_dir' option is used, the directory specified will contain the
output files directly (they won't be output in a subdirectory having the family
name).

=item B<-email> (string):

e-mail address that will receive an e-mail once the pipeline has finished its
job.

=item B<-rap> (flag):

Use RAP instead of SDI for rooting.

=back

Defaults:
 bootstrap           = -4
 alignment           = mafft
 Improve_alignment   = trimal
 tree                = phyml
 rooting             = rap

=cut

# CODE START
#############

++$| if !$|;
umask 0022;

# display welcome message
print <<"___1291_WELCOME_MESSAGE___";

######################################
#                                    #
#        Welcome to $GREENPHYL_NAME        #
#                                    #
######################################
                                 v0.7

___1291_WELCOME_MESSAGE___

# options processing
my ($help, $man, $load, $input_filename, $input_family_name, $query, $use_rap,
    $skip_masking_failure, $skip_masking, $output_dir, $email, $show_config,
    $clear_previous_analyses);

# default: use rap
$use_rap = 1;

# parse options and print usage if there is a syntax error.
GetOptions(
    'help|h|?'         => \$help,
    'man'              => \$man,
    'debug'            => \$g_debug,
    'resume|r=s'       => \$g_resume,
    'autoresume|a'     => \$g_autoresume,
    'end|e=s'          => \$g_end,
    'infile|in|i=s'    => \$input_filename,
    'family|f=s'       => \$input_family_name,
    'rap!'             => \$use_rap,
    'skip_bad_masking' => \$skip_masking_failure,
    'skip_masking'     => \$skip_masking,
    'output_dir=s'     => \$output_dir,
    'email=s'          => \$email,
    'show_config'      => \$show_config,
    'clear'            => \$clear_previous_analyses,
) or pod2usage(1);
if ($help) { pod2usage(0); }
if ($man) { pod2usage( -verbose => 2 ); }

if ($show_config)
{
    GreenPhylConfig::displayConfig();
    exit(0);
}

# check if parameters were not specified using the dash-name syntax
$input_family_name = shift() unless $input_family_name;
$input_filename    = shift() unless $input_filename;

# options checking
if (not defined $input_filename)
{
    printError("$GREENPHYL_NAME: no fasta file specified!");
    pod2usage(1);
}
elsif (!-e $input_filename)
{
    printError("$GREENPHYL_NAME: fasta_in_file does not exists! [$input_filename]");
    pod2usage(1);
}

if (($g_resume ne $STEP_ALIGNMENT)
    && ($g_resume ne $STEP_HMM)
    && ($g_resume ne $STEP_MASKING)
    && ($g_resume ne $STEP_PHYLOGENY)
    && ($g_resume ne $STEP_TREE_ROOTING)
    && ($g_resume ne $STEP_RESAMPLE_INFERENCE_OF_ORTHOLOGY))
{
    printError("$GREENPHYL_NAME: invalid step specified (for resume)!");
    pod2usage(1);
}

if (($g_end ne $STEP_ALIGNMENT)
    && ($g_end ne $STEP_HMM)
    && ($g_end ne $STEP_MASKING)
    && ($g_end ne $STEP_PHYLOGENY)
    && ($g_end ne $STEP_TREE_ROOTING)
    && ($g_end ne $STEP_RESAMPLE_INFERENCE_OF_ORTHOLOGY))
{
    printError("$GREENPHYL_NAME: invalid step specified (for end)!");
    pod2usage(1);
}



#sequence count
$g_sequences_count = `grep '>' $input_filename | wc -l`;

# init output directory...
if ($output_dir)
{
    # remove trailing slash
    $g_output_dir = $output_dir; # init output path name
    $g_output_dir =~ s/\/+$//;
    $g_output_dir .= '/' . $input_family_name
}
else
{
    $g_output_dir = $GreenPhylConfig::OUTPUT_PATH . '/' . $input_family_name; # init output path name
}
if (!-e $g_output_dir)
{
    if (system("mkdir -p '$g_output_dir'"))
    {
        warn "WARNING: failed to create directory '$g_output_dir' with the appropriate access rights!\n"; # init access rights
    }
    else
    {
       chmod(02775, "$g_output_dir");
    }
}

# clear previous analyses if needed
if ($clear_previous_analyses)
{
    printDebug("Clear previous analyses");
    if (system("rm -f $g_output_dir/*"))
    {
        warn "WARNING: failed to clear directory '$g_output_dir'!\n$!\n";
    }
}

system("cp $input_filename $g_output_dir/$input_family_name$INPUT_FILE_EXTENSION"); # copy input file
$input_filename = "$g_output_dir/$input_family_name$INPUT_FILE_EXTENSION";
# extract file name seed
($g_main_filename) = ($input_filename =~ m/^(.*[\/\\]?[^.\/\\]+)\.[^\/\\]+$/); # extract file name without extension (get name character until the first dot)
if (not $g_main_filename)
{
    ($g_main_filename) = "$g_output_dir/$input_family_name";
}
if (not $g_main_filename)
{
    printError("$GREENPHYL_NAME: failed to extract the main part of the file name!");
}

if ($g_debug)
{
    GreenPhylConfig::displayConfig();
    printDebug("File name seed: $g_main_filename");
}

if ($g_autoresume)
{
    print "NOTE: Auto-Resume Mode ON\n";
}

print "Start working on $input_filename...\n";

# encode species names
my $dictionary_filename = $g_main_filename . $DICTIONARY_FILE_EXTENSION;
my $encoded_filename = &convertIDs({'conversion' => 'encode', 'format' => 'fasta'}, $input_filename, $dictionary_filename);
printDebug("input_file: $input_filename\nencoded_file: $encoded_filename");

# it allows to restart analysis at different steps if prevous steps were performed successfully
my %step_name_to_step_number = (
                         $STEP_ALIGNMENT                       => 0,
                         $STEP_HMM                             => 1,
                         $STEP_MASKING                         => 2,
                         $STEP_PHYLOGENY                       => 3,
                         $STEP_TREE_ROOTING                    => 4,
                         $STEP_RESAMPLE_INFERENCE_OF_ORTHOLOGY => 5,
                     );

if ($g_resume ne $STEP_ALIGNMENT)
{
    printDebug("Resume to step $g_resume (" . $step_name_to_step_number{$g_resume} . ")");
}

if ($g_end ne $STEP_RESAMPLE_INFERENCE_OF_ORTHOLOGY)
{
    printDebug("Will stop after step $g_end (" . $step_name_to_step_number{$g_end} . ")");
}

## MAFFT
my $alignment_filename = '';
if (($step_name_to_step_number{$STEP_ALIGNMENT} >= $step_name_to_step_number{$g_resume})
    && (!$g_autoresume || (!-s getMAFFTOutput($encoded_filename))))
{
    $alignment_filename = &executeMAFFT($encoded_filename);
    $g_autoresume = 0; # reset auto-resume to re-run next steps if auto-resume was set
}
else
{
    printDebug("Skip alignment");
    $alignment_filename = &getMAFFTOutput($encoded_filename);
}

## HMMBUILD
my $hmm_thread;
if (($step_name_to_step_number{$STEP_HMM} >= $step_name_to_step_number{$g_resume})
    && ($step_name_to_step_number{$STEP_HMM} <= $step_name_to_step_number{$g_end})
    && (!$g_autoresume || (!-s getHMMBUILDOutput($encoded_filename))))
{
    ($hmm_thread) = threads->create(
            sub
            {
                my ($l_stdout, $l_stderr, $stdoutput_h, $stderror_h);
                open($stdoutput_h, '>', \$l_stdout) or confess "Failed to redirect STDOUT! $!\n";
                open($stderror_h, '>', \$l_stderr) or confess "Failed to redirect STDERR! $!\n";
                $g_stdout = $stdoutput_h;
                $g_stderr = $stderror_h;
                #executeHMMBUILD($alignment_filename);
                close($stdoutput_h);
                close($stderror_h);
                $gs_thread_messages .= "STDOUT:\n" . ($l_stdout || '') . "\nSTDERR:\n" . ($l_stderr || '') . "\n";
            }
        );
}
else
{
    printDebug("Skip HMM");
    #getHMMBUILDOutput($alignment_filename);
}

## MASKING
my $masked_filtered_alignment_filename = '';
my $filtered_fasta_alignment_filename  = '';
my $filtered_seq_filename              = '';
my $output_mask_filename1              = '';
my $output_mask_filename2              = '';
my $output_mask_filename3              = '';
my $bsp_filename1                      = '';
my $bsp_filename2                      = '';
my $masking_html_filename1             = '';
my $masking_html_filename2             = '';
my $filtered_html_filename             = '';
my $filtered_sequence_ratio            = 0;
if ($skip_masking)
{
    ($masked_filtered_alignment_filename,
     $filtered_fasta_alignment_filename,
     $filtered_seq_filename,
     $output_mask_filename1,
     $output_mask_filename2,
     $output_mask_filename3,
     $bsp_filename1,
     $bsp_filename2,
     $masking_html_filename1,
     $masking_html_filename2,
     $filtered_html_filename,
    ) = &getMaskingAndFilteringOutputs($alignment_filename);
    if (system("rm $masked_filtered_alignment_filename"))
    {
        cluck "WARNING: Failed to remove '$masked_filtered_alignment_filename'\n";
    }
    # open source alignment file
    my $input_fasta_align = Bio::AlignIO->new(
        '-file'   => $alignment_filename,
        '-format' => 'fasta',
        '-alphabet' => 'protein',
    );
    # create PHYLIP version of the alignment
    my $output_phylip_align = Bio::AlignIO->new(
        '-file'   => ">$masked_filtered_alignment_filename",
        '-format' => 'phylip',
        '-alphabet' => 'protein',
    );
    while (my $aln = $input_fasta_align->next_aln)
    {
        $output_phylip_align->write_aln($aln);
    }
}
elsif (($step_name_to_step_number{$STEP_MASKING} >= $step_name_to_step_number{$g_resume})
    && ($step_name_to_step_number{$STEP_MASKING} <= $step_name_to_step_number{$g_end})
    && (!$g_autoresume
        || (!-s (getMaskingAndFilteringOutputs($alignment_filename))[0])
        || (!-s (getMaskingAndFilteringOutputs($alignment_filename))[1])))
{
    if ($skip_masking_failure)
    {
        printDebug("If an error occurs during masking, the original aligment will be used.");
    }
    else
    {
        printDebug("Abort on error during masking.");
    }

    try
    {
        # my $refined_alignment_filename = &executeRefiner($alignment_filename);
        # my $GreenPhylConfig::HMM_ALIGN_CMDment_filename     = &executeHmm($refined_alignment_filename);
        ($masked_filtered_alignment_filename,
         $filtered_fasta_alignment_filename,
         $filtered_seq_filename,
         $output_mask_filename1,
         $output_mask_filename2,
         $output_mask_filename3,
         $bsp_filename1,
         $bsp_filename2,
         $masking_html_filename1,
         $masking_html_filename2,
         $filtered_html_filename,
         $filtered_sequence_ratio,
        ) = executeMaskingAndFiltering($alignment_filename);
        if (-s $filtered_seq_filename)
        {
            my $decoded_filename = &convertIDs({'conversion' => 'decode'}, $filtered_seq_filename, $dictionary_filename);
            unlink($filtered_seq_filename);
            rename($decoded_filename, $filtered_seq_filename);
        }
        # truncate names to 20 characters in HTML outputs
        sub truncate_name
        {
            my ($match, $dictionary) = @_;
            $match =~ s/<\/span>\s{10}//;
            # check if seq name should be truncated
            if (20 < length($match))
            {
                return substr($dictionary->{$match}, 0, 9) . '..' . substr($dictionary->{$match}, length($match) - 9) . '</span>';
            }
            else
            {
                return sprintf('%-20s</span>', substr($dictionary->{$match}, 0, 20));
            }
        }
        my $match = 'SEQ\d{7}<\/span>\s{10}';
        # decode SEQXXX names in HTML files
        if (-s $masking_html_filename1)
        {
            my $decoded_filename = &convertIDs({'conversion' => 'decode'}, $masking_html_filename1, $dictionary_filename, undef, $match, \&truncate_name );
            unlink($masking_html_filename1);
            rename($decoded_filename, $masking_html_filename1);
        }
        if (-s $masking_html_filename2)
        {
            my $decoded_filename = &convertIDs({'conversion' => 'decode'}, $masking_html_filename2, $dictionary_filename, undef, $match, \&truncate_name );
            unlink($masking_html_filename2);
            rename($decoded_filename, $masking_html_filename2);
        }
        if (-s $filtered_html_filename)
        {
            my $decoded_filename = &convertIDs({'conversion' => 'decode'}, $filtered_html_filename, $dictionary_filename, undef, $match, \&truncate_name);
            unlink($filtered_html_filename);
            rename($decoded_filename, $filtered_html_filename);
        }
        
        # check filtration level and stop if the threshold is reached
        if ($filtered_sequence_ratio > $GreenPhylConfig::FILTERING_THRESHOLD)
        {
            warn "Filtering threshold ($GreenPhylConfig::FILTERING_THRESHOLD) has been reached! Aborting pipeline!\n";
        }

        #compare aln from mafft and after masking
        #printDebug("Quality alignment assessment");
        #checkQualityAlignment($alignment_filename, $masked_filtered_alignment_filename);
    }
    otherwise
    {
        my $error = shift();
        if ($skip_masking_failure)
        {
            cluck $error . "\nWARNING: Masking error ignored: original alignment will be used instead of the curated one (for '$g_output_dir')!\n";
                ($masked_filtered_alignment_filename,
                 $filtered_fasta_alignment_filename,
                 $filtered_seq_filename,
                 $output_mask_filename1,
                 $output_mask_filename2,
                 $output_mask_filename3,
                 $bsp_filename1,
                 $bsp_filename2,
                 $masking_html_filename1,
                 $masking_html_filename2,
                 $filtered_html_filename,
                ) = &getMaskingAndFilteringOutputs($alignment_filename);
            if (system("rm $masked_filtered_alignment_filename"))
            {
                confess "ERROR: Failed to remove '$masked_filtered_alignment_filename'\n";
            }
            if (system("cp $alignment_filename $masked_filtered_alignment_filename"))
            {
                confess "ERROR: failed to copy alignment file '$alignment_filename' to '$masked_filtered_alignment_filename'!\n";
            }
        }
        else
        {
            confess $error;
        }
    };

    # creates file for jalview
    &convertIDs({'conversion' => 'decode'}, $filtered_fasta_alignment_filename, $dictionary_filename, "$filtered_fasta_alignment_filename.aln");
    $g_autoresume = 0; # reset auto-resume to re-run next steps if auto-resume was set
}
else
{
    printDebug("Skip masking");
    ($masked_filtered_alignment_filename,
     $filtered_fasta_alignment_filename,
     $filtered_seq_filename,
     $output_mask_filename1,
     $output_mask_filename2,
     $output_mask_filename3,
     $bsp_filename1,
     $bsp_filename2,
     $masking_html_filename1,
     $masking_html_filename2,
     $filtered_html_filename,
    ) = &getMaskingAndFilteringOutputs($alignment_filename);
}

### PHYML
my ($phylogeny_tree_filename, $bootstrap_trees_filename, $distance_matrix_filename, $root_status, $decoded_tree_filename) = ('', '', '', '', '');
if (($step_name_to_step_number{$STEP_PHYLOGENY} >= $step_name_to_step_number{$g_resume})
    && ($step_name_to_step_number{$STEP_PHYLOGENY} <= $step_name_to_step_number{$g_end})
    && (!$g_autoresume
        || (!-s (getPhyMLOutputs($masked_filtered_alignment_filename))[0])
        || (!$use_rap && (!-s (getPhyMLOutputs($masked_filtered_alignment_filename))[2])) ))
{
    # use 0 instead of 1 to skip rooting by retree
    $root_status = 0;
    ($phylogeny_tree_filename, $bootstrap_trees_filename, $distance_matrix_filename) = &executePhyML($masked_filtered_alignment_filename, $root_status);
    $decoded_tree_filename = &convertIDs({'conversion' => 'decode'}, $phylogeny_tree_filename, $dictionary_filename);
    $g_autoresume = 0; # reset auto-resume to re-run next steps if auto-resume was set
}
else
{
    printDebug("Skip phylogeny");
    ($phylogeny_tree_filename, $bootstrap_trees_filename, $distance_matrix_filename, $decoded_tree_filename) = &getPhyMLOutputs($masked_filtered_alignment_filename);
}

### SDI/RAP
my $rooted_tree_filename      = '';
my $rooted_trees_filename     = '';
my $reconcilied_tree_filename = '';

if (($step_name_to_step_number{$STEP_TREE_ROOTING} >= $step_name_to_step_number{$g_resume})
    && ($step_name_to_step_number{$STEP_TREE_ROOTING} <= $step_name_to_step_number{$g_end})
    && (!$g_autoresume
        || ($use_rap?((!-s (getRapGreenOutput($masked_filtered_alignment_filename))[0]) || (!-s (getRapGreenOutput($masked_filtered_alignment_filename))[1]))
                    :(!-s getSDIBootstrapOutput($masked_filtered_alignment_filename)))))
{
    # to be or not to be (rooted), that is the question...
    $root_status = ($root_status) ? 0 : 1;

    # decode species names and root tree
    if ($use_rap)
    {
        # execute Rap Green
        ($rooted_tree_filename, $reconcilied_tree_filename) = &executeRapGreen($decoded_tree_filename);
    }
    else
    {
        # run SDI
        $rooted_tree_filename       = &executeSDI($decoded_tree_filename, $root_status);

        # decode species names and root bootstrapped tree
        my $decoded_trees_filename  = &convertIDs({'conversion' => 'decode'}, $bootstrap_trees_filename, $dictionary_filename);
        $rooted_trees_filename      = &executeSDI($decoded_trees_filename, $root_status);
    }
    $g_autoresume = 0; # reset auto-resume to re-run next steps if auto-resume was set
}
else
{
    printDebug("Skip rooting");
    if ($use_rap)
    {
        ($rooted_tree_filename, $reconcilied_tree_filename)  = &getRapGreenOutput($phylogeny_tree_filename . $OUTPUT_FILE_EXTENSION);
    }
    else
    {
        $rooted_tree_filename  = &getSDIOutput($phylogeny_tree_filename . $OUTPUT_FILE_EXTENSION);
        $rooted_trees_filename = &getSDIBootstrapOutput($phylogeny_tree_filename . $OUTPUT_FILE_EXTENSION)
    }
}



### RIO
if ($use_rap)
{
    printDebug("Using rap --> skip RIO");
}
elsif (5 > $GreenPhylConfig::BOOTSTRAP_COUNT)
{
    printWarning("Not enough bootstrap ($GreenPhylConfig::BOOTSTRAP_COUNT): skipping RIO step!");
}
else
{
    # note: no auto-resume for RIO: too complicated to check, just run again the step...
    # use RIO
    if (($step_name_to_step_number{$STEP_RESAMPLE_INFERENCE_OF_ORTHOLOGY} >= $step_name_to_step_number{$g_resume})
        && ($step_name_to_step_number{$STEP_RESAMPLE_INFERENCE_OF_ORTHOLOGY} <= $step_name_to_step_number{$g_end}))
    {
        # decode species names
        my $decoded_bootstrap_filename = &convertIDs({'conversion' => 'decode'}, $bootstrap_trees_filename, $dictionary_filename);
        my $decoded_distmat_filename   = &convertIDs({'conversion' => 'decode'}, $distance_matrix_filename, $dictionary_filename);

        #
        my $bootstrap_fh;
        open($bootstrap_fh, "<$decoded_bootstrap_filename")
            or confess "Failed to open bootstrap trees file ($decoded_bootstrap_filename)!";
        my $bootstrap_trees = join('', <$bootstrap_fh>);
        close($bootstrap_fh);

        my $bootstrap_trees_filename = $rooted_trees_filename; #output from sdi_r

        &executeRIO($dictionary_filename, $rooted_tree_filename, $bootstrap_trees_filename, $decoded_distmat_filename);
        $g_autoresume = 0; # reset auto-resume to re-run next steps if auto-resume was set
    }
}

# wait for HMM thread to end before exiting
if ($hmm_thread)
{
    $hmm_thread->join();
}

# display thread messages
print "\nThread messages:\n";
print ($gs_thread_messages || 'none');
print "\n";

# give write access to group
foreach my $filename (glob("$g_output_dir/*"))
{
    if (!-d "$g_output_dir/$filename")
    {
        chmod(0664, "$g_output_dir/$filename");
        # set group
        #+val: disabled because we use setgid flag (ie. chmod g+s on parent directory which belongs to greenphyl group)
        #chown([getpwnam(getlogin())]->[2], $GREENPHYL_GROUP_ID, $filename);
    }
}

# send email if needed
if ($email)
{
    send_mail(
        'email'   => $email,
        'subject' => '[GreenPhyl] Pipeline $$: done',
        'body'    => "Pipeline ($$) done!",
    );
}

exit(0);

# CODE END
###########




=pod

=head1 AUTHORS

Valentin GUIGNON (CIRAD, Bioversity-France), valentin.guignon@cirad.fr, V.Guignon@cgiar.org
Mathieu ROUARD (Bioversity-France), m.rouard@cgiar.org
Matthieu CONTE (Bioversity-France), m.conte@cgiar.org
Jean-François DUFAYARD (CIRAD), jeanfrancois.dufayard@gmail.com

Date 18/07/2011

=head1 SEE ALSO

GreenPhyl publications

=cut
