#!/usr/bin/env perl

package CTAT_fus_pipe_runner;

use strict;
use warnings;
use Carp;

use lib ($ENV{EUK_MODULES});
use Pipeliner;
use Process_cmd;
use File::Basename;

####
sub run_fusion_pipe {
    my ($left_fq, $right_fq, $analysis_id) = @_;

    unless (-s $left_fq && -s $right_fq && $analysis_id) {
        confess "Error, set params accordingly";
    }
    my $CTAT_GENOME_LIB = $ENV{CTAT_GENOME_LIB};
    unless ($CTAT_GENOME_LIB) {
        die "Error, must set CTAT_GENOME_LIB env var";
    }
    
    my $capture_dir = "retain_outputs_${analysis_id}";
   
    my $checkpoint_dir = "checkpoints_dir";
    unless (-d $checkpoint_dir) {
        mkdir $checkpoint_dir or die "Error, cannot mkdir $checkpoint_dir";
    }

    unless (-d $capture_dir) {
        mkdir $capture_dir or die "Error, cannot mkdir $capture_dir";
    }
    
    my $pipeliner = new Pipeliner(-verbose => 2);
    $pipeliner->set_checkpoint_dir($checkpoint_dir);

    ###################
    ## Run STAR-Fusion
    
    my $cmd = "/home/unix/bhaas/GITHUB/CTAT_FUSIONS/STAR-Fusion/STAR-Fusion --left_fq $left_fq --right_fq $right_fq -O STAR-Fusion --genome_lib_dir $CTAT_GENOME_LIB";
    $pipeliner->add_commands( new Command($cmd, "star-fusion.ok"));
            
    $pipeliner->add_commands( new Command("gzip -c STAR-Fusion/Chimeric.out.junction > STAR-Fusion/Chimeric.out.junction.gz", "gzip_chimeric_out.ok"));
    $pipeliner->add_commands( new Command("samtools view -Sb STAR-Fusion/Chimeric.out.sam > STAR-Fusion/Chimeric.out.bam", "chimeric_sam_to_bam.ok"));
    

    # capture outputs:
    &append_capture_outputs_cmds($pipeliner, $capture_dir, ["STAR-Fusion/star-fusion.fusion_candidates.final", 
                                                            "STAR-Fusion/star-fusion.fusion_candidates.final.abridged.FFPM",
                                                            "STAR-Fusion/Chimeric.out.junction.gz",
                                                            "STAR-Fusion/Chimeric.out.bam"
                                 ] );
    
    
    ######################
    ## Run FusionInspector
    
    $cmd = "/home/unix/bhaas/GITHUB/CTAT_FUSIONS/FusionInspector/FusionInspector --left_fq $left_fq --right_fq $right_fq --out_dir FusionInspector --out_prefix finspector --prep_for_IGV --fusions STAR-Fusion/star-fusion.fusion_candidates.final.abridged.FFPM --genome_lib_dir $CTAT_GENOME_LIB --aligner_path /home/unix/bhaas/GITHUB/STAR/bin/Linux_x86_64_static/STAR "; 
    
    $pipeliner->add_commands( new Command($cmd, "FusionInspector.ok") );
    
    &append_capture_outputs_cmds($pipeliner, $capture_dir, ["FusionInspector/finspector.fusion_predictions.final",
                                                            "FusionInspector/finspector.fusion_predictions.final.abridged.FFPM",
                                                            "FusionInspector/finspector.fa",
                                                            "FusionInspector/finspector.gtf",
                                                            "FusionInspector/finspector.bed",
                                                            "FusionInspector/finspector.spanning_reads.bam",
                                                            "FusionInspector/finspector.spanning_reads.bam.bai",
                                                            "FusionInspector/finspector.spanning_reads.bam.bed",
                                                            "FusionInspector/finspector.junction_reads.bam",
                                                            "FusionInspector/finspector.junction_reads.bam.bai",
                                                            "FusionInspector/finspector.junction_reads.bam.bed",
                                                            "FusionInspector/finspector.igv.FusionJuncSpan",
                                                            "FusionInspector/cytoBand.txt"] );
    
    
    
     
    $cmd = "/home/unix/bhaas/GITHUB/CTAT_FUSIONS/DISCASM/DISCASM --aligned_bam STAR-Fusion/Aligned.sortedByCoord.out.bam --chimeric_junctions STAR-Fusion/Chimeric.out.junction --left_fq $left_fq --right_fq $right_fq --out_dir DISCASM_OI_ASM --denovo_assembler Oases --normalize_reads";
    $pipeliner->add_commands( new Command($cmd, "DISCASM.ok") );

    $pipeliner->add_commands(new Command("gzip DISCASM_OI_ASM/oases_out_dir/oases.transcripts.fa", "gzip_discasm_fasta.ok"));
    
    my $norm_left_fq = "DISCASM_OI_ASM/" . basename($left_fq) . ".extracted.fq.normalized_K25_C50_pctSD200.fq";
    my $norm_right_fq = "DISCASM_OI_ASM/" . basename($right_fq) . ".extracted.fq.normalized_K25_C50_pctSD200.fq";
    
    $pipeliner->add_commands(new Command("gzip $norm_left_fq $norm_right_fq", "gzip_norm_fq.ok"));
                
    &append_capture_outputs_cmds($pipeliner, $capture_dir, ["DISCASM_OI_ASM/oases_out_dir/oases.transcripts.fa.gz",
                                                            "$norm_left_fq.gz",
                                                            "$norm_right_fq.gz"] );
    
    

    $pipeliner->run();
        


}
####
sub append_capture_outputs_cmds {
    my ($pipeliner, $capture_dir, $outputs_to_capture_aref) = @_;

    foreach my $output_to_capture (@$outputs_to_capture_aref) {
        my $cmd = "cp $output_to_capture $capture_dir/";
        
        # create checkpoint based on output filename
        my $token = $output_to_capture;
        $token =~ s/\W/_/g;
        my $chkpt_file = "$token.ok";
        
        $pipeliner->add_commands( new Command($cmd, $chkpt_file) );
    }

    return;
}
    






main: {

    if (basename($0) eq basename(__FILE__)) {
        
        my $usage = "\n\n\tusage: $0 left.fq.gz right.fq.gz analysis_id\n\n";
        
        my $left_fq = $ARGV[0] or die $usage;
        my $right_fq = $ARGV[1] or die $usage;
        my $analysis_id = $ARGV[2] or die $usage;
        
        &run_fusion_pipe($left_fq, $right_fq, $analysis_id);
        

        exit(0);
        
    }
}


1; #EOM
