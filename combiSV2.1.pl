#!/usr/bin/env perl
######################################################
#         SOFTWARE COPYRIGHT NOTICE AGREEMENT        #
#  Copyright (C) {2020-2021}  {Nicolas Dierckxsens}  #
#              All Rights Reserved                   #
#         See file LICENSE for details.              #
######################################################
#                    combiSV 2.1
#           nicolasdierckxsens@hotmail.com
use strict;
use Getopt::Long;

print "\n\n-----------------------------------------------";
print "\ncombiSV\n";
print "Version 2.1\n";
print "Author: Nicolas Dierckxsens, (c) 2020-2021\n";
print "-----------------------------------------------\n\n";

my $input_pbsv = "";
my $input_sniffles = "";
my $input_cutesv = "";
my $input_nanovar = "";
my $input_svim = "";
my $output_file = "";
my $output_file2 = "";
my $high_recall = "";
my $input_nanosv = "";
my $count_tools = "2";
my $min_coverage = '3';
my $length_match = '800';
my $min_SV_length = '50';

GetOptions (
            "pbsv=s" => \$input_pbsv,
            "sniffles=s" => \$input_sniffles,
            "cutesv=s" => \$input_cutesv,
            "nanovar=s" => \$input_nanovar,
            "svim=s" => \$input_svim,
            "nanosv=s" => \$input_nanosv,
            "c=s" => \$min_coverage,
            "o=s" => \$output_file,
            ) or die "Incorrect usage!\n";

if ($input_sniffles eq "" && $input_cutesv eq "" && $input_pbsv eq "" && $input_svim eq "")
{
    print "\n\nUsage:     perl combiSV2.1.pl -pbsv <pbsv_output.vcf> -sniffles <sniffles_output.vcf> -cutesv <cutesv_output.vcf> -nanovar <nanovar_output.vcf> -svim <svim_output.vcf> -nanosv <nanosv_output.vcf>\n\n";
    print "\nOPTIONAL ARGUMENTS\n";
    print "-c     minimum coverage of variation allele [default = 3]\n";
    print "-o     name of the output files\n";
}

my $total_tools = '0';
if ($min_coverage eq "")
{
    $min_coverage = '3';
}
if ($min_coverage =~ m/\d+/)
{
}
else
{
    die "\n\nWARNING: minimum coverage of variation allele has to be an integer\n\n";
}
if ($min_coverage < 2)
{
    $min_coverage = '2';
    print "\n\nWARNING: minimum coverage of variation allele too low, has been set to 2\n\n";
}
if ($input_pbsv ne "")
{
    open(INPUT_PBSV, $input_pbsv) or die "\n\nCan't open pbsv's vcf file $input_pbsv, $!\n\n";
    $total_tools++;
}
if ($input_sniffles ne "")
{
    open(INPUT_SNIFFLES, $input_sniffles) or die "\n\nCan't open Sniffles' vcf file $input_sniffles, $!\n\n";
    $total_tools++;
}
if ($input_cutesv ne "")
{
    open(INPUT_CUTESV, $input_cutesv) or die "\n\nCan't open cuteSV's vcf file $input_cutesv, $!\n\n";
    $total_tools++;
}
if ($input_nanovar ne "")
{
    open(INPUT_NANOVAR, $input_nanovar) or die "\n\nCan't open NanoVar's vcf file $input_nanovar, $!\n\n";
    $total_tools++;
}
if ($input_svim ne "")
{
    open(INPUT_SVIM, $input_svim) or die "\n\nCan't open SVIM's vcf file $input_svim, $!\n\n";
    $total_tools++;
}
if ($input_nanosv ne "")
{
    open(INPUT_NANOSV, $input_nanosv) or die "\n\nCan't open NanoSV's vcf file $input_nanosv, $!\n\n";
    $total_tools++;
}

if ($input_sniffles eq "" && $input_cutesv eq "" && $input_pbsv eq "" && $input_svim eq "")
{
    die "\n\nError: A Sniffles, pbsv, SVIM or cuteSV input is mandatory.\n\n";
}

if ($output_file eq "")
{
    $output_file = "combiSV";
    $output_file2 = "simplified_combiSV";
}
else
{
    $output_file2 = $output_file;
    $output_file = "simplified_".$output_file2;  
}
my $last_four = substr $output_file, -4, 4;
if ($last_four ne ".vcf")
{
    $output_file .= ".vcf";
    $output_file2 .= ".vcf";
}

if ($high_recall eq "")
{
    $high_recall = "1";
}
if ($high_recall ne "" && $high_recall ne "1" && $high_recall ne "2")
{
    die "\n\nIncorrect usage of -s parameter, should be '1' or '2'\n\n";
}

if ($total_tools > '5')
{
    $count_tools = '3';
    $high_recall = '2';
}
elsif ($total_tools <= 5)
{
    $count_tools = '2';
    $high_recall = '2';
}


open(COMBINED, ">" .$output_file) or die "\nCan't open file $output_file, $!\n";
open(COMBINED2, ">" .$output_file2) or die "\nCan't open file $output_file2, $!\n";

my %SVs;
my %count;
my %pbsv;
my %sniffles;
my %cutesv;
my %nanovar;
my %mixed_types;
my $lenght_margin1 = 0.70;
my $lenght_margin2 = 1.3;

my $count = '0';

if ($input_pbsv ne "")
{
    while (my $line = <INPUT_PBSV>)
    {
        chomp($line);
        my $first_nuc = substr $line, 0, 1;
        if ($first_nuc ne "#")
        {
            my @list = split /\t/, $line;
            
            my $info = $list[7];
            my $SVLEN = "";
            my $type = "";
            my @info = split /;/, $info;
            my @HAP = split /:/, $list[9];
            my @HAP2 = split /,/, $HAP[1];
            my @HAP3 = split /,/, $HAP[3];
            
            my $chr = $list[0];
            if ($list[0] =~ m/chr(\d+|X|Y)/)
            {
                $chr = $1;
            }
            
            foreach my $info_tmp (@info)
            {
                my $first_five = substr $info_tmp, 0, 5;
                if ($info_tmp =~ m/SVLEN=>*-*(\d+)/)
                {
                    $SVLEN = $1;
                }
                elsif ($first_five eq "SVTYP")
                {
                    $type = substr $info_tmp, 7;
                }
            }
            if ($SVLEN eq "")
            {
                $SVLEN = ".";
            }
            if ($SVLEN < 0)
            {
                $SVLEN *= -1;
            }
            my $haplo = $HAP[0];
            
            my $converted_line = $chr."\t".$list[1]."\t".$SVLEN."\t".$type."\t".$haplo;
    
            if ($type ne "BND" && $type ne "cnv" && $list[6] eq "PASS" && ($type ne "DUP" || $high_recall eq "2" || $haplo eq "1/1") && ($HAP2[1] >= $min_coverage ||
               ($HAP2[1] >= $min_coverage-1 && $HAP2[1]/$HAP[2] > 0.3)) && $haplo ne "0/0" && $HAP2[1]/$HAP[2] > 0.3 && ($SVLEN eq "." || $SVLEN > 45))
            {                  
                $SVs{$chr}{$list[1]} = $converted_line;
                $count{$chr}{$list[1]} = '1';
                if ($type eq "INV" || (($type eq "DEL" || $haplo eq "1/1") && $high_recall eq "2") && ($HAP2[1]+$HAP2[2]) > 9 && (($HAP3[2] > 0 && $HAP3[3] > 0) || $HAP3[2] eq "" || $HAP3[3] eq ""))
                {
                    $count{$chr}{$list[1]} = '2';
                }
                if (($type eq "INV" || $type eq "INS") && $haplo eq "1/1" && $total_tools < 5)
                {
                    $count{$chr}{$list[1]} = $count_tools;
                }
                $pbsv{$chr}{$list[1]} = $line;
            }
        }
    }
}

my %pbsv2;
my %cutesv2;
my %sniffles2;
my %nanovar2;
my %svim2;

if ($input_cutesv ne "")
{
CUTESV: while (my $line = <INPUT_CUTESV>)
    {
        chomp($line);      
        my $first = substr $line, 0, 1;
        my $prev_match2 = "";
        
        if ($first ne "#")
        {
            my @list = split /\t/, $line;
            my $pos = $list[1];
            my $length;
            my $type;
            my $v = '1';
            my $min = '1';
            
            my $info = $list[7];      
            my @info = split /;/, $info;
            my @HAP = split /:/, $list[9];
            
            foreach my $info_tmp (@info)
            {
                my $first_five = substr $info_tmp, 0, 5;
                if ($first_five eq "SVLEN")
                {
                    $length = substr $info_tmp, 6;
                }
                elsif ($first_five eq "SVTYP")
                {
                    $type = substr $info_tmp, 7;
                }
            }
            if ($length eq "")
            {
                $length = ".";
            }
            if ($length < 0)
            {
                $length *= -1;
            }
            my $chr = $list[0];
            if ($list[0] =~ m/chr(\d+|X|Y)/)
            {
                $chr = $1;
            }
            my $converted_line = $chr."\t".$list[1]."\t".$length."\t".$type."\t".$HAP[0];
            
            if (exists ($cutesv{$chr}{$list[1]}))
            {
                next CUTESV;
            }
            
            if ($list[6] eq "PASS" && ($length eq "." || $length > 45) && ($HAP[2] >= $min_coverage || ($HAP[2] >= $min_coverage-1 && $HAP[2]/($HAP[2]+$HAP[1]) > 0.3))
                && $HAP[2]/($HAP[2]+$HAP[1]) > 0.12 && $type ne "BND" && $HAP[0] ne "./.")
            {               
                my $prev_match = "";
                if (exists($pbsv2{$chr}{$list[1]}))
                {
                    $prev_match = "yes";
                }
                
                if (exists($pbsv{$chr}{$list[1]}) && exists($SVs{$chr}{$list[1]}) && $prev_match eq "")
                {
                    my $count_tmp = $count{$chr}{$list[1]};

                    my @list2 = split /\t/, $SVs{$chr}{$list[1]};
                    my $new_pos = $list2[1];
                    my $new_length = $list2[2];
                    my $new_type = $list2[3];
                    my $new_haplo = $list2[4];
                    
                    $pbsv2{$chr}{$list[1]} = undef;
    
                    $count{$chr}{$list[1]} = $count_tmp+1;
                    if ($type eq "INS" && $list2[3] ne $type && $type ne ".")
                    {
                        $new_type = $type;        
                    }
                    if ($HAP[0] ne "." && ($type eq "INV" || $type eq "INS"))
                    {
                        $new_haplo = $HAP[0];
                    }
                    if ($type eq "INV" || $list2[2] eq ".")
                    {                                                                    
                        if ($length ne ".")
                        {
                            $new_length = $length;
                        }
                        if ($HAP[0] ne ".")
                        {
                            $new_haplo = $HAP[0];
                        }
                    }
                    delete $SVs{$chr}{$pos};
                    delete $count{$chr}{$pos};
                    my $line2 = $list2[0]."\t".$new_pos."\t".$new_length."\t".$new_type."\t".$new_haplo;
                    $SVs{$list2[0]}{$new_pos} = $line2;
                    $count{$chr}{$new_pos} = $count_tmp+1;
                    $cutesv{$chr}{$list[1]} = $line;
                    
                    if ($type eq "INS")
                    {
                        $count{$chr}{$new_pos} = $count_tools;
                    }
                    next CUTESV;
  
                }
                elsif ($input_pbsv ne "")
                {
                    my $score_prev = "";
                    my $pos_prev = "";
POS_ALMOST2da:                  
                    while ($v < $length_match)
                    {
POS_ALMOST2d:           my $pos_tmp = ($min*$v)+$pos;
     
                        if (exists($pbsv2{$chr}{$pos_tmp}))
                        {
                            $prev_match2 = "yes";
                        }
                        
                        if (exists($pbsv{$chr}{$pos_tmp}) && exists($SVs{$chr}{$pos_tmp}) && $prev_match2 eq "")
                        {                          
                            my $count_tmp = $count{$chr}{$pos_tmp};
                            my @list2 = split /\t/, $SVs{$chr}{$pos_tmp};
                            
                            if (exists($cutesv{$chr}{$pos_tmp}) && $list2[3] eq $type && ($list2[2] > $length*$lenght_margin1 && $list2[2] < $length*$lenght_margin2))
                            {
                                next CUTESV;
                            }
                            else
                            {                                               
                                my $score = '0';
                                if ($list2[3] eq $type || ($type eq "INS" && $list2[3] eq "DUP") || ($type eq "DUP" && $list2[3] eq "INS"))
                                {
                                    $score += 1;
                                }
                                my $g = (800-$v)/800;
                                $score += $g;
                                my $score_tmp = '0';
                                if ($list2[2] ne "." && $length ne ".")
                                {
                                    if ($list2[2] >= $length)
                                    {
                                        $score_tmp = 1-(($list2[2]-$length)/$length);                                           
                                        if ($score_tmp < 0)
                                        {
                                            $score_tmp = '0';
                                        }
                                    }
                                    else
                                    {
                                        $score_tmp = 1-(($length-$list2[2])/$length);
                                        if ($score_tmp < 0)
                                        {
                                            $score_tmp = '0';
                                        }
                                    }
                                }
                                $score += $score_tmp;

                                if ($score_prev ne "" && $score_prev > $score && $score_prev ne "no")
                                {
                                    $pos_tmp = $pos_prev;
                                    $count_tmp = $count{$chr}{$pos_tmp};
                                    @list2 = split /\t/, $SVs{$chr}{$pos_tmp};
                                }

                                if ($score > 2 || ($score_prev ne "" && $score_prev > $score) || $score_prev eq "no")
                                {
                                    
                                    if ($score_prev eq "no")
                                    {
                                        print $SVs{$chr}{$pos_tmp}."\n";
                                        print $line."\n";
                                    }
                                    my $new_pos = $pos_tmp;
                                    my $new_length = $list2[2];
                                    my $new_type = $list2[3];
                                    my $new_haplo = $list2[4];
                                    
                                    if ($type eq "INS" && $list2[3] ne $type)
                                    {                              
                                        $new_type = $type;
                                    }
                                    if ($HAP[0] ne "." && ($type eq "INV" || $type eq "INS"))
                                    {
                                        $new_haplo = $HAP[0];
                                    }
                                    if ($type eq "INV" || $list2[2] eq ".")
                                    {                                                                    
                                        if ($length ne ".")
                                        {
                                            $new_length = $length;
                                        }
                                        $new_pos = $pos;
                                    }
        
                                    delete $SVs{$chr}{$pos_tmp};
                                    delete $count{$chr}{$pos_tmp};
                                    my $line2 = $list2[0]."\t".$new_pos."\t".$new_length."\t".$new_type."\t".$new_haplo;
                                    $SVs{$list2[0]}{$new_pos} = $line2;
                                    $count{$chr}{$new_pos} = $count_tmp+1;
                                    $cutesv{$chr}{$new_pos} = $line;
                                    $pbsv2{$chr}{$new_pos} = undef;
                                    
                                    if ($type eq "INS")
                                    {
                                        $count{$chr}{$new_pos} = $count_tools;
                                    }
                                    
                                    if ($new_pos eq $pos)
                                    {
                                        if (exists($pbsv{$chr}{$pos_tmp}))
                                        {
                                            my $tmp = $pbsv{$chr}{$pos_tmp};
                                            delete $pbsv{$chr}{$pos_tmp};
                                            $pbsv{$chr}{$new_pos} = $tmp;
                                        }
                                        if (exists($cutesv{$chr}{$pos_tmp}))
                                        {
                                            my $tmp = $cutesv{$chr}{$pos_tmp};
                                            delete $cutesv{$chr}{$pos_tmp};
                                            $cutesv{$chr}{$new_pos} = $tmp;
                                        }
                                    }
                                    next CUTESV;
                                }
                                elsif ($min eq '1')
                                {                   
                                    $min = '-1';
                                    if ($score > $score_prev)
                                    {
                                        $score_prev = $score;
                                        $pos_prev = $pos_tmp;
                                    }
                                    goto POS_ALMOST2d;
                                }
                                else
                                {
                                    $min = '1';
                                    if ($score > $score_prev)
                                    {
                                        $score_prev = $score;
                                        $pos_prev = $pos_tmp;
                                    }
                                }
                                $v++;
                                goto POS_ALMOST2d;               
                            }
                        }
                        elsif ($min eq '1')
                        {                   
                            $min = '-1';
                            goto POS_ALMOST2d;
                        }
                        else
                        {
                            $min = '1';
                        }
                        $v++;
                    }
                    if ($score_prev ne "" && $score_prev ne "no" && $score_prev > 1.6)
                    {
                        $score_prev = "no";
                        $v = '1';
                        #goto POS_ALMOST2da;
                    }
                }
                if ($HAP[2]/($HAP[2]+$HAP[1]) > 0.2 && ($total_tools eq "2" || ($type ne "INV" && $type ne "BND")))
                {
                    $SVs{$chr}{$list[1]} = $converted_line;
                    $count{$chr}{$list[1]} = '1';
                    $cutesv{$chr}{$list[1]} = $line;
                    if ($type ne "INV" && $type ne "BND" && $type ne "DUP" && $HAP[0] ne "" && $HAP[0] ne "./." && $prev_match2 eq "")
                    {
                        $count{$chr}{$list[1]} = '2';
                    }
                }         
            }
        }
    }
}

undef %pbsv2;
undef %cutesv2;

if ($input_sniffles ne "")
{
SNIFFLES: while (my $line = <INPUT_SNIFFLES>)
    {
        chomp($line);
        my $first = substr $line, 0, 1;
        if ($first ne "#")
        {
            my @list = split /\t/, $line;
            my $pos = $list[1];
            my $length;
            my $type;
            my $v = '1';
            my $min = '1';
            
            my $info = $list[7];      
            my @info = split /;/, $info;
            my @HAP = split /:/, $list[9];
            
            foreach my $info_tmp (@info)
            {
                my $first_five = substr $info_tmp, 0, 5;
                if ($first_five eq "SVLEN")
                {
                    $length = substr $info_tmp, 6;
                }
                elsif ($first_five eq "SVTYP")
                {
                    $type = substr $info_tmp, 7;
                }
            }
            if ($length eq "")
            {
                $length = ".";
            }
            if ($length < 0)
            {
                $length *= -1;
            }
            my $chr = $list[0];
            if ($list[0] =~ m/chr(\d+|X|Y)/)
            {
                $chr = $1;
            }
            my $converted_line = $chr."\t".$list[1]."\t".$length."\t".$type."\t".$HAP[0];
            
            if (exists ($sniffles{$chr}{$list[1]}))
            {
                next SNIFFLES;
            }
            
            if ($list[6] eq "PASS" && ($length eq "." || $length > 45) && ($HAP[2] >= $min_coverage || ($HAP[2] >= $min_coverage-1 && $HAP[2]/($HAP[2]+$HAP[1]) > 0.3))
                && $HAP[2]/($HAP[2]+$HAP[1]) > 0.1 && ($HAP[0] ne "0/0" || $input_cutesv  eq ""))
            {
                my $prev_match_pbsv = "";
                my $prev_match_cutesv = "";
                if (exists($pbsv2{$chr}{$list[1]}))
                {
                    $prev_match_pbsv = "yes";
                }
                if (exists($cutesv2{$chr}{$list[1]}))
                {
                    $prev_match_cutesv = "yes";
                }
                
                if ((exists($pbsv{$chr}{$list[1]}) && $prev_match_pbsv eq "") || (exists($cutesv{$chr}{$list[1]}) && $prev_match_cutesv eq ""))
                {
                    my @list2 = split /\t/, $SVs{$chr}{$list[1]};
                    
                    if (exists($pbsv{$chr}{$list[1]}))
                    {
                        $pbsv2{$chr}{$list[1]} = undef;
                    }
                    if (exists($cutesv{$chr}{$list[1]}))
                    {
                        $cutesv2{$chr}{$list[1]} = undef;
                    }
                    
                    my $no_pbsv = "";
                    if (exists($pbsv{$chr}{$list[1]}))
                    {                             
                    }
                    else
                    {
                        $no_pbsv = "yes";
                    }
                    
                    my $count_tmp = $count{$chr}{$list[1]};        
                    
                    my $new_pos = $pos;
                    my $new_length = $list2[2];
                    my $new_type = $list2[3];
                    my $new_haplo = $list2[4];
    
                    if (($list2[3] eq $type || $list2[3] eq "BND" || $list2[3] eq "." || ($type eq "." && $list2[3] ne "BND") || ($list2[3] ne "DEL" && $type ne "DEL" && $list2[3] ne "BND")))
                    {
                        $count{$chr}{$list[1]} = $count_tmp+1;
                        
                        if ($no_pbsv ne "")
                        {
                            $new_length = $length;
                        }
                        
                        if ($type eq "INS" && $list2[3] ne $type && $type ne ".")
                        {
                            $new_type = $type; 
                        }
                        if ($list2[3] eq "BND" || ($list2[3] eq "INV" && ($type eq "DEL" || $type eq "INS")) || ($list2[3] eq "DUP" && $type eq "INS") || $list2[3] eq ".")
                        {                 
                            if ($length ne ".")
                            {
                                $new_length = $length;
                            }
                            if ($type ne ".")
                            {
                                $new_type = $type;
                            }
                        }
                        if ($type eq "INV" || $list2[2] eq ".")
                        {                                                                    
                            if ($length ne ".")
                            {
                                $new_length = $length;
                            }                                         
                        }
                        
                        delete $SVs{$chr}{$pos};
                        delete $count{$chr}{$pos};
                        
                        my $line2 = $list2[0]."\t".$new_pos."\t".$new_length."\t".$new_type."\t".$new_haplo;
                        $SVs{$list2[0]}{$new_pos} = $line2;
                        $count{$chr}{$new_pos} = $count_tmp+1;
                        $sniffles{$chr}{$new_pos} = $line;
                        
                        if ($type eq "INS")
                        {
                            $count{$chr}{$new_pos} = $count_tools;
                        }
                        next SNIFFLES;
                    }    
                }
                elsif ($input_pbsv ne "" || $input_cutesv ne "")
                {
                    my $score_prev = "";
                    my $pos_prev = "";
POS_ALMOST2ca:                    
                    while ($v < $length_match)
                    {
POS_ALMOST2c:           my $pos_tmp = ($min*$v)+$pos;

                        my $prev_match_pbsv2 = "";
                        my $prev_match_cutesv2 = "";
                        if (exists($pbsv2{$chr}{$pos_tmp}))
                        {
                            $prev_match_pbsv2 = "yes";
                        }
                        if (exists($cutesv2{$chr}{$pos_tmp}))
                        {
                            $prev_match_cutesv2 = "yes";
                        }
                        
                        if ((exists($pbsv{$chr}{$pos_tmp}) && $prev_match_pbsv2 eq "") || (exists($cutesv{$chr}{$pos_tmp}) && $prev_match_cutesv2 eq ""))
                        {
                            my @list2 = split /\t/, $SVs{$chr}{$pos_tmp};
                            my $count_tmp = $count{$chr}{$pos_tmp};

                            my $no_pbsv = "";
                            if (exists($pbsv{$chr}{$pos_tmp}))
                            {                             
                            }
                            else
                            {
                                $no_pbsv = "yes";
                            }

                            if (exists ($sniffles{$chr}{$pos_tmp}) && $list2[3] eq $type && ($list2[2] > $length*$lenght_margin1 && $list2[2] < $length*$lenght_margin2))
                            {
                                next SNIFFLES;
                            }
                            else
                            {                  
                                if ($list2[3] eq $type || $list2[3] eq "BND" || $list2[3] eq "." || ($type eq "." && $list2[3] ne "BND") ||
                                      ($list2[3] ne "DEL" && $type ne "DEL" && $list2[3] ne "BND"))
                                {
                                    my $score = '0';
                                    if ($list2[3] eq $type || ($type eq "INS" && $list2[3] eq "DUP") || ($type eq "DUP" && $list2[3] eq "INS"))
                                    {
                                        $score += 1;
                                    }
                                    my $g = (800-$v)/800;
                                    $score += $g;
                                    my $score_tmp = '0';
                                    if ($list2[2] ne "." && $length ne ".")
                                    {
                                        if ($list2[2] >= $length)
                                        {
                                            $score_tmp = 1-(($list2[2]-$length)/$length);                                           
                                            if ($score_tmp < 0)
                                            {
                                                $score_tmp = '0';
                                            }
                                        }
                                        else
                                        {
                                            $score_tmp = 1-(($length-$list2[2])/$length);
                                            if ($score_tmp < 0)
                                            {
                                                $score_tmp = '0';
                                            }
                                        }
                                    }
                                    $score += $score_tmp;
        
                                    if ($score_prev ne "" && $score_prev > $score && $score_prev ne "no")
                                    {
                                        $pos_tmp = $pos_prev;
                                        $count_tmp = $count{$chr}{$pos_tmp};
                                        @list2 = split /\t/, $SVs{$chr}{$pos_tmp};
                                    }
        
                                    if ($score > 2 || ($score_prev ne "" && $score_prev > $score) || $score_prev eq "no")
                                    {
                                        
                                        if (exists($pbsv{$chr}{$pos_tmp}))
                                        {
                                            $pbsv2{$chr}{$pos_tmp} = undef;
                                        }
                                        if (exists($cutesv{$chr}{$pos_tmp}))
                                        {
                                            $cutesv2{$chr}{$pos_tmp} = undef;
                                        }
                                        
                                        my $new_pos = $pos_tmp;
                                        my $new_length = $list2[2];
                                        my $new_type = $list2[3];
                                        my $new_haplo = $list2[4];
                                        
                                        if ($no_pbsv eq "yes")
                                        {
                                            $new_pos = $pos;
                                            $new_length = $length;
                                        }
                                        
                                        if ($type eq "INS" && $list2[3] ne $type && $type ne ".")
                                        {                              
                                            $new_type = $type;                                   
                                        }
                                        if ($list2[3] eq "BND" || ($list2[3] eq "INV" && ($type eq "DEL" || $type eq "INS")) || $list2[3] eq "." || ($type eq "INS" && $list2[3] ne $type))
                                        {        
                                            if ($length ne ".")
                                            {
                                                $new_length = $length;
                                            }
                                            if ($type ne ".")
                                            {
                                                $new_type = $type;
                                            }
                                        }
                                        if ($type eq "INV" || $list2[2] eq ".")
                                        {                                                                    
                                            if ($length ne ".")
                                            {
                                                $new_length = $length;
                                            }
                                            $new_pos = $pos;
                                        }
            
                                        delete $SVs{$chr}{$pos_tmp};
                                        delete $count{$chr}{$pos_tmp};
                                        my $line2 = $list2[0]."\t".$new_pos."\t".$new_length."\t".$new_type."\t".$new_haplo;
                                        $SVs{$list2[0]}{$new_pos} = $line2;
                                        $count{$chr}{$new_pos} = $count_tmp+1;
                                        $sniffles{$chr}{$pos_tmp} = $line;
                                        
                                        if ($type eq "INS")
                                        {
                                            $count{$chr}{$new_pos} = $count_tools;
                                        }
                                        if ($new_pos eq $pos)
                                        {
                                            if (exists($pbsv{$chr}{$pos_tmp}))
                                            {
                                                my $tmp = $pbsv{$chr}{$pos_tmp};
                                                delete $pbsv{$chr}{$pos_tmp};
                                                $pbsv{$chr}{$new_pos} = $tmp;
                                            }
                                            if (exists($cutesv{$chr}{$pos_tmp}))
                                            {
                                                my $tmp = $cutesv{$chr}{$pos_tmp};
                                                delete $cutesv{$chr}{$pos_tmp};
                                                $cutesv{$chr}{$new_pos} = $tmp;
                                            }
                                            if (exists($sniffles{$chr}{$pos_tmp}))
                                            {
                                                my $tmp = $sniffles{$chr}{$pos_tmp};
                                                delete $sniffles{$chr}{$pos_tmp};
                                                $sniffles{$chr}{$new_pos} = $tmp;
                                            }
                                        }
                                        next SNIFFLES;
                                    }
                                    elsif ($min eq '1')
                                    {                   
                                        $min = '-1';
                                        if ($score > $score_prev)
                                        {
                                            $score_prev = $score;
                                            $pos_prev = $pos_tmp;
                                        }
                                        goto POS_ALMOST2ca;
                                    }
                                    else
                                    {
                                        $min = '1';
                                        if ($score > $score_prev)
                                        {
                                            $score_prev = $score;
                                            $pos_prev = $pos_tmp;
                                        }
                                    }
                                    $v++;
                                    goto POS_ALMOST2ca;
                                }
                                next SNIFFLES;
                            }
                        }
                        elsif ($min eq '1')
                        {                   
                            $min = '-1';
                            goto POS_ALMOST2c;
                        }
                        else
                        {
                            $min = '1';
                        }
                        $v++;
                    }
                    if ($score_prev ne "" && $score_prev ne "no" && $score_prev > 1.6)
                    {
                        $score_prev = "no";
                        $v = '1';
                        goto POS_ALMOST2ca;
                    }
                }

                $SVs{$chr}{$list[1]} = $converted_line;
                $count{$chr}{$list[1]} = '1';
                $sniffles{$chr}{$list[1]} = $line;
                if (($type eq "INS" || $type eq "DEL") && $HAP[0] ne "0/0" && $HAP[0] ne "./." && $HAP[0] ne "1/1" && ($HAP[2]+$HAP[1]) > 4)
                {
                    $count{$chr}{$list[1]} = $count_tools;
                }
                elsif ($high_recall eq "2" && $HAP[0] ne "0/0" && $HAP[0] ne "1/1" && ($HAP[2]+$HAP[1]) > 4 && $input_cutesv eq "")
                {
                    $count{$chr}{$list[1]} = '2';
                    $SVs{$chr}{$list[1]} = $converted_line;
                    $sniffles{$chr}{$list[1]} = $line;
                }
                elsif ($high_recall eq "2" && $HAP[0] ne "0/0" && $HAP[0] ne "1/1" && $input_cutesv eq "")
                {
                    $count{$chr}{$list[1]} = '2';
                    $SVs{$chr}{$list[1]} = $converted_line;
                    $sniffles{$chr}{$list[1]} = $line;
                }
            }
        }
    }
}

undef %pbsv2;
undef %cutesv2;

if ($input_nanovar ne "")
{
NANOVAR: while (my $line = <INPUT_NANOVAR>)
    {
        chomp($line);
        my $score_prev = "";
        my $pos_prev = "";
        my $first_nuc = substr $line, 0, 1;
        if ($first_nuc ne "#")
        {
            my @list = split /\t/, $line;
            my $pos = $list[1];
            my $type;
            my $length;
            my $v = '1';
            my $min = '1';
            
            my $info = $list[7];      
            my @info = split /;/, $info;
            my @HAP = split /:/, $list[9];
            my @HAP_REF = split /,/, $HAP[1];
            my @HAP2 = split /,/, $HAP[2];
            my $not_save = "";
            
            foreach my $info_tmp (@info)
            {
                my $first_five = substr $info_tmp, 0, 5;
                if ($info_tmp =~ m/SVLEN=>*-*(\d+)/)
                {
                    $length = $1;
                }
                elsif ($first_five eq "SVTYP")
                {
                    $type = substr $info_tmp, 7;
                }
            }
            if ($length eq "")
            {
                $length = ".";
            }
            if ($length < 0)
            {
                $length *= -1;
            }
            my $chr = $list[0];
            if ($list[0] =~ m/chr(\d+|X|Y)/)
            {
                $chr = $1;
            }
            my $converted_line = $chr."\t".$list[1]."\t".$length."\t".$type."\t".$HAP[0];
            
            if (exists ($nanovar{$chr}{$list[1]}))
            {
                next NANOVAR;
            }
            
            if (($HAP2[1] >= $min_coverage || ($HAP2[1] >= $min_coverage-1 && $HAP2[0]/($HAP_REF[0]+$HAP2[0]) > 0.3)) && $list[6] eq "PASS" && $type ne "."
                && $type ne "BND" && ($HAP2[0]/($HAP_REF[0]+$HAP2[0]) > 0.12 || $HAP_REF[0] eq ".") && ($length eq "." || $length > 45) && ($HAP[0] ne "./." || $count_tools < 4))
            {    
                my @list2 = split /\t/, $SVs{$chr}{$list[1]};
                
                my $prev_match_pbsv = "";
                my $prev_match_cutesv = "";
                my $prev_match_sniffles = "";
                if (exists($pbsv2{$chr}{$list[1]}))
                {
                    $prev_match_pbsv = "yes";
                }
                if (exists($cutesv2{$chr}{$list[1]}))
                {
                    $prev_match_cutesv = "yes";
                }
                if (exists($sniffles2{$chr}{$list[1]}))
                {
                    $prev_match_sniffles = "yes";
                }
                
                if ((exists($pbsv{$chr}{$list[1]}) && $prev_match_pbsv eq "") || (exists($cutesv{$chr}{$list[1]}) && $prev_match_cutesv eq "") ||
                    (exists($sniffles{$chr}{$list[1]}) && $prev_match_sniffles eq ""))
                {
                    if (exists($pbsv{$chr}{$list[1]}))
                    {
                        $pbsv2{$chr}{$list[1]} = undef;
                    }
                    if (exists($cutesv{$chr}{$list[1]}))
                    {
                        $cutesv2{$chr}{$list[1]} = undef;
                    }
                    if (exists($sniffles{$chr}{$list[1]}))
                    {
                        $sniffles2{$chr}{$list[1]} = undef;
                    }
                                     
                    my $count_tmp = $count{$chr}{$list[1]};

                    if (($list2[3] eq $type || $list2[3] eq "." || $type eq "." || ($list2[3] ne "DEL" && $type ne "DEL") || $type eq "INV" || $list2[3] eq "INV")
                        && ($list2[3] ne "BND" || exists($sniffles{$chr}{$list[1]}) || exists($cutesv{$chr}{$list[1]})))
                    {
                        $count{$chr}{$list[1]} = $count_tmp+1;
                        
                        my $new_pos = $pos;
                        my $new_length = $list2[2];
                        my $new_type = $list2[3];
                        my $new_haplo = $list2[4];
                        
                        if (($list2[3] eq "BND" || ($list2[3] eq "INV" && ($type eq "DEL" || $type eq "INS")) || $list2[3] eq ".") && $type ne "DUP")
                        {                 
                            if ($type ne ".")
                            {
                                $new_type = $type;
                            }   
                        }
                        if ($type eq "INV")
                        {
                            if ($length ne ".")
                            {
                                $new_length = $length;
                            }
                            if ($type ne ".")
                            {
                                $new_type = $type;
                            }
                            if ($HAP[0] ne "./.")
                            {
                                $new_haplo = $HAP[0];
                            }
                        }
                        if ($list2[2] eq ".")
                        {
                            $new_length = $length;
                        }
                        if (exists($pbsv{$chr}{$pos}))
                        {}
                        else
                        {
                            if ($HAP[0] ne "./.")
                            {
                                $new_haplo = $HAP[0];
                            }
                        }
                        if ($new_haplo eq "./.")
                        {                             
                            $new_haplo = $HAP[0];                           
                        }  

                        delete $SVs{$chr}{$pos};
                        delete $count{$chr}{$pos};
                        my $line2 = $list2[0]."\t".$new_pos."\t".$new_length."\t".$new_type."\t".$new_haplo;
                        $SVs{$list2[0]}{$new_pos} = $line2;
                        $count{$chr}{$new_pos} = $count_tmp+1;
                        $nanovar{$chr}{$new_pos} = $line;      
                        
                        if (($HAP[0] eq "0/1" || $HAP[0] eq "1/0") && $type eq "INS" && $length < 2000)
                        {
                            $count{$chr}{$new_pos} = $count_tmp+2;
                        }
                    }
                    next NANOVAR;
                }
                else
                {
                    
POS_ALMOST2fa:                    
                    while ($v < $length_match)
                    {
POS_ALMOST2f:           my $pos_tmp = ($min*$v)+$pos;

                        my $prev_match_pbsv2 = "";
                        my $prev_match_cutesv2 = "";
                        my $prev_match_sniffles2 = "";
                        if (exists($pbsv2{$chr}{$pos_tmp}))
                        {
                            $prev_match_pbsv2 = "yes";
                        }
                        if (exists($cutesv2{$chr}{$pos_tmp}))
                        {
                            $prev_match_cutesv2 = "yes";
                        }
                        if (exists($sniffles2{$chr}{$pos_tmp}))
                        {
                            $prev_match_sniffles2 = "yes";
                        }
                        
                        if ((exists($pbsv{$chr}{$pos_tmp}) && $prev_match_pbsv2 eq "") || (exists($cutesv{$chr}{$pos_tmp}) && $prev_match_cutesv2 eq "") ||
                            (exists($sniffles{$chr}{$pos_tmp}) && $prev_match_sniffles2 eq ""))
                        {
                            my $count_tmp = $count{$chr}{$pos_tmp};
                            if (exists($nanovar{$chr}{$pos_tmp}) && $high_recall eq "gg")
                            {
                                if ($min eq '1')
                                {                   
                                    $min = '-1';
                                    goto POS_ALMOST2f;
                                }
                                else
                                {
                                    $min = '1';
                                }
                            }
                            else
                            {                              
                                my @list2 = split /\t/, $SVs{$chr}{$pos_tmp};
                                
                                if (($list2[3] eq $type || ($list2[3] eq "BND" && $type eq "INV") || $list2[3] eq "." || $type eq "." || ($list2[3] ne "DEL" && $type ne "DEL")
                                     || $type eq "INV" || $list2[3] eq "INV") && ($list2[3] ne "BND" || exists($sniffles{$chr}{$pos_tmp}) || exists($cutesv{$chr}{$pos_tmp})))             
                                {         
                                    $count{$chr}{$pos_tmp} = $count_tmp+1;
                                    
                                    my $score = '0';
                                    if ($list2[3] eq $type || ($type eq "INS" && $list2[3] eq "DUP") || ($type eq "DUP" && $list2[3] eq "INS"))
                                    {
                                        $score += 1;
                                    }
                                    my $g = (800-$v)/800;
                                    $score += $g;
                                    my $score_tmp = '0';
                                    if ($list2[2] ne "." && $length ne ".")
                                    {
                                        if ($list2[2] >= $length)
                                        {
                                            $score_tmp = 1-(($list2[2]-$length)/$length);                                           
                                            if ($score_tmp < 0)
                                            {
                                                $score_tmp = '0';
                                            }
                                        }
                                        else
                                        {
                                            $score_tmp = 1-(($length-$list2[2])/$length);
                                            if ($score_tmp < 0)
                                            {
                                                $score_tmp = '0';
                                            }
                                        }
                                    }
                                    $score += $score_tmp;
        
                                    if ($score_prev ne "" && $score_prev > $score && $score_prev ne "no")
                                    {
                                        $pos_tmp = $pos_prev;
                                        $count_tmp = $count{$chr}{$pos_tmp};
                                        @list2 = split /\t/, $SVs{$chr}{$pos_tmp};
                                    }
        
                                    if ($score > 2 || ($score_prev ne "" && $score_prev > $score) || $score_prev eq "no")
                                    {                                                                 
                                        if (exists($pbsv{$chr}{$pos_tmp}))
                                        {
                                            $pbsv2{$chr}{$pos_tmp} = undef;
                                        }
                                        if (exists($cutesv{$chr}{$pos_tmp}))
                                        {
                                            $cutesv2{$chr}{$pos_tmp} = undef;
                                        }
                                        if (exists($sniffles{$chr}{$pos_tmp}))
                                        {
                                            $sniffles2{$chr}{$pos_tmp} = undef;
                                        }
                                        
                                        my $new_pos = $pos_tmp;
                                        my $new_length = $list2[2];
                                        my $new_type = $list2[3];
                                        my $new_haplo = $list2[4];               
                                        
                                        if (($list2[3] eq "BND" || $list2[3] eq ".") && $type ne "DUP")
                                        {                 
                                            if ($length ne ".")
                                            {
                                                $new_length = $length;
                                            }
                                            if ($type ne ".")
                                            {
                                                $new_type = $type;
                                            }                            
                                        }
                                        if ($type eq "INV")
                                        {
                                            $new_pos = $pos;
                                            
                                            if (exists($pbsv{$chr}{$pos_tmp}))
                                            {
                                                my $tmp = $pbsv{$chr}{$pos_tmp};
                                                delete $pbsv{$chr}{$pos_tmp};
                                                $pbsv{$chr}{$new_pos} = $tmp;
                                            }
                                            if (exists($cutesv{$chr}{$pos_tmp}))
                                            {
                                                my $tmp = $cutesv{$chr}{$pos_tmp};
                                                delete $cutesv{$chr}{$pos_tmp};
                                                $cutesv{$chr}{$new_pos} = $tmp;
                                            }
                                            if (exists($sniffles{$chr}{$pos_tmp}))
                                            {
                                                my $tmp = $sniffles{$chr}{$pos_tmp};
                                                delete $sniffles{$chr}{$pos_tmp};
                                                $sniffles{$chr}{$new_pos} = $tmp;
                                            }
                                            if (exists($nanovar{$chr}{$pos_tmp}))
                                            {
                                                my $tmp = $nanovar{$chr}{$pos_tmp};
                                                delete $nanovar{$chr}{$pos_tmp};
                                                $nanovar{$chr}{$new_pos} = $tmp;
                                            }
                                            
                                            if ($length ne ".")
                                            {
                                                $new_length = $length;
                                            }
                                            if ($type ne ".")
                                            {
                                                $new_type = $type;
                                            }
                                            if ($HAP[0] ne "./.")
                                            {
                                                $new_haplo = $HAP[0];
                                            }
                                        }
                                        if ($list2[2] eq ".")
                                        {
                                            $new_length = $length;
                                        }
                                        if (exists($pbsv{$chr}{$pos_tmp}))
                                        {}
                                        else
                                        {
                                            if ($HAP[0] ne "./.")
                                            {
                                                $new_haplo = $HAP[0];
                                            }
                                        }
                                        if ($new_haplo eq "./.")
                                        {                             
                                            $new_haplo = $HAP[0];                           
                                        }  
            
                                        delete $SVs{$chr}{$pos_tmp};
                                        delete $count{$chr}{$pos_tmp};
                                        my $line2 = $list2[0]."\t".$new_pos."\t".$new_length."\t".$new_type."\t".$new_haplo;
                                        $SVs{$list2[0]}{$new_pos} = $line2;
                                        $count{$chr}{$new_pos} = $count_tmp+1;
                                        $nanovar{$chr}{$new_pos} = $line;
                                        
                                        if (($HAP[0] eq "0/1" || $HAP[0] eq "1/0") && $type eq "INS" && $length < 2000)
                                        {
                                            $count{$chr}{$new_pos} = $count_tmp+2;
                                        }
                                        next NANOVAR;
                                    }
                                    elsif ($min eq '1')
                                    {                   
                                        $min = '-1';
                                        if ($score > $score_prev)
                                        {
                                            $score_prev = $score;
                                            $pos_prev = $pos_tmp;
                                        }
                                        goto POS_ALMOST2fa;
                                    }
                                    else
                                    {
                                        $min = '1';
                                        if ($score > $score_prev)
                                        {
                                            $score_prev = $score;
                                            $pos_prev = $pos_tmp;
                                        }
                                    }
                                    $v++;
                                    goto POS_ALMOST2fa;
                                }
                            }
                            next NANOVAR;
                        }
                        elsif ($min eq '1')
                        {                   
                            $min = '-1';
                            goto POS_ALMOST2f;
                        }
                        else
                        {
                            $min = '1';
                        }
                        $v++;
                    }
                    
                }               
        
                if ($type ne "BND" && $HAP[0] ne "./.")
                {
                    $SVs{$chr}{$list[1]} = $converted_line;
                    $nanovar{$chr}{$list[1]} = $line;
                    $count{$chr}{$list[1]} = '1';
                    if (($HAP[0] eq "0/1" || $HAP[0] eq "1/0") && $type eq "INS" && $pos_prev eq "")
                    {
                        $count{$chr}{$list[1]} = 2;
                    }
                }                             
            }
        }
    }
}

undef %pbsv2;
undef %cutesv2;
undef %sniffles2;

my %svim;
if ($input_svim ne "")
{
SVIM: while (my $line = <INPUT_SVIM>)
    {
        chomp($line);
        my $first_nuc = substr $line, 0, 1;
        if ($first_nuc ne "#")
        {
            my @list = split /\t/, $line;
            my $pos = $list[1];
            my $length;
            my $type;
            my $v = '1';
            my $min = '1';
            
            my $info = $list[7];      
            my @info = split /;/, $info;
            my @HAP = split /:/, $list[9];
            
            my @HAP_INFO = split /:/, $list[8];
            my $support = '0';
        
            foreach my $info_tmp (@info)
            {
                my $first_five = substr $info_tmp, 0, 5;
                if ($first_five eq "SVLEN")
                {
                    $length = substr $info_tmp, 6;
                }
                elsif ($first_five eq "SVTYP")
                {
                    $type = substr $info_tmp, 7;
                    if ($type eq "DUP:TANDEM" || $type eq "DUP:INT")
                    {
                        $type = "DUP";
                    }
                }
                elsif ($first_five eq "SUPPO")
                {
                    $support = substr $info_tmp, 8;
                }
            }
            if ($length eq "")
            {
                $length = ".";
            }
            if ($length < 0)
            {
                $length *= -1;
            }
            my $chr = $list[0];
            if ($list[0] =~ m/chr(\d+|X|Y)/)
            {
                $chr = $1;
            }
            my $hapi = $HAP[0];
            
            if (exists ($svim{$chr}{$list[1]}))
            {
                next SVIM;
            }
            my $ratio = '0';
            
            if ($HAP[1] > 0)
            {
                $ratio = $support/($HAP[1])
            }
            
            if ($list[6] eq "PASS" && ($support >= $min_coverage || ($support >= $min_coverage-1 && $ratio > 0.3)) && $type ne "BND" && $type ne "DUP" && $hapi ne "0/0" && $hapi ne "./." &&
                ($ratio > 0.2 || $HAP_INFO[1] eq "CN" || $HAP[1] eq ".") && ($length eq "." || $length > 45))
            {          
                my $converted_line = $chr."\t".$list[1]."\t".$length."\t".$type."\t".$hapi; 
                my $prev_match_pbsv = "";
                my $prev_match_cutesv = "";
                my $prev_match_sniffles = "";
                my $prev_match_nanovar = "";
                
                if (exists($pbsv2{$chr}{$list[1]}))
                {
                    $prev_match_pbsv = "yes";
                }
                if (exists($cutesv2{$chr}{$list[1]}))
                {
                    $prev_match_cutesv = "yes";
                }
                if (exists($sniffles2{$chr}{$list[1]}))
                {
                    $prev_match_sniffles = "yes";
                }
                if (exists($nanovar2{$chr}{$list[1]}))
                {
                    $prev_match_nanovar = "yes";
                }
                
                if ((exists($pbsv{$chr}{$list[1]}) && $prev_match_pbsv eq "") || (exists($cutesv{$chr}{$list[1]}) && $prev_match_cutesv eq "") ||
                    (exists($sniffles{$chr}{$list[1]}) && $prev_match_sniffles eq "") || (exists($nanovar{$chr}{$list[1]}) && $prev_match_nanovar eq ""))
                {
                    if (exists($pbsv{$chr}{$list[1]}))
                    {
                        $pbsv2{$chr}{$list[1]} = undef;
                    }
                    if (exists($cutesv{$chr}{$list[1]}))
                    {
                        $cutesv2{$chr}{$list[1]} = undef;
                    }
                    if (exists($sniffles{$chr}{$list[1]}))
                    {
                        $sniffles2{$chr}{$list[1]} = undef;
                    }
                    if (exists($nanovar{$chr}{$list[1]}))
                    {
                        $nanovar2{$chr}{$list[1]} = undef;
                    }
                    
                    my $count_tmp = $count{$chr}{$list[1]};
                    my @list2 = split /\t/, $SVs{$chr}{$list[1]};
                    
                    my $new_pos = $pos;
                    my $new_length = $list2[2];
                    my $new_type = $list2[3];
                    my $new_haplo = $list2[4];
                    
                    if ($new_haplo eq "./.")
                    {                             
                        $new_haplo = $hapi;                           
                    }  
                                 
                    delete $SVs{$chr}{$pos};
                    delete $count{$chr}{$pos};
                    my $line2 = $list2[0]."\t".$new_pos."\t".$new_length."\t".$new_type."\t".$new_haplo;
                    $SVs{$list2[0]}{$new_pos} = $line2;
                    $count{$chr}{$new_pos} = $count_tmp+1;
                    $svim{$chr}{$new_pos} = $line;
                    
                    if ($hapi eq "1/1")
                    {
                        $count{$chr}{$new_pos} = $count_tools;
                    }
            
                    next SVIM;
                }               
                else
                {
                    my $score_prev = "";
                    my $pos_prev = "";
POS_ALMOST2ga:                   
                    while ($v < $length_match)
                    {
POS_ALMOST2g:           my $pos_tmp = ($min*$v)+$pos;
                      
                        my $prev_match_pbsv2 = "";
                        my $prev_match_cutesv2 = "";
                        my $prev_match_sniffles2 = "";
                        my $prev_match_nanovar2 = "";
                        
                        if (exists($pbsv2{$chr}{$pos_tmp}))
                        {
                            $prev_match_pbsv2 = "yes";
                        }
                        if (exists($cutesv2{$chr}{$pos_tmp}))
                        {
                            $prev_match_cutesv2 = "yes";
                        }
                        if (exists($sniffles2{$chr}{$pos_tmp}))
                        {
                            $prev_match_sniffles2 = "yes";
                        }
                        if (exists($nanovar2{$chr}{$pos_tmp}))
                        {
                            $prev_match_nanovar2 = "yes";
                        }

                        if ((exists($pbsv{$chr}{$pos_tmp}) && $prev_match_pbsv2 eq "") || (exists($cutesv{$chr}{$pos_tmp}) && $prev_match_cutesv2 eq "") ||
                            (exists($sniffles{$chr}{$pos_tmp}) && $prev_match_sniffles2 eq "") || (exists($nanovar{$chr}{$pos_tmp}) && $prev_match_nanovar2 eq ""))
                        {
                            my $count_tmp = $count{$chr}{$pos_tmp};
                            my @list2 = split /\t/, $SVs{$chr}{$pos_tmp};
                            
                            if ($count_tmp > 20)
                            {
                                if ($min eq '1')
                                {                   
                                    $min = '-1';
                                    goto POS_ALMOST2g;
                                }
                                else
                                {
                                    $min = '1';
                                }
                            }
                            elsif (exists($svim{$chr}{$pos_tmp}) && $list2[3] eq $type && ($list2[2] > $length*$lenght_margin1 && $list2[2] < $length*$lenght_margin2))
                            {
                                next SVIM;
                            }
                            elsif ($list2[3] ne "DUP")
                            {                       
                                my $score = '0';
                                if ($list2[3] eq $type || ($type eq "INS" && $list2[3] eq "DUP") || ($type eq "DUP" && $list2[3] eq "INS"))
                                {
                                    $score += 1;
                                }
                                my $g = (800-$v)/800;
                                $score += $g;
                                my $score_tmp = '0';
                                if ($list2[2] ne "." && $length ne ".")
                                {
                                    if ($list2[2] >= $length)
                                    {
                                        $score_tmp = 1-(($list2[2]-$length)/$length);                                           
                                        if ($score_tmp < 0)
                                        {
                                            $score_tmp = '0';
                                        }
                                    }
                                    else
                                    {
                                        $score_tmp = 1-(($length-$list2[2])/$length);
                                        if ($score_tmp < 0)
                                        {
                                            $score_tmp = '0';
                                        }
                                    }
                                }
                                $score += $score_tmp;
    
                                if ($score_prev ne "" && $score_prev > $score && $score_prev ne "no")
                                {
                                    $pos_tmp = $pos_prev;
                                    $count_tmp = $count{$chr}{$pos_tmp};
                                    @list2 = split /\t/, $SVs{$chr}{$pos_tmp};
                                }
    
                                if ($score > 2 || ($score_prev ne "" && $score_prev > $score) || $score_prev eq "no")
                                {                                
                                    if (exists($pbsv{$chr}{$pos_tmp}))
                                    {
                                        $pbsv2{$chr}{$pos_tmp} = undef;
                                    }
                                    if (exists($cutesv{$chr}{$pos_tmp}))
                                    {
                                        $cutesv2{$chr}{$pos_tmp} = undef;
                                    }
                                    if (exists($sniffles{$chr}{$pos_tmp}))
                                    {
                                        $sniffles2{$chr}{$pos_tmp} = undef;
                                    }
                                    if (exists($nanovar{$chr}{$pos_tmp}))
                                    {
                                        $nanovar2{$chr}{$pos_tmp} = undef;
                                    }
                            
                                    my $new_pos = $pos_tmp;
                                    my $new_length = $list2[2];
                                    my $new_type = $list2[3];
                                    my $new_haplo = $list2[4];
                                    
                                    if (exists($sniffles{$chr}{$pos_tmp}))
                                    {}
                                    elsif (exists($cutesv{$chr}{$pos_tmp}))
                                    {}
                                    elsif (exists($pbsv{$chr}{$pos_tmp}))
                                    {}
                                    elsif ($type eq "INV")
                                    {
                                        $new_pos = $pos;
                                    }
                                    if ($new_haplo eq "./.")
                                    {                             
                                        $new_haplo = $hapi;                           
                                    }
                                                    
                                    delete $SVs{$chr}{$pos_tmp};
                                    delete $count{$chr}{$pos_tmp};
                                    my $line2 = $list2[0]."\t".$new_pos."\t".$new_length."\t".$new_type."\t".$new_haplo;
                                    $SVs{$list2[0]}{$new_pos} = $line2;
                                    $count{$chr}{$new_pos} = $count_tmp+1;
                                    $svim{$chr}{$new_pos} = $line;
                                    
                                    if ($hapi eq "1/1")
                                    {
                                        $count{$chr}{$new_pos} = $count_tools;
                                    }
                                    
                                    if ($new_pos eq $pos)
                                    {
                                        if (exists($pbsv{$chr}{$pos_tmp}))
                                        {
                                            my $tmp = $pbsv{$chr}{$pos_tmp};
                                            delete $pbsv{$chr}{$pos_tmp};
                                            $pbsv{$chr}{$new_pos} = $tmp;
                                        }
                                        if (exists($cutesv{$chr}{$pos_tmp}))
                                        {
                                            my $tmp = $cutesv{$chr}{$pos_tmp};
                                            delete $cutesv{$chr}{$pos_tmp};
                                            $cutesv{$chr}{$new_pos} = $tmp;
                                        }
                                        if (exists($sniffles{$chr}{$pos_tmp}))
                                        {
                                            my $tmp = $sniffles{$chr}{$pos_tmp};
                                            delete $sniffles{$chr}{$pos_tmp};
                                            $sniffles{$chr}{$new_pos} = $tmp;
                                        }
                                        if (exists($nanovar{$chr}{$pos_tmp}))
                                        {
                                            my $tmp = $nanovar{$chr}{$pos_tmp};
                                            delete $nanovar{$chr}{$pos_tmp};
                                            $nanovar{$chr}{$new_pos} = $tmp;
                                        }
                                        if (exists($svim{$chr}{$pos_tmp}))
                                        {
                                            my $tmp = $svim{$chr}{$pos_tmp};
                                            delete $svim{$chr}{$pos_tmp};
                                            $svim{$chr}{$new_pos} = $tmp;
                                        }
                                    }                                 
                                    next SVIM;
                                }
                                elsif ($min eq '1')
                                {                   
                                    $min = '-1';
                                    if ($score > $score_prev)
                                    {
                                        $score_prev = $score;
                                        $pos_prev = $pos_tmp;
                                    }
                                    goto POS_ALMOST2ga;
                                }
                                else
                                {
                                    $min = '1';
                                    if ($score > $score_prev)
                                    {
                                        $score_prev = $score;
                                        $pos_prev = $pos_tmp;
                                    }
                                }
                                $v++;
                                goto POS_ALMOST2ga;
                            }          
                        }
                        elsif ($min eq '1')
                        {                   
                            $min = '-1';
                            goto POS_ALMOST2g;
                        }
                        else
                        {
                            $min = '1';
                        }
                        $v++;
                    }
                }

                $svim{$chr}{$list[1]} = $line;
                $SVs{$chr}{$list[1]} = $converted_line;
                if ($hapi eq "1/1")
                {
                    $count{$chr}{$list[1]} = $count_tools;
                }
                else
                {
                    $count{$chr}{$list[1]} = '1';
                }
            }
        }
    }
}
undef %pbsv2;
undef %cutesv2;
undef %sniffles2;
undef %nanovar2;

my %nanosv;
my $nanosv_converted = "";
if ($input_nanosv ne "")
{
NANOSV: while (my $line = <INPUT_NANOSV>)
    {
        chomp($line);
        my $first_nuc = substr $line, 0, 1;
        if ($first_nuc ne "#")
        {
            my @list = split /\t/, $line;
            my $pos = $list[1];
            my $length;
            my $type;
            my $v = '1';
            my $min = '1';
            
            my $qual = $list[6];
            my @qual = split /;/, $qual;
            my $info = $list[7];      
            my @info = split /;/, $info;
            my @HAP = split /:/, $list[9];
            my @HAP_REF = split /,/, $HAP[1];
            my @HAP2 = split /,/, $HAP[2];
            
            foreach my $info_tmp (@info)
            {
                my $first_five = substr $info_tmp, 0, 5;
                if ($first_five eq "SVLEN")
                {
                    $length = substr $info_tmp, 6;
                }
                elsif ($first_five eq "SVTYP")
                {
                    $type = substr $info_tmp, 7;
                }
            }
            if ($length eq "")
            {
                $length = ".";
            }
            if ($length < 0)
            {
                $length *= -1;
            }
            my $chr = $list[0];
            if ($list[0] =~ m/chr(\d+|X|Y)/)
            {
                $chr = $1;
            }
            my $hapi = $HAP[0];
            
            if ($nanosv_converted eq "yes")
            {
                $pos = $list[1];
                $length = $list[2];
                $type = $list[3];
                $hapi = $list[4];
            }

            my $converted_line = $chr."\t".$list[1]."\t".$length."\t".$type."\t".$hapi;
            
            my $qual_check = "";
            my $CIPOS = "";
            my $CIEND = "";
            my $mapq = "";
            my $cluster = "";
            foreach my $qual_tmp (@qual)
            {
                if ($qual_tmp eq "LowQual")
                {
                    $qual_check = "no";
                }
                if ($qual_tmp eq "MapQual")
                {
                    $mapq = "yes";
                }
                if ($qual_tmp eq "SVcluster")
                {
                    $cluster = "yes";
                }
                if ($qual_tmp eq "CIPOS")
                {
                    $CIPOS = "yes";
                }
                if ($qual_tmp eq "CIEND")
                {
                    $CIEND = "yes";
                }
            }
            if ($CIEND ne "yes" && $CIPOS ne "yes" && $mapq eq "yes" && $cluster eq "yes")
            {
                $qual_check = "no";
            }
            
            if (($HAP2[1] >= $min_coverage+1 ) && ($HAP2[0]/($HAP_REF[0]+$HAP2[0]) > 0.2 || $HAP_REF[0] eq ".") && $hapi ne "0/0" && ($length eq "." || $length > 45) && $qual_check ne "no")
            {                
                my $prev_match_pbsv = "";
                my $prev_match_cutesv = "";
                my $prev_match_sniffles = "";
                my $prev_match_nanovar = "";
                my $prev_match_svim = "";
                
                if (exists($pbsv2{$chr}{$list[1]}))
                {
                    $prev_match_pbsv = "yes";
                }
                if (exists($cutesv2{$chr}{$list[1]}))
                {
                    $prev_match_cutesv = "yes";
                }
                if (exists($sniffles2{$chr}{$list[1]}))
                {
                    $prev_match_sniffles = "yes";
                }
                if (exists($nanovar2{$chr}{$list[1]}))
                {
                    $prev_match_nanovar = "yes";
                }
                if (exists($svim2{$chr}{$list[1]}))
                {
                    $prev_match_svim = "yes";
                }
                
                if ((exists($pbsv{$chr}{$list[1]}) && $prev_match_pbsv eq "") || (exists($cutesv{$chr}{$list[1]}) && $prev_match_cutesv eq "") ||
                    (exists($sniffles{$chr}{$list[1]}) && $prev_match_sniffles eq "") || (exists($nanovar{$chr}{$list[1]}) && $prev_match_nanovar eq "")
                    || (exists($svim{$chr}{$list[1]}) && $prev_match_svim eq ""))
                   
                {
                    if (exists($pbsv{$chr}{$list[1]}))
                    {
                        $pbsv2{$chr}{$list[1]} = undef;
                    }
                    if (exists($cutesv{$chr}{$list[1]}))
                    {
                        $cutesv2{$chr}{$list[1]} = undef;
                    }
                    if (exists($sniffles{$chr}{$list[1]}))
                    {
                        $sniffles2{$chr}{$list[1]} = undef;
                    }
                    if (exists($nanovar{$chr}{$list[1]}))
                    {
                        $nanovar2{$chr}{$list[1]} = undef;
                    }
                    if (exists($svim{$chr}{$list[1]}))
                    {
                        $svim2{$chr}{$list[1]} = undef;
                    }
                
                    my $count_tmp = $count{$chr}{$list[1]};
                    my @list3 = split /\t/, $SVs{$chr}{$list[1]};
                    
                    if ($list3[3] eq "DEL" && exists($sniffles{$list[0]}{$list[1]}) && $high_recall eq "1")
                    {
                    }
                    else
                    {
                        my @list2 = split /\t/, $SVs{$chr}{$pos};
                        my $new_pos = $pos;
                        my $new_length = $list2[2];
                        my $new_type = $list2[3];
                        my $new_haplo = $list2[4];
                        

                        if (exists($nanovar{$chr}{$pos}))
                        {}
                        elsif ($list2[3] eq "INV" && $length ne ".")
                        {
                            $new_length = $length;
                        }
                        if (exists($pbsv{$chr}{$list[1]}))                                 
                        {}
                        elsif ($hapi ne "./.")
                        {                             
                            $new_haplo = $hapi;                           
                        } 
                                        
                        delete $SVs{$chr}{$pos};
                        delete $count{$chr}{$pos};
                        my $line2 = $list2[0]."\t".$new_pos."\t".$new_length."\t".$new_type."\t".$new_haplo;
                        $SVs{$list2[0]}{$new_pos} = $line2;
                        $count{$chr}{$new_pos} = $count_tmp+1;
                        $nanosv{$chr}{$new_pos} = $line;
                    }
                    next NANOSV;
                }
                else
                {
                    my $score_prev = "";
                    my $pos_prev = "";
POS_ALMOST2ha:                   
                    while ($v < $length_match)
                    {
POS_ALMOST2h:           my $pos_tmp = ($min*$v)+$pos;
                        
                        my $prev_match_pbsv2 = "";
                        my $prev_match_cutesv2 = "";
                        my $prev_match_sniffles2 = "";
                        my $prev_match_nanovar2 = "";
                        my $prev_match_svim2 = "";
                        
                        if (exists($pbsv2{$chr}{$pos_tmp}))
                        {
                            $prev_match_pbsv2 = "yes";
                        }
                        if (exists($cutesv2{$chr}{$pos_tmp}))
                        {
                            $prev_match_cutesv2 = "yes";
                        }
                        if (exists($sniffles2{$chr}{$pos_tmp}))
                        {
                            $prev_match_sniffles2 = "yes";
                        }
                        if (exists($nanovar2{$chr}{$pos_tmp}))
                        {
                            $prev_match_nanovar2 = "yes";
                        }
                        if (exists($svim2{$chr}{$pos_tmp}))
                        {
                            $prev_match_svim2 = "yes";
                        }
                        
                        if ((exists($pbsv{$chr}{$pos_tmp}) && $prev_match_pbsv2 eq "") || (exists($cutesv{$chr}{$pos_tmp}) && $prev_match_cutesv2 eq "") ||
                            (exists($sniffles{$chr}{$pos_tmp}) && $prev_match_sniffles2 eq "") || (exists($nanovar{$chr}{$pos_tmp}) && $prev_match_nanovar2 eq "")
                            || (exists($svim{$chr}{$pos_tmp}) && $prev_match_svim2 eq ""))
                        {
                            my $count_tmp = $count{$chr}{$pos_tmp};
                            my @list3 = split /\t/, $SVs{$chr}{$list[1]};
                            
                            if ($list3[3] eq "DEL" && exists($sniffles{$list[0]}{$pos_tmp}) && $high_recall eq "1")
                            {
                            }
                            else
                            {
                                if ($count_tmp > 20)
                                {
                                    if ($min eq '1')
                                    {                   
                                        $min = '-1';
                                        goto POS_ALMOST2h;
                                    }
                                    else
                                    {
                                        $min = '1';
                                    }
                                }
                                else
                                {                                            
                                    my @list2 = split /\t/, $SVs{$chr}{$pos_tmp};
                                    
                                    my $score = '0';
                                    if ($list2[3] eq $type || ($type eq "INS" && $list2[3] eq "DUP") || ($type eq "DUP" && $list2[3] eq "INS"))
                                    {
                                        $score += 1;
                                    }
                                    my $g = (800-$v)/800;
                                    $score += $g;
                                    my $score_tmp = '0';
                                    if ($list2[2] ne "." && $length ne ".")
                                    {
                                        if ($list2[2] >= $length)
                                        {
                                            $score_tmp = 1-(($list2[2]-$length)/$length);                                           
                                            if ($score_tmp < 0)
                                            {
                                                $score_tmp = '0';
                                            }
                                        }
                                        else
                                        {
                                            $score_tmp = 1-(($length-$list2[2])/$length);
                                            if ($score_tmp < 0)
                                            {
                                                $score_tmp = '0';
                                            }
                                        }
                                    }
                                    $score += $score_tmp;
        
                                    if ($score_prev ne "" && $score_prev > $score && $score_prev ne "no")
                                    {
                                        $pos_tmp = $pos_prev;
                                        $count_tmp = $count{$chr}{$pos_tmp};
                                        @list2 = split /\t/, $SVs{$chr}{$pos_tmp};
                                    }
        
                                    if ($score > 1.5 || ($score_prev ne "" && $score_prev > $score) || $score_prev eq "no")
                                    {                                                                 
                                        if (exists($pbsv{$chr}{$pos_tmp}))
                                        {
                                            $pbsv2{$chr}{$pos_tmp} = undef;
                                        }
                                        if (exists($cutesv{$chr}{$pos_tmp}))
                                        {
                                            $cutesv2{$chr}{$pos_tmp} = undef;
                                        }
                                        if (exists($sniffles{$chr}{$pos_tmp}))
                                        {
                                            $sniffles2{$chr}{$pos_tmp} = undef;
                                        }
                                        if (exists($nanovar{$chr}{$pos_tmp}))
                                        {
                                            $nanovar2{$chr}{$pos_tmp} = undef;
                                        }
                                        if (exists($svim{$chr}{$pos_tmp}))
                                        {
                                            $svim2{$chr}{$pos_tmp} = undef;
                                        }
                                        
                                        my $new_pos = $pos_tmp;
                                        my $new_length = $list2[2];
                                        my $new_type = $list2[3];
                                        my $new_haplo = $list2[4];
                                        
    
                                        if (exists($nanovar{$chr}{$pos_tmp}))
                                        {}
                                        elsif ($list2[3] eq "INV" && $length ne ".")
                                        {
                                            $new_length = $length;
                                        }
                                        if (exists($pbsv{$chr}{$list[1]}))                                 
                                        {}
                                        elsif ($hapi ne "./.")
                                        {                             
                                            $new_haplo = $hapi;
                                        }
                                        if (exists($pbsv{$chr}{$list[1]}))                                 
                                        {}
                                        else
                                        {
                                            $new_pos = $pos;
                                        }
                                                        
                                        delete $SVs{$chr}{$pos_tmp};
                                        delete $count{$chr}{$pos_tmp};
                                        my $line2 = $list2[0]."\t".$new_pos."\t".$new_length."\t".$new_type."\t".$new_haplo;
                                        $SVs{$list2[0]}{$new_pos} = $line2;
                                        $count{$chr}{$new_pos} = $count_tmp+1;
                                        $nanosv{$chr}{$new_pos} = $line;
                                        
                                        if ($new_pos eq $pos)
                                        {
                                            if (exists($pbsv{$chr}{$pos_tmp}))
                                            {
                                                my $tmp = $pbsv{$chr}{$pos_tmp};
                                                delete $pbsv{$chr}{$pos_tmp};
                                                $pbsv{$chr}{$new_pos} = $tmp;
                                            }
                                            if (exists($cutesv{$chr}{$pos_tmp}))
                                            {
                                                my $tmp = $cutesv{$chr}{$pos_tmp};
                                                delete $cutesv{$chr}{$pos_tmp};
                                                $cutesv{$chr}{$new_pos} = $tmp;
                                            }
                                            if (exists($sniffles{$chr}{$pos_tmp}))
                                            {
                                                my $tmp = $sniffles{$chr}{$pos_tmp};
                                                delete $sniffles{$chr}{$pos_tmp};
                                                $sniffles{$chr}{$new_pos} = $tmp;
                                            }
                                            if (exists($nanovar{$chr}{$pos_tmp}))
                                            {
                                                my $tmp = $nanovar{$chr}{$pos_tmp};
                                                delete $nanovar{$chr}{$pos_tmp};
                                                $nanovar{$chr}{$new_pos} = $tmp;
                                            }
                                            if (exists($svim{$chr}{$pos_tmp}))
                                            {
                                                my $tmp = $svim{$chr}{$pos_tmp};
                                                delete $svim{$chr}{$pos_tmp};
                                                $svim{$chr}{$new_pos} = $tmp;
                                            }
                                            if (exists($nanosv{$chr}{$pos_tmp}))
                                            {
                                                my $tmp = $nanosv{$chr}{$pos_tmp};
                                                delete $nanosv{$chr}{$pos_tmp};
                                                $nanosv{$chr}{$new_pos} = $tmp;
                                            }
                                        }       
                                                                               
                                        next NANOSV;
                                    }
                                    elsif ($min eq '1')
                                    {                   
                                        $min = '-1';
                                        if ($score > $score_prev)
                                        {
                                            $score_prev = $score;
                                            $pos_prev = $pos_tmp;
                                        }
                                        goto POS_ALMOST2ha;
                                    }
                                    else
                                    {
                                        $min = '1';
                                        if ($score > $score_prev)
                                        {
                                            $score_prev = $score;
                                            $pos_prev = $pos_tmp;
                                        }
                                    }
                                    $v++;
                                    goto POS_ALMOST2ha;
                                }
                            }
                        }
                        elsif ($min eq '1')
                        {                   
                            $min = '-1';
                            goto POS_ALMOST2h;
                        }
                        else
                        {
                            $min = '1';
                        }
                        $v++;
                    }
                }
            }
        }
    }
}   
    
    
#Print SVs------------------------------------------------------------------

my $output_sniffles = "Sniffles_".$output_file2;
my $sniffles_count = keys %sniffles;
if ($sniffles_count > 0)
{
    open(SNIFFLES, ">" .$output_sniffles) or die "\nCan't open file $output_sniffles, $!\n";
}
my $output_pbsv = "pbsv_".$output_file2;
my $pbsv_count = keys %pbsv;
if ($pbsv_count > 0)
{
    open(PBSV, ">" .$output_pbsv) or die "\nCan't open file $output_pbsv, $!\n";
}
my $output_nanovar = "NanoVar_".$output_file2;
my $nanovar_count = keys %nanovar;
if ($nanovar_count > 0)
{
    open(NANOVAR, ">" .$output_nanovar) or die "\nCan't open file $output_nanovar, $!\n";
}
my $output_svim = "SVIM_".$output_file2;
my $svim_count = keys %svim;
if ($svim_count > 0)
{
    open(SVIM, ">" .$output_svim) or die "\nCan't open file $output_svim, $!\n";
}
my $output_nanosv = "NanoSV_".$output_file2;
my $nanosv_count = keys %nanosv;
if ($nanosv_count > 0)
{
    open(NANOSV, ">" .$output_nanosv) or die "\nCan't open file $output_nanosv, $!\n";
}

my $output_cutesv = "cuteSV_".$output_file2;
my $cutesv_count = keys %cutesv;
if ($cutesv_count > 0)
{
    open(CUTESV, ">" .$output_cutesv) or die "\nCan't open file $output_cutesv, $!\n";
}
my $datetime = localtime();   
print COMBINED2 "##fileformat=VCFv4.2\n";
print COMBINED2 "##fileDate=".$datetime."\n";
print COMBINED2 "##source=combiSV-v2.1\n";
print COMBINED2 "##ALT=<ID=DEL,Description=\"Deletion\">\n";
print COMBINED2 "##ALT=<ID=INV,Description=\"Inversion\">\n";
print COMBINED2 "##ALT=<ID=DUP,Description=\"Duplication\">\n";
print COMBINED2 "##ALT=<ID=INS,Description=\"Insertion\">\n";
##ALT=<ID=BND,Description="Breakend">
print COMBINED2 "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
print COMBINED2 "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n";
print COMBINED2 "##INFO=<ID=SVCALLERS,Number=.,Type=String,Description=\"SV callers that support this SV\">\n";
print COMBINED2 "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
print COMBINED2 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n";

print COMBINED "#CHROM\tPOS\tSVLENGTH\tTYPE\tVARHAP\n";
my %no_number;

my $INS_count = '0';
my $DEL_count = '0';
my $DUP_count = '0';
my $INV_count = '0';
my $BND_count = '0';
my $id_count = '1';

foreach my $chr2 (sort {$a <=> $b} keys %SVs)
{                             
    if ($chr2 =~ m/^\d+$/)
    {
        foreach my $pos2 (sort {$a <=> $b} keys %{$SVs{$chr2}})  
        {
            my $SV_callers = "";
            my @list = split /\t/, $SVs{$chr2}{$pos2};
            if ($count{$chr2}{$pos2} >= $count_tools || ($count{$chr2}{$pos2} > 1 && $list[3] eq "INV" && $list[4] eq "1/1"))
            {
                $count++;               
                if ($list[3] eq "INS")
                {
                    $INS_count++;
                }
                elsif ($list[3] eq "DEL")
                {
                    $DEL_count++;
                }
                elsif ($list[3] eq "DUP")
                {
                    $DUP_count++;
                }
                elsif ($list[3] eq "INV")
                {
                    $INV_count++;
                }
                elsif ($list[3] eq "BND")
                {
                    $BND_count++;
                }
                
                if ($sniffles_count > 0)
                {
                    if (exists($sniffles{$chr2}{$pos2}))
                    {
                        print SNIFFLES $sniffles{$chr2}{$pos2}."\n";
                        $SV_callers = "Sniffles";
                    }
                }
                if ($pbsv_count > 0)
                {
                    if (exists($pbsv{$chr2}{$pos2}))
                    {
                        print PBSV $pbsv{$chr2}{$pos2}."\n";
                        if ($SV_callers eq "")
                        {
                            $SV_callers = "pbsv";
                        }
                        else
                        {
                            $SV_callers .= ",pbsv";
                        }                       
                    }
                }
                if ($cutesv_count > 0)
                {
                    if (exists($cutesv{$chr2}{$pos2}))
                    {
                        print CUTESV $cutesv{$chr2}{$pos2}."\n";
                        if ($SV_callers eq "")
                        {
                            $SV_callers = "cutesv";
                        }
                        else
                        {
                            $SV_callers .= ",cutesv";
                        }                       
                    }
                }
                if ($nanovar_count > 0)
                {
                    if (exists($nanovar{$chr2}{$pos2}))
                    {
                        print NANOVAR $nanovar{$chr2}{$pos2}."\n";
                        if ($SV_callers eq "")
                        {
                            $SV_callers = "NanoVar";
                        }
                        else
                        {
                            $SV_callers .= ",NanoVar";
                        }     
                    }
                }
                if ($svim_count > 0)
                {
                    if (exists($svim{$chr2}{$pos2}))
                    {
                        print SVIM $svim{$chr2}{$pos2}."\n";
                        if ($SV_callers eq "")
                        {
                            $SV_callers = "SVIM";
                        }
                        else
                        {
                            $SV_callers .= ",SVIM";
                        }   
                    }
                }
                if ($nanosv_count > 0)
                {
                    if (exists($nanosv{$chr2}{$pos2}))
                    {
                        print NANOSV $nanosv{$chr2}{$pos2}."\n";
                        if ($SV_callers eq "")
                        {
                            $SV_callers = "NanoSV";
                        }
                        else
                        {
                            $SV_callers .= ",NanoSV";
                        }   
                    }
                }
                
                print COMBINED $SVs{$chr2}{$pos2}."\n";
                my @list_tmp = split /\t/, $SVs{$chr2}{$pos2};
                print COMBINED2 $list_tmp[0]."\t".$list_tmp[1]."\tid.".$id_count."\tN\t<".$list_tmp[3].">\t.\tPASS\tSVTYPE=".$list_tmp[3].
                ";SVLEN=".$list_tmp[2].";SVCALLERS=".$SV_callers."\tGT\t".$list_tmp[4]."\n";
                
                $id_count++;
                delete $SVs{$chr2}{$pos2};
            }
        }  
    }
}
foreach my $chr2 (sort {$a <=> $b} keys %SVs)
{                             
    foreach my $pos2 (sort {$a <=> $b} keys %{$SVs{$chr2}})  
    {
        my $SV_callers = "";
        my @list = split /\t/, $SVs{$chr2}{$pos2};
        if ($count{$chr2}{$pos2} >= $count_tools || ($count{$chr2}{$pos2} > 1 && $list[3] eq "INV" && $list[4] eq "1/1"))
        {
            $count++;
            
            if ($list[3] eq "INS")
            {
                $INS_count++;
            }
            elsif ($list[3] eq "DEL")
            {
                $DEL_count++;
            }
            elsif ($list[3] eq "DUP")
            {
                $DUP_count++;
            }
            elsif ($list[3] eq "INV")
            {
                $INV_count++;
            }
            elsif ($list[3] eq "BND")
            {
                $BND_count++;
            }
            
            if ($sniffles_count > 0)
            {
                if (exists($sniffles{$chr2}{$pos2}))
                {
                    print SNIFFLES $sniffles{$chr2}{$pos2}."\n";
                    $SV_callers = "Sniffles";
                }
            }
            if ($pbsv_count > 0)
            {
                if (exists($pbsv{$chr2}{$pos2}))
                {
                    print PBSV $pbsv{$chr2}{$pos2}."\n";
                    if ($SV_callers eq "")
                    {
                        $SV_callers = "pbsv";
                    }
                    else
                    {
                        $SV_callers .= ",pbsv";
                    }       
                }
            }
            if ($cutesv_count > 0)
            {
                if (exists($cutesv{$chr2}{$pos2}))
                {
                    print CUTESV $cutesv{$chr2}{$pos2}."\n";
                    if ($SV_callers eq "")
                    {
                        $SV_callers = "cutesv";
                    }
                    else
                    {
                        $SV_callers .= ",cutesv";
                    }                       
                }
            }
            if ($nanovar_count > 0)
            {
                if (exists($nanovar{$chr2}{$pos2}))
                {
                    print NANOVAR $nanovar{$chr2}{$pos2}."\n";
                    if ($SV_callers eq "")
                    {
                        $SV_callers = "NanoVar";
                    }
                    else
                    {
                        $SV_callers .= ",NanoVar";
                    }     
                }
            }
            if ($svim_count > 0)
            {
                if (exists($svim{$chr2}{$pos2}))
                {
                    print SVIM $svim{$chr2}{$pos2}."\n";
                    if ($SV_callers eq "")
                    {
                        $SV_callers = "SVIM";
                    }
                    else
                    {
                        $SV_callers .= ",SVIM";
                    }   
                }
            }
            if ($nanosv_count > 0)
            {
                if (exists($nanosv{$chr2}{$pos2}))
                {
                    print NANOSV $nanosv{$chr2}{$pos2}."\n";
                    if ($SV_callers eq "")
                    {
                        $SV_callers = "NanoSV";
                    }
                    else
                    {
                        $SV_callers .= ",NanoSV";
                    }   
                }
            }
                
            print COMBINED $SVs{$chr2}{$pos2}."\n";
            my @list_tmp = split /\t/, $SVs{$chr2}{$pos2};
            print COMBINED2 $list_tmp[0]."\t".$list_tmp[1]."\tid.".$id_count."\tN\t<".$list_tmp[3].">\t.\tPASS\tSVTYPE=".$list_tmp[3].
            ";SVLEN=".$list_tmp[2].";SVCALLERS=".$SV_callers."\tGT\t".$list_tmp[4]."\n";
            
            $id_count++;
            delete $SVs{$chr2}{$pos2};
        } 
    }
}
print "Combined SVs :  ".$count."\n";
print "Insertions   :  ".$INS_count."\n";
print "Deletions    :  ".$DEL_count."\n";
print "Duplications :  ".$DUP_count."\n";
print "Inversions   :  ".$INV_count."\n";
print "BND          :  ".$BND_count."\n\n";


close SNIFFLES;
close PBSV;
close CUTESV;
close NANOVAR;
close SVIM;
close NANOSV;
close INPUT_SNIFFLES;
close INPUT_CUTESV;
close INPUT_SVIM;
close INPUT_NANOVAR;
close INPUT_PBSV;
close COMBINED;
close COMBINED2;
close INPUT_NANOSV;
