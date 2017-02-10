#!/usr/bin/perl 

use strict;
#use warnings;

open (FILE, $ARGV[0]);

my (@line, $fftw, $fail, $omp, $val_n, $val_k, $n, $k);

# printing only certain n
if( $#ARGV >= 1 ) {
    $val_n = $ARGV[1]
}
if( $#ARGV == 2 ) {
    $val_k = $ARGV[2]
}

$fail = 0;
while (<FILE>) {
    chomp;
    @line = split;
    if( $_ =~ "Time to run FFTW" ) {
        $fftw = $line[5];
        if( $val_n ) { 
            if( $val_n eq $n ) { 
                if( $val_k ) {
                    if( $val_k eq $k ) { 
                        print " fftw: ",$fftw," $omp\n";
                    }
                } else {
                    print " fftw: ",$fftw," $omp\n";
                }
            }
        } else {
            print " fftw: ",$fftw," $omp\n";
        }
        $fail = 0;
    }
    if( $line[2] =~ /n=/ ) {
        if( $fail == 1 ) { 
            if( !$val_n ) { print "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 $omp\n"} 
        }
        chop $line[2];
        chop $line[3];
        $n = substr $line[2], 2;
        $k = substr $line[3], 2;
        if( $val_n ) { 
            if( $val_n eq $n ) { 
                if( $val_k )
                {
                    if( $val_k eq $k ) {
                        print $n," ",$k," ";
                    }
                } else {
                        print $n," ",$k," ";
                }
            }
        } else {
            print $n," ",$k," ";
        }
        $fail = 1;
    }
    #if( $#line == 8 and $line[0] !~ /\%/ and $line[0] != "Projected") {
    if( $#line == 7 and $line[0] !~ /\%/ and $_ !~ "Projected") {
        if( $val_n ) { 
            if( $val_n eq $n ) { 
                if( $val_k ) { 
                    if( $val_k eq $k ) { 
                        print $line[0]," ",$line[1]," ",$line[2]," ",$line[3]," ",$line[4], " ",$line[5], " ",$line[6], " ",$line[7];
                    }
                } else {
                    print $line[0]," ",$line[1]," ",$line[2]," ",$line[3]," ",$line[4], " ",$line[5], " ",$line[6], " ",$line[7];
                }
            }
        } else {
            print $line[0]," ",$line[1]," ",$line[2]," ",$line[3]," ",$line[4], " ",$line[5], " ",$line[6], " ",$line[7];
        }
    }
    if( $line[0] =~ /^OMP/ ) {
        $omp = $line[2];
    }
}
