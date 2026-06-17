#!/usr/bin/perl
#
# Make module dependencies
#

$outfile="moddepend";
unlink $outfile;
open(OUTPUT,">$outfile");

foreach $file (@ARGV){
    $uses="";
    @contains=();
    open(INPUT,"<$file") || die "Can't open $file\n";
    while(<INPUT>){

# Find modules referenced in this file
	if(/^ +use ([A-Za-z_0-9\-]+)/i){
	    $uses="$uses $1";
	}

# Find modules contained in this file
	if(/^ +module ([A-Za-z_0-9\-]+)/i){
	    push(@contains,$1);
	}

    }
    close(INPUT);

# If module is both contained and referenced in this file, we don't
# need to worry about it
    foreach $mod (@contains){
	$uses =~ s/$mod//ig;
    }
# Remove mpi module
    $uses =~ s/mpi//ig;
    if( ! $uses){ next;}

    $file =~ s/\.(f|F|f90|F90) *$/.o/g;
    $uses=lc $uses;
    $uses =~ s/( |^)([A-Za-z_0-9\-]+)/$1modules\/$2.mod/ig;
    print OUTPUT "$file: $uses\n\n";

}
close(OUTPUT);
