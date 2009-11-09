#!/usr/bin/perl
#
# Reads a Grid Component Registry with definitions of INTERNAL, IMPORT and
# EXPORT states and generate the necessary code fragments and resource files
# with the appropriate MAPL calls for setting up such states, producing
# output, etc.
#
# REVISION HISTORY:
# 16Aug2006  da Silva  Initial version borrowing code from fdp.
# 04Oct2007  da Silva  Mods for new MAPL.
#
#........................................................................

use File::Basename;
use Getopt::Std;         # command line options

$Iam = "mapl_acg";

# Defaults and such
# -----------------
  init();

# Main loop over file
# -------------------
$pFlag = 0;
$iFlag = 0;
$xFlag = 0;

LINE: while (<>) {

    chomp($_);

    next LINE if (/^\#/ );  # skip comment lines

#   Component name?
#   ---------------
    if ( /COMP_NAME\:/ ) {
        @tokens = split();
        $gcname = $tokens[1] unless $opt_n; # command line takes precedence
        $gcName = $gcname unless ( $gcName = $opt_N );
	next LINE;
    }
         
#   Beginning of table? 
#   -------------------
    if ( /^\<InternalSpec/ || /^\<$gcname\:\:InternalSpec/ ) {
	$pFlag = 1;
	$pName = getVal ( 'name', $_ );
	$pCols = getVal ( 'cols', $_ );
	@pCols = split(',',$pCols);
	next LINE;
    }
    if ( /^\<ImportSpec/ || m/^\<$gcname\:\:ImportSpec/ ) {
	$iFlag = 1;
	$iName = getVal ( 'name', $_ );
	$iCols = getVal ( 'cols', $_ );
	@iCols = split(',',$iCols);
	next LINE;
    }

#   ------------------------------------------------
#   Note: prefered syntax is 
#         <ExportSpec name="XX" cols="..."> 
#         ... table contents ...
#         </ExportSpec>
#   although we still support the old syntax
#   for backward compatibility. 
#         <XX::ExportSpec cols="..."> 
#         ... table contents ...
#         </XX::ExportSpec>
#   ------------------------------------------------

    if ( /^\<ExportSpec/ ) {
	$xFlag = 1;
	$xName = getVal ( 'name', $_ );
	$xCols = getVal ( 'cols', $_ );
	@xCols = split(',',$xCols);
	next LINE;
    }
    if ( /^\<$gcname\:\:ExportSpec/ ) { # backward compatibility
	$xFlag = 1;
	$xName = "$gcname";
	$xCols = "short_name,units,dims,vlocation,stat,refresh_interval,averaging_interval,num_subtiles,long_name";
	@xCols = split(',',$xCols);
	next LINE;
    }

#   End of table?
#   ------------
    if ( /^\<\/InternalSpec\>/ || /^\<\/$gcname\:\:InternalSpec/ ) {
	$pFlag = 0;
	next LINE;
    }
    if ( /^\<\/ImportSpec\>/ || /^\<\/$gcname\:\:ImportSpec/ ) {
	$iFlag = 0;
	next LINE;
    }
    if ( /^\<\/ExportSpec\>/ ||/^\<\/$gcname\:\:ExportSpec/ ) {
	$xFlag = 0;
	next LINE;
    }

#   Internal table record
#   ---------------------
    if ( $pFlag ) {

        push @pTable, getRec('Internal',$_,@pCols);
	next LINE;

#   Import table record
#   -------------------
    } elsif ( $iFlag ) {

        push @iTable, getRec('Import',$_,@iCols);
	next LINE;

#   Export table record
#   -------------------
    } elsif ( $xFlag ) {

        push @xTable, getRec('Export',$_,@xCols);
	next LINE;
       
    } 

}

#   Parsing error
#   -------------
    die "Internal Table not terminated" if ( $pFlag );
    die   "Import Table not terminated" if ( $iFlag );
    die   "Export Table not terminated" if ( $xFlag );

#   Set output file names based on component name
#   ---------------------------------------------
    FileNames();

#   Concatenate all tables for GetPointer fragment
#   ----------------------------------------------
    push @aTable, @pTable;
    push @aTable, @iTable;
    push @aTable, @xTable;

#   Write code fragments with State Specs
#   -------------------------------------
    write_InternalSpec() if ( @pTable );
    write_ImportSpec()   if ( @iTable );
    write_ExportSpec()   if ( @xTable );

#   Write fragment with Declarations and MAPL_GetPointer()
#   ------------------------------------------------------
    write_GetPointer()   if ( @aTable );

#   Write fragments for history 
#   ---------------------------
    write_History()      if ( @xTable );

#   Write skeleton Grid Component
#   -----------------------------
    write_skeletonGC() if ( $opt_s );

exit(0);

#.........................................................................

sub init {

# Command line options
# --------------------
  getopts('c:CEFvmdn:N:i:x:g:h:Gs');

  if ( $opt_G ) {
      $MAPL = "GEOS";
  } else {
      $MAPL = "MAPL";
  }

  usage() if ( $opt_h || $#ARGV < 0 );

 $ifilen = $ARGV[0] unless ( $ifilen = $opt_i );
($base,$path,$suffix) = fileparse($ifilen,'\..*');

$gcname  = $base unless ($gcname = $opt_n );

$gcname =~ s/_REGISTRY//g;
$gcname =~ s/_Registry//g;
$gcname =~ s/_registry//g;

$gcName = $gcname unless ( $gcName = $opt_N );

$mangle = $opt_m;

}

#.........................................................................

sub FileNames {

$pFilen = "$gcname"."_InternalSpec___.h" unless ( $pFilen=$opt_p );
$iFilen = "$gcname"."_ImportSpec___.h"   unless ( $iFilen=$opt_i );
$xFilen = "$gcname"."_ExportSpec___.h"   unless ( $xFilen=$opt_x );
$gFilen = "$gcname"."_GetPointer___.h"   unless ( $gFilen=$opt_g );
$hFilen = "$gcname"."_History___.rc"     unless ( $hFilen=$opt_h );
$gcFilen = "$gcname"."_GridCompMod___.F90"  unless ( $gcFilen=$opt_c );

}

#.........................................................................

sub getVal { # Needs improvement - requires something like: name="whatever" 
        my $varn = shift;
        my $line = shift;
        my @toks = split ( "$varn", $line ); 
	my $rhs  = $toks[1];
        @toks = split('"',"$rhs");
        return $toks[1];
}

#.........................................................................

sub getRec {
        my $Type  = shift;
        my $line = shift;
        my @cols = @_;
        my @tokens = split('\|', $line);
        my $i = 0; 
        my %record = undef;
        $record{Type} = $Type;
        foreach $key ( @cols ) {
	    $record{$key} = trim($tokens[$i]);
	    $i++;
	}
        return \%record;
    }

#.........................................................................

sub trim {
    my $str = shift;
    my @words = split(' ',$str);
    $str = join ' ', @words;
    return $str;
}

#.........................................................................

sub write_InternalSpec {

    open(FILE,">$pFilen") or die "cannot open InternalSpec file $pFilen";
    print "Building INTERNAL Spec file $pFilen\n" if ( $opt_v ); 
    Preamble('!');
    $n = @pTable - 1;
    for $i ( 0..$n ) {
        write_SpecRecord ( "Internal", $pTable[$i]); 
    }
    close(FILE);
}

#.........................................................................

sub write_ImportSpec {

    open(FILE,">$iFilen") or die "cannot open ImportSpec file $iFilen";
    print "Building IMPORT Spec file $iFilen\n" if ( $opt_v ); 
    Preamble('!');
    $n = @iTable - 1;
    for $i ( 0..$n ) {
        write_SpecRecord ( "Import", $iTable[$i]); 
    }
    close(FILE);
}

#.........................................................................

sub write_ExportSpec {

    open(FILE,">$xFilen") or die "cannot open ExportSpec file";

    print "Building EXPORT Spec file $xFilen\n" if ( $opt_v ); 

    Preamble('!');

    $n = @xTable - 1;
    for $i ( 0..$n ) {

        %r = %{$xTable[$i]};

        $name  = $r{short_name};
        $long  = $r{long_name};
        $units = $r{units};

        if ( $r{dims} eq "xyz" ) { $dims = "MAPL_DimsHorzVert"; }
	else {                     $dims = "MAPL_DimsHorzOnly"; }

        if ( $r{vlocation} eq "C" ) { $vloc = "MAPL_VLocationCenter"; }
     elsif ( $r{vlocation} eq "E" ) { $vloc = "MAPL_VLocationEdge"; }
	else                   { $vloc = "MAPL_VLocationNone"; }

        print FILE <<EOF;

     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = '$name',  &
        LONG_NAME          = '$long',  &
        UNITS              = '$units', &
        DIMS               = $dims,    &
        VLOCATION          = $vloc,    &
EOF

     $stat = $r{stat};

     if ( "$stat" eq "B" ) { $stat = 'MAPL_BundleItem'; } # what else?


     if ( "$stat" ne "" ) {
        print FILE <<EOF;
        STAT               = $stat,    &
EOF
     }

     $rfr  = $r{refresh_interval};
     if ( "$rfr" ne "" ) {
        print FILE <<EOF;
        REFRESH_INTERVAL   = $rfr,     &
EOF
     }

     $avg  = $r{averaging_interval};
     if ( "$avg" ne "" ) {
        print FILE <<EOF;
        AVERAGING_INTERVAL = $avg,     &
EOF
     }

     $subt = $r{num_subtitles};
     if ( "$subt" ne "" ) {
        print FILE <<EOF;
        NUM_SUBTILES       = $subt,    &
EOF
     }

#    End of item
#    -----------
     print FILE <<EOF
                                                       RC=STATUS  )
     VERIFY_(STATUS)

EOF

    }

    close(FILE);

}

#.........................................................................

sub write_SpecRecord {

    my $which  = shift;
    my $record = shift;

    $method = "MAPL_Add"."$which"."Spec";
    %r = %{$record};

    $name  = $r{short_name};
    $long  = $r{long_name};
    $units = $r{units};

    $long =~ s/__ENSEMBLE__//g; # filter this out

    if ( $r{dims} eq "xyz" ) { $dims = "MAPL_DimsHorzVert"; }
    else {                     $dims = "MAPL_DimsHorzOnly"; }

       if ( $r{vlocation} eq "C" ) { $vloc = "MAPL_VLocationCenter"; }
    elsif ( $r{vlocation} eq "E" ) { $vloc = "MAPL_VLocationEdge"; }
    else                           { $vloc = "MAPL_VLocationNone"; }

    print FILE <<EOF;

     call $method(GC, &
        SHORT_NAME         = '$name',  &
        LONG_NAME          = '$long',  &
        UNITS              = '$units', &
        DIMS               = $dims,    &
        VLOCATION          = $vloc,    &
EOF

     $stat = $r{stat};

     if ( "$stat" eq "B" ) { $stat = 'MAPL_BundleItem'; } # what else?


     if ( "$stat" ne "" ) {
        print FILE <<EOF;
        STAT               = $stat,    &
EOF
     }

     $rfr  = $r{refresh_interval};
     if ( "$rfr" ne "" ) {
        print FILE <<EOF;
        REFRESH_INTERVAL   = $rfr,     &
EOF
     }

     $avg  = $r{averaging_interval};
     if ( "$avg" ne "" ) {
        print FILE <<EOF;
        AVERAGING_INTERVAL = $avg,     &
EOF
     }

     $subt = $r{num_subtitles};
     if ( "$subt" ne "" ) {
        print FILE <<EOF;
        NUM_SUBTILES       = $subt,    &
EOF
     }

#    End of item
#    -----------
     print FILE <<EOF
                                                       RC=STATUS  )
     VERIFY_(STATUS)

EOF

}


#.........................................................................

sub write_GetPointer {

    if ( $opt_F ) { write_GetPointer_FlatArray(); }
    else          { write_GetPointer_ChemArray(); }

}

#.........................................................................

sub write_GetPointer_FlatArray {

    open(FILE,">$gFilen") or die "cannot open ExportSpec file";

    print "Building F90 code fragment file $gFilen (flat arrays)\n" 
        if ( $opt_v ); 

    Preamble('!');

#   Declarations:
#   -------------
    print FILE <<EOF;

!       Local arrays referencing the Import/Export states
!       -------------------------------------------------
EOF

    $n = @aTable - 1;
VAR:  for $i ( 0..$n ) {

        %r = %{$aTable[$i]};

        $Type  = $r{Type};
        $name  = $r{short_name};
        $long  = $r{long_name};

        if ( $r{dims} eq "xyz" ) { $rank = "(:,:,:)"; }
	else {                     $rank = "(:,:)  "; }

        if ( $opt_E ) {
          if ( $long =~ / __ENSEMBLE__/ ) {
	    $long =~ s/__ENSEMBLE__//g; # filter this out
          } else {
            next VAR;
          }
	}

        print FILE <<EOF;
        real, pointer, dimension$rank :: $name ! $Type: $long        
EOF
        print "#define $name \n" if ( $opt_d );

       }

        print FILE <<EOF;

!       Get pointers to data in state
!       -----------------------------
EOF

#   Get the pointer - non-binned variables
#   --------------------------------------
    $n = @aTable - 1;
    for $i ( 0..$n ) {

        %r = %{$aTable[$i]};
        $name  = $r{short_name};
        $long  = $r{long_name};

#       Special handle instance variables
#       ---------------------------------
        if ( $opt_E ) {
          if ( $long =~ / __ENSEMBLE__/ ) {
	    print FILE <<EOF;
        call MAPL_GetPointer ( EXPORT, $name,  '$name'//iNAME, RC=STATUS )
        VERIFY_(STATUS)
EOF
	  }
        } else {
	  print FILE <<EOF;
        call MAPL_GetPointer ( EXPORT, $name,  '$name', RC=STATUS )
        VERIFY_(STATUS)
EOF
	}
    }

    close(FILE);

}

#.........................................................................

sub write_GetPointer_ChemArray {

    open(FILE,">$gFilen") or die "cannot open ExportSpec file";

    if ( $opt_C ) {

	print "Building F90 code fragment file $gFilen (Chem Arrays)\n" 
              if ( $opt_v ); }

    else {

	print "Building F90 code fragment file $gFilen\n" 
              if ( $opt_v ); 

    }

    Preamble('!');

#   First pass: look for binned variables
#   -------------------------------------
    $n = @aTable - 1;
    VAR: for $i ( 0..$n ) {

        %r = %{$aTable[$i]};

        $name  = $r{short_name};
        $long  = $r{long_name};
        $dims  = $r{dims};

        if ( $opt_E ) {
          if ( $long =~ / __ENSEMBLE__/ ) {
	    $long =~ s/__ENSEMBLE__//g; # filter this out
          } else {
            next VAR;
          }
	}

#       Binned variables, a little trickier
#       -----------------------------------
        if ( $long =~ / Bin / ) {

	    @tokens = split('Bin',$long);
            $title = $tokens[0];
            $nBin = $tokens[1];
            $nBin =~ s/ //g;           

            $name =~ s/$nBin//g;
            $nBin = $nBin + 0;

#           For this work bins must be in ascending order
#           --------------------------------------------
            if ( "$dims" eq "xyz" ) { $v3Bin{$name} = $nBin; }
            else                    { $v2Bin{$name} = $nBin; }

            $vBin{$name} = $nBin; 
            $tBin{$name} = $title; 

            $xName[$i] = $name;
            $xBin[$i] = $nBin;

        } else {

            $xName[$i] = $name;  # needed for -C option
        }
	
	$Type{$name} = $r{Type}; # for later reference

    }

#   Declarations: bin sizes
#   -----------------------
    print FILE <<EOF;

!     Bin sizes
!     ---------
EOF
    foreach $name ( keys %vBin ) {
	$n = $vBin{$name};
        print FILE <<EOF;
      integer, parameter              :: NBIN_$name = $n ! $tBin{$name}
EOF
      print "#define $name \n" if ( $opt_d );
      }

#   Declarations: binned variables
#   ------------------------------
    print FILE <<EOF;

!     Bin-indexed Chem Arrays
!     -----------------------
EOF

    foreach $name ( keys %vBin ) {
      $n = $vBin{$name};
      $Type = $Type{$name};

      print FILE <<EOF;
      type(Chem_Array), target        ::    $name(NBIN_$name) ! $Type: $tBin{$name}
      type(Chem_Array), pointer       :: ptr$name(:)  ! $Type: $tBin{$name}
EOF
      print "#define ptr$name \n" if ( $opt_d );

  }


#   Declarations: non-binned variables
#   ----------------------------------
    print FILE <<EOF;

!     Local array referencing the Import/Export states
!     ------------------------------------------------
EOF

    $n = @aTable - 1;
    VAR2: for $i ( 0..$n ) {

        %r = %{$aTable[$i]};

        $Type  = $r{Type};
        $name  = $r{short_name};
        $long  = $r{long_name};

        if ( $r{dims} eq "xyz" ) { $rank = "(:,:,:)"; }
	else {                     $rank = "(:,:)  "; }

        if ( $opt_E ) {
          if ( $long =~ / __ENSEMBLE__/ ) {
	    $long =~ s/__ENSEMBLE__//g; # filter this out
          } else {
            next VAR2;
          }
	}

#       Binned variables
#       ----------------
        unless ( $long =~ / Bin / ) {

          if ( $opt_C ) {
	      print FILE <<EOF;
      type(Chem_Array), target        ::    $name ! $Type: $long
      type(Chem_Array), pointer       :: ptr$name ! $Type: $long
EOF
          } else {         
              print FILE <<EOF;
      real, pointer, dimension$rank :: $name ! $Type: $long
EOF
          }
          print "#define ptr$name \n" if ( $opt_d );
          }
    
    }

        print FILE <<EOF;

!     Get pointers to data in state
!     -----------------------------
EOF

#   Get the pointer - non-binned variables
#   --------------------------------------
    $n = @aTable - 1;
    VAR3: for $i ( 0..$n ) {

        %r = %{$aTable[$i]};

        $name  = $r{short_name};
        $long  = $r{long_name};
        $dims  = $r{dims};

        if ( $r{dims} eq "xyz" ) { $rank = "3d"; }
	else {                     $rank = "2d"; }

#                          ---------------
#                          Doing ensembles
#                          ---------------
        if ( $opt_E ) {

          if ( $long =~ / __ENSEMBLE__/ ) {
	    $long =~ s/__ENSEMBLE__//g; # filter this out
          } else {
            next VAR3;
          }

#         Binned variables are Chem Arrays
#         --------------------------------
	  if ( $long =~ / Bin / ) {
	      print FILE "\n      ptr$xName[$i] => $xName[$i]   ! $long\n"
                                                       if ( $xBin[$i] == 1 ); 
	      print FILE <<EOF;
      call MAPL_GetPointer ( EXPORT, $xName[$i]($xBin[$i])%data$rank,  '$name'//iNAME, RC=STATUS )
      VERIFY_(STATUS)
EOF

#         If desired, even non-binned variables are Chem Arrays
#         -----------------------------------------------------
          } elsif ( $opt_C ) {
             print FILE <<EOF;

      ptr$xName[$i] => $xName[$i]   ! $long
      call MAPL_GetPointer ( EXPORT, $xName[$i]%data$rank,  '$name'//iNAME, RC=STATUS )
      VERIFY_(STATUS)
EOF

          } else {
            print FILE <<EOF;
      call MAPL_GetPointer ( EXPORT, $name,  '$name'//iNAME, RC=STATUS )
      VERIFY_(STATUS)
EOF
          }

#                             -------------------
#                             Not doing ensembles
#                             -------------------
      } else {

#         Binned variables are Chem Arrays
#         --------------------------------
	  if ( $long =~ / Bin / ) {
	      print FILE "\n      ptr$xName[$i] => $xName[$i]   ! $long\n"
                                                       if ( $xBin[$i] == 1 ); 
	      print FILE <<EOF;
      call MAPL_GetPointer ( EXPORT, $xName[$i]($xBin[$i])%data$rank,  '$name', RC=STATUS )
      VERIFY_(STATUS)
EOF

#         If desired, even non-binned variables are Chem Arrays
#         -----------------------------------------------------
          } elsif ( $opt_C ) {
             print FILE <<EOF;

      ptr$xName[$i] => $xName[$i]   ! $long
      call MAPL_GetPointer ( EXPORT, $xName[$i]%data$rank,  '$name', RC=STATUS )
      VERIFY_(STATUS)
EOF

          } else {
            print FILE <<EOF;
      call MAPL_GetPointer ( EXPORT, $name,  '$name', RC=STATUS )
      VERIFY_(STATUS)
EOF
          }

      } # if doing ensembles...

    } # variable loop

    close(FILE);

}

#.........................................................................

sub write_History {

    open(FILE,">$hFilen") or die "cannot open ExportSpec file";

    print "Building History Spec fragment file $hFilen\n" if ( $opt_v ); 

    $Name = "$gcName"."::"."$name";

    Preamble('#');

#   2D quantities
#   -------------

        print FILE <<EOF;

  list(#)%filename:   '/dev/null/%s.$gcName.sfc.%y4%m2%d2_%h2%n2z',
  list(#)%format:     'CFIO',
  list(#)%mode:       'time-averaged',
  list(#)%frequency:  030000,
  list(#)%duration:   030000,
EOF

     $label =   'list(#)%fields:    ';
     $n = @aTable - 1;
ONE: for $i ( 0..$n ) {

        %r = %{$aTable[$i]};
        $name  = $r{short_name};
        $Name = "$gcName"."::"."$name";

        if ( "$r{dims}" eq "xyz" ) { 
	    $has3D = 1;
            next ONE;
        }

        if ( $opt_m ) {
	    print FILE "  $label '$Name'     , '$gcName'      ,  '$name'   ,\n";
        } else {
	    print FILE "  $label '$name'     , '$gcName'      ,\n";
        }

    $label =   '                   '

    }

# 3D quantities
# -------------
  if ( $has3D ) {

        print FILE <<EOF;

  list(#)%filename:   '/dev/null/%s.$gcName.eta.%y4%m2%d2_%h2%n2z',
  list(#)%format:     'CFIO',
  list(#)%mode:       'time-averaged',
  list(#)%frequency:  030000,
  list(#)%duration:   030000,
EOF

     $label =   'list(#)%fields:    ';
     $n = @aTable - 1;
VAR: for $i ( 0..$n ) {

        %r = %{$aTable[$i]};
        $name  = $r{short_name};
        $Name = "$gcName"."::"."$name";

        next VAR if ( "$r{dims}" eq "xy" );

        if ( $opt_m ) {
	    print FILE "  $label '$Name'     , '$gcName'      ,  '$name'   ,\n";
        } else {
	    print FILE "  $label '$name'     , '$gcName'      ,\n";
        }

    $label =   '                   '

    }

 }

    close(FILE);

}

#.........................................................................

sub write_skeletonGC {

    open(FILE,">$gcFilen") or die "cannot open ExportSpec file";

    print "Building Skeleton GridComponent file $gcFilen\n" if ( $opt_v ); 

    $rf90code = gcF90code(); # load template code (see below)

    chomp($today = `date '+%d%b%Y'`); # for revision history

    foreach $line ( @$rf90code ) {
	$line =~ s/MAPL_/GEOS_/g if ( $opt_G );
	$line =~ s/TODAY/$today/g;
	$line =~ s/GCNAME/$gcname/g;
	print FILE "$line";
    }

    close(FILE);

}

#.........................................................................

sub gcF90code {

    my @f90code = q {

#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GCNAME_GridCompMod - The GCNAME Grid Component
!
! !INTERFACE:
!
   Module GCNAME_GridCompMod
!
! !USES:
!
   use ESMF_Mod
   use MAPL_Mod

   Implicit NONE
   Private

! !PUBLIC MEMBER FUNCTIONS:

   Public SetServices
!
! !DESCRIPTION: 
!
!  This module implements the GCNAME Grid Component.
!
! !REVISION HISTORY:
!
!  TODAY  mapl_acg  Automatically generated code.
!
!EOP
!-------------------------------------------------------------------------

! Optional, non-MAPL internal state 
! ---------------------------------
  type GCNAME_State
     private
     integer :: place_holder
   end type GCNAME_State

  type GCNAME_WRAP
     type (GCNAME_State), pointer :: PTR => null()
  end type GCNAME_WRAP

CONTAINS

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetServices --- Sets IRF services for the GCNAME Grid Component
!
! !INTERFACE:

   subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION: Sets Initialize, Run and Finalize services. 
!
! !REVISION HISTORY:
!
!  TODAY  mapl_acg  Automatically generated code.
!
!EOP
!-------------------------------------------------------------------------

!   ErrLog Variables
!   ----------------
    character(len=ESMF_MAXSTR)      :: IAm = 'SetServices'
    integer                         :: STATUS
    character(len=ESMF_MAXSTR)      :: COMP_NAME

!   Optional, non-MAPL internal state 
!   ---------------------------------
    type (GCNAME_State), pointer    :: LegacyState  ! internal, that is
    type (GCNAME_wrap)              :: wrap

!   Configuration info
!   ------------------
    type (ESMF_Config)              :: CF

!                              ------------


!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // '::' // trim(Iam)


!                       ------------------------
!                       ESMF Functional Services
!                       ------------------------

!   Set the Initialize, Run, Finalize entry points
!   ----------------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_SETINIT,  Initialize_,  &
                                      RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint ( GC,  ESMF_SETRUN,   Run_,        &
                                      RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint ( GC,  ESMF_SETFINAL, Finalize_,  &
                                      RC=STATUS)
    VERIFY_(STATUS)

!   Wrap internal legacy state for storing in GC
!   --------------------------------------------
    allocate ( LegacyState, stat=STATUS )
    VERIFY_(STATUS)
    wrap%ptr => LegacyState

!   Store internal state in GC
!   --------------------------
    call ESMF_UserCompSetInternalState ( GC, 'GCNAME_State', wrap, STATUS )
    VERIFY_(STATUS)


!                         ------------------
!                         MAPL Data Services
!                         ------------------

!BOP

! !IMPORT STATE:
#  include "GCNAME_ImportSpec___.h"

! !EXPORT STATE:
#  include "GCNAME_ExportSpec___.h"

! !INTERNAL STATE:
#  include "GCNAME_InternalSpec___.h"

!EOP

!   Set the Profiling timers
!   ------------------------
    call MAPL_GenericStateClockAdd ( GC, name = "RUN", RC=STATUS )
    VERIFY_(STATUS)

!   Generic Set Services
!   --------------------
    call MAPL_GenericSetServices ( GC, RC=STATUS )
    VERIFY_(STATUS)

!   All done
!   --------
    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Initialize_ --- Initialize the GCNAME Grid Component
!
! !INTERFACE:
!

   subroutine Initialize_ ( GC, IMPORT, EXPORT, CLOCK, rc )

! !USES:

   implicit NONE

! !ARGUMENTS:

   type(ESMF_GridComp), intent(inout) :: GC      ! Grid Component
   type(ESMF_State), intent(inout)    :: IMPORT  ! Import State
   type(ESMF_State), intent(inout)    :: EXPORT  ! Export State
   type(ESMF_Clock),  intent(inout)   :: CLOCK   ! The clock
   integer, intent(out)               :: RC      ! Error return code:
                                                 !  0 - all is well
                                                 !  1 - 

! !DESCRIPTION: 
!
!  This routine initializes the GCNAME Grid Component.
!
! !REVISION HISTORY:
!
!  TODAY  mapl_acg  Automatically generated code.
!
!EOP
!-------------------------------------------------------------------------

!  ErrLog Variables
!  ----------------
   character(len=ESMF_MAXSTR)      :: IAm = 'Initialize_'
   integer                         :: STATUS
   character(len=ESMF_MAXSTR)      :: COMP_NAME

!  Optional, non-MAPL internal state 
!  ---------------------------------
   type (GCNAME_State), pointer  :: LegacyState   ! internal, that is
   type (GCNAME_wrap)            :: wrap

!  Configuration info
!  ------------------
   type (ESMF_Config)              :: CF

!                                 ---

!  Get pointer to states
!  ---------------------
#  include "GCNAME_GetPointer___.h"

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
   VERIFY_(STATUS)
   Iam = trim(COMP_NAME) // '::' // Iam

!  Get my internal legacy state
!  ----------------------------
   call ESMF_UserCompGetInternalState (GC, 'GCNAME_State', WRAP, STATUS)
   VERIFY_(STATUS)
   LegacyState => wrap%ptr

!  Initialize MAPL Generic
!  ------------------------
   call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS )
   VERIFY_(STATUS)

!  Add your code here...

   RETURN_(ESMF_SUCCESS)

   end subroutine Initialize_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Run_ --- Run the GCNAME Grid Component
!
! !INTERFACE:
!

   subroutine Run_ ( GC, IMPORT, EXPORT, CLOCK, rc )

! !USES:

   implicit NONE

! !ARGUMENTS:

   type(ESMF_GridComp), intent(inout) :: GC      ! Grid Component
   type(ESMF_State), intent(inout)    :: IMPORT  ! Import State
   type(ESMF_State), intent(inout)    :: EXPORT  ! Export State
   type(ESMF_CLOCK),  intent(inout)   :: CLOCK   ! The clock
   integer, intent(out)               :: RC      ! Error return code:
                                                 !  0 - all is well
                                                 !  1 - 

! !DESCRIPTION: 
!
!  This routine implements the run method for the GCNAME Grid Component.
!
! !REVISION HISTORY:
!
!  TODAY  mapl_acg  Automatically generated code.
!
!EOP
!-------------------------------------------------------------------------

!  ErrLog Variables
!  ----------------
   character(len=ESMF_MAXSTR)     :: IAm = 'Run_'
   integer                        :: STATUS
   character(len=ESMF_MAXSTR)     :: COMP_NAME

!  Optional, non-MAPL internal state 
!  ---------------------------------
   type (GCNAME_State), pointer  :: LegacyState   ! internal, that is
   type (GCNAME_wrap)            :: wrap

!  Configuration info
!  ------------------
   type (ESMF_Config)            :: CF

!                                 ---

!  Get pointer to states
!  ---------------------
#  include "GCNAME_GetPointer___.h"

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
   VERIFY_(STATUS)
   Iam = trim(COMP_NAME) // '::' // Iam

!  Get my internal legacy state
!  ----------------------------
   call ESMF_UserCompGetInternalState (GC, 'GCNAME_State', WRAP, STATUS)
   VERIFY_(STATUS)
   LegacyState => wrap%ptr

!  Add your code here...

!  Run MAPL Generic
!  ----------------
   call MAPL_GenericRun ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS )
   VERIFY_(STATUS)

!  Add more of your code here...

   RETURN_(ESMF_SUCCESS)

   end subroutine Run_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Finalize_ --- Finalize the GCNAME Grid Component
!
! !INTERFACE:
!

   subroutine Finalize_ ( GC, IMPORT, EXPORT, CLOCK, rc )

! !USES:

   implicit NONE


! !ARGUMENTS:

   type(ESMF_GridComp), intent(inout) :: GC      ! Grid Component
   type(ESMF_State), intent(inout)    :: IMPORT  ! Import State
   type(ESMF_State), intent(inout)    :: EXPORT  ! Export State
   type(ESMF_Clock),  intent(inout)   :: CLOCK   ! The clock
   integer, intent(out)               :: RC      ! Error return code:
                                                 !  0 - all is well
                                                 !  1 - 

! !DESCRIPTION: 
!
!  This routine finalizes the GCNAME Grid Component.
!
! !REVISION HISTORY:
!
!  TODAY  mapl_acg  Automatically generated code.
!
!
!EOP
!-------------------------------------------------------------------------

!  ErrLog Variables
!  ----------------
   character(len=ESMF_MAXSTR)      :: IAm = 'Finalize_'
   integer                         :: STATUS
   character(len=ESMF_MAXSTR)      :: COMP_NAME

!  Optional, non-MAPL internal state 
!  ---------------------------------
   type (GCNAME_State), pointer    :: LegacyState   ! internal, that is
   type (GCNAME_wrap)              :: wrap

!  Configuration info
!  ------------------
   type (ESMF_Config)             :: CF

!                                 ---

!  Get pointer to states
!  ---------------------
#  include "GCNAME_GetPointer___.h"

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
   VERIFY_(STATUS)
   Iam = trim(COMP_NAME) // '::' // Iam

!  Get my internal legacy state
!  ----------------------------
   call ESMF_UserCompGetInternalState (GC, 'GCNAME_State', WRAP, STATUS)
   VERIFY_(STATUS)
   LegacyState => wrap%ptr

!  Add your code here...

!  Finalize MAPL Generic
!  ---------------------
   call MAPL_GenericFinalize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS )
   VERIFY_(STATUS)

!  Add more of your code here...

   RETURN_(ESMF_SUCCESS)

   end subroutine Finalize_

 end module GCNAME_GridCompMod

};

    return \@f90code;

}

#.........................................................................

sub Preamble {

    my $c = shift;

    print FILE <<EOF;
$c                          -------------------
$c                          W  A  R  N  I  N  G
$c                          -------------------
$c
$c   This code fragment is automatically generated by $Iam.
$c   Please DO NOT edit it. Any modification made in here will be overwritten
$c   next time this file is auto-generated. Instead, enter your additions
$c   or deletions in the $gcname_Registry.rc file. 
$c

EOF

}

#.........................................................................

sub usage {

    $GetPointer = "$MAPL"."_GetPointer";
   print <<"EOF";

NAME
     mapl_acg - generates code fragments for setting up MAPL states
          
SYNOPSIS

     mapl_acg.pl [OPTIONS]  Registry_filename
          
DESCRIPTION

     Reads a Grid Component Registry with definitions of INTERNAL,
     IMPORT and EXPORT states and generate the necessary code
     fragments and resource files. The following code fragments/
     resorce file can be generated:

        Code Fragment Type     Default Output File Name
        -------------------    ------------------------
        Internal Spec Setup    gcNAME_InternalSpec___.h
        Import   Spec Setup    gcNAME_ImportSpec___.h
        Export   Spec Setup    gcNAME_ExportSpec___.h
        $GetPointer()      gcNAME_GetPointer___.h
        History streams        gcNAME_History___.rc
        Grid Component         gcNAME_GridCompMod___.F90

     The Grid Component name "gcNAME" is either specified with he -n
     option, defined in the Registry or derived from the input file
     name.  Notice that files matching "*___.h" and "*___.rc" are
     ignored by the CVS ESMA repository at NASA/GSFC and are not
     meant to be checked in.

     The F90 code fragments contain the appropriate MAPL calls for
     setting up import/export/external states. In addition, it also
     generates a code fragment (*_GetPointer.h) for automatically
     declaring and retrieving fields from these states in the form of
     simple F90 arrays, or "arrays of arrays" - read on. You can also
     generate a skeleton for constructing ESMF Grid Components from
     scratch, see below.

OPTIONS
     -c             file name for skeleton Grid Component code
     -C             make all arrays Chem Arrays
     -d             print variables #defines on the screen.
     -E             ensemble mode 
     -F             make all arrays flat FORTRAN arrays
     -G             Use GEOS_ preffix instead of MAPL_ in output
     -g             file name for code fragment with calls to $GetPointer()
     -h             file name for code fragment with History.rc fragmet
     -i             file name for code fragment with Import   state spec
     -m             Mangle variable name in history stream
     -N name        component name to be used in history stream
     -n name        name of grid component; the default is taken from the
                    rewsource file or derived from the input file name
     -p             file name for code fragment with Internal state spec
     -s             write skeleton Grid Component code
     -v             verbose mode
     -x             file name for code fragment with Export state spec

CREATING NEW GRID COMPONENTS

     To aid the construction of ESMF Grid Components from scratch a
     skeleton grid component can be generated with the option "-s".
     By the default, the output file is named "gcNAME_GridCompMod___.F90"
     and it is recomended that it be manually renamed "gcNAME_GridCompMod.F90"
     so that it can be checked in the CVS repository. (This default is
     designed to prevent one from accidentally overwriting an existing
     Grid Component file.) The recommended sequence of steps for
     creating a new MAPL-compliant Grid Component named XX is this:

     1) Prepare a registry file "XX_Registry.rc" defining all fields 
        comprising your Import/Export/Internal states.

     2) Add targets to the GNUmakefile for automatically creating the
        corresponding code fragments. Typically, one would include the 
        following:
        
        ACGS :=  XX_InternalSpec___.h XX_ImportSpec___.h XX_ExportSpec___.h\
                 XX_GetPointer___.h   XX_History___.h

        $(ACGS) : XX_Registry.rc $(ACG)
               @$(ACG) $(ACG_FLAGS) -C -G XX_Registry.rc

        $(DEPS) : $(ACGS)

     3) The first time around, one should manually create the skeleton 
        grid component:
        
        % mapl_acg.pl -v -s -G XX_Registry.rc
        % mv XX_GridCompMod___.F90 XX_GridCompMod.F90 

     4) Then, edit "XX_GridCompMod.F90" and include your new code.
        This skeleton includes hooks for setting up a LegacyState
        for those codes that have internal state information that
        cannot be expressed as ESMF states. Notice that while legacy
        states are allowed, they are not recommended. If you do not
        need legacy states, you can remove any reference to these
        states from your Grid Component. (Leaving the existing
        template in place will cause no harm).

ARRAY OF ARRAYS (CHEM_ARRAYS)

     To be written. You can safely ignore this if you are not involved
     with the GOCART development.

REGISTRY FILE FORMAT

     To be written. At this point your best bet is to look at the
     Sample_Registry.rc file included with this script (under
     MAPL_Base).  Keep in mind that the cols="..." definition is
     extremely important; all column names should be lower case, comma
     separated, and match exactly the name of the arguments in the
     MAPL_Add*Spec() routine. The tables can have fewer columns than
     the number of optional arguments in the corresponding MAPL
     routine; obviously, when columns are ommited the missing
     arguments cannot be specified and the default values will be
     taken.

SEE ALSO
     The MAPL User's Guide.

AUTHOR
     Arlindo da Silva, NASA/GSFC.

EOF

  exit(1)

 }
