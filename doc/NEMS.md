NEMS Directory Structure
========================

The NEMS directory contains the source code and test scripts for the
NEMS.  Most of the documentation is in the `doc` subdirectory or in
the `../doc/` directory.

What is NOT in the NEMS is the source code to the components, the
choice of build environment and some aspects of running the system.
Those parts of the NEMS are in a directory above this, at the
so-called "app" level.  Further documentation, specific to the app, is
also at the app level.

Within NEMS resides:

* `exe` - NEMS.x and other executables built from `src`
* `src` - main program for NEMS
 * `atmos` - The ATM component code and its internal state.  This will soon be removed
  * `gfs` - All files directly related to the global spectal model
  * `nmm` - All files related to the non-hydrostatic multi-scale model
  * `fim` - all files related to the FIM
  * `gen` - All files related to the GEN
 * `share` - Assorted source shared by all model cores.
 * `ENS_Cpl` - The Ensemble coupler directory.
 * `post` - contains subroutines to conenct NEMS and the Unified Post
   Processor (UPP), no UPP source code included
 * `chem/gocart` - The GOCART source files.
 * `conf` - various compliation specifications

* `tests` - test execution scripts.  Most of these will soon be removed or
   moved to the app level.
* `compsets` - configuration for the NEMSCompsetRun.  This will soon be
   moved to the app level.
* `doc` - documentation.
* `NEMSAppBuilder` - a script to build NEMS, as discussed elsewhere in the
  documentation
* `NEMSCompsetRun` - script to run NEMS