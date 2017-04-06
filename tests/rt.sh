#! /bin/bash

# NOTE: This script is a bash script.  It uses bash-specific features
# and cannot run correctly under other shells.

function die {
    echo "$*" 1>&2
    exit 1
}

function usage {
  set +x
  echo
  echo "Usage: $0 -c <subset> | -f | -s | -t <subset> | -n /path/to/baseline | -h"
  echo
  echo "  -c <subset> = create new baseline results for <subset>"
  echo "  -f  run full suite of regression tests"
  echo "  -s  run standard suite of regression tests (same as -t standard)"
  echo "  -t <subset> = runs specified subset of tests"
  echo "  -n /path/to/baseline = specify path to REGRESSION_TEST directory"
  echo "  -h  display this help"
  echo "  -r PLATFORM:/path/to/run"
  echo "      rerun past suite without regenerating it"
  echo "  -p project = set the project to use for cpu time"
  echo
  echo "Common <subset>s: gfs, nmm, slg, wam, debug"
  echo "See ../../compsets/all.input for a full list"
  echo
  if [[ ! -z "$*" ]] ; then
      echo "SCRIPT ABORTING: $*" 1>&2
  fi
  exit 1
}

if [[ -d NEMS/tests ]] ; then
    cd NEMS/tests
elif [[ -s tests/rtgen ]] ; then
    cd tests
elif [[ ! -s rtgen ]] ; then
    die Run this script from the NEMS/tests directory.
fi

set_info=''
cmd='./rtgen -S '
baseline=NO
rerun=NO
RUNDIR=''
PLATFORM_NAME=''

export PYTHONPATH=$( pwd -P )/produtil/ush${PYTHONPATH:+:$PYTHONPATH}

while getopts ":c:fst:n:hr:p:" opt; do
    case $opt in
        r)
            if [[ $OPTARG =~ ^([a-zA-Z][a-zA-Z_.0-9]*):(.+)$ ]] ; then
                PLATFORM_NAME="${BASH_REMATCH[1]}"
                RUNDIR="${BASH_REMATCH[2]}"
                rerun=YES
            else
                echo "${BASH_REMATCH[@]}"
                usage "Rerun argument must be PLATFORM_NAME:RUNDIR"
            fi
            ;;
        p)
            cmd="$cmd -p $OPTARG" 
            ;;
        c)
            if [[ ! -z "$set_info" ]] ; then
                usage "Only one of -c, -s, -t, or -f can be used."
            fi
            set_info="$OPTARG"
            cmd="$cmd -b"
            baseline=YES
            ;;
        s)
            if [[ ! -z "$set_info" ]] ; then
                usage "Only one of -c, -s, -t, or -f can be used."
            fi
            set_info='standard'
            ;;
        f)  
            if [[ ! -z "$set_info" ]] ; then
                usage "Only one of -c, -s, -t, or -f can be used."
            fi
            set_info=' '
            ;;
        h)
            usage
            ;;
        t)
            if [[ ! -z "$set_info" ]] ; then
                usage "Only one of -c, -s, -t, or -f can be used."
            fi
            set_info="$OPTARG"
            ;;
        n)
            cmd="$cmd -n  $OPTARG"
            ;;
        \?)
            usage "Unknown option -$OPTARG"
            ;;
        :)
            usage "Option -$OPTARG requires an argument."
            ;;
    esac
done

if [[ ! -z "$set_info" ]] ; then
    cmd="$cmd $set_info"
fi

if [[ -z "$set_info" && "$rerun" == NO ]] ; then
    usage "At least one of -n, -c, -t, -f or -s must be specified."
fi

if [[ ! -z "$model" && ! -z "$fullstd" ]] ; then
    cmd="$cmd 'union($model,$fullstd)'"
elif [[ ! -z "$model" || ! -z "$fullstd" ]] ; then
    cmd="$cmd '$model$fullstd'"
fi

# Arcane bash magic below.  This is what it does:
#
#  $cmd
#    ^-- Runs rtgen
#
#  $cmd < /dev/null
#          ^-- Workaround for Jet bug: do not send stdin on a batch node
#
#  $cmd < /dev/null | tee >(cat >&2)
#                         ^--- Copy stdout to stderr
#
#  $cmd < /dev/null | tee >(cat >&2) | grep 'RUNDIR=' | tail -1
#            Get the RUNDIR variable ---^---------------^
#
# The result is that $cmd's stdout will go to the terminal AND to the
# "grep | tail".  The last line will end up in $result
#
# result='RUNDIR=/path/to/rundir  PLATFORM_NAME=wcoss.cray'

SUCCESS=PASS

if [[ "$rerun" == NO ]] ; then
    # When we're not rerunning, we need to generate a new workflow area.

    echo "rt.sh: will run $cmd"

    result=$( $cmd < /dev/null | tee >(cat >&2)  | grep 'RUNDIR=' | tail -1 )
    eval $result
    
    if [[ -z "$RUNDIR" ]] ; then
        die "ERROR: The rtgen script failed (no RUNDIR).  See above for details."
    fi
    if [[ -z "$PLATFORM_NAME" ]] ; then
        die "ERROR: The rtgen script failed (no PLATFORM_NAME).  See above for details."
    fi
else
    : # if we get here, $cmd is not used
fi

cmd="$RUNDIR/rtrun --loop -v"
echo "rt.sh: will run $cmd"

$cmd # runs for a long time
status="$?"

if [[ "$status" != 0 ]] ; then
    SUCCESS=FAIL
fi

LOGDIR="../../log/report-$PLATFORM_NAME-log"
REPORT="$LOGDIR/rtreport.txt"

if [[ "$baseline" == NO ]] ; then
    echo "rt.sh: time to generate the report"
    if [[ ! -d "$RUNDIR" ]] ; then
        echo "rt.sh: run area does not exist: $RUNDIR"
        die "rt.sh: workflow did not run there"
    fi
    echo "rt.sh: copy build logs to $LOGDIR"
    mkdir -p "$LOGDIR"
    if ( ! cp -fp "$RUNDIR/tmp/log/build"* "$LOGDIR/." ) ; then
        die "rt.sh: could not copy build logs.  Missing file?  I/O error?"
    fi
    cmd="$RUNDIR/rtreport > $REPORT"
    echo "rt.sh: run $cmd"
    "$RUNDIR/rtreport" > "$REPORT"
    repstatus="$?"
    echo "rt.sh: status $repstatus"
    if [[ "$repstatus" == 0 ]] ; then
        if ( ! grep 'REGRESSION TEST WAS SUCCESSFUL' "$REPORT" > /dev/null ) ; then
            echo "rt.sh: report says at least one test failed."
            echo "rt.sh: for details, look in $( pwd )/$REPORT"
            SUCCESS=FAIL
        else
            echo "rt.sh: report says test succeeded."
        fi
    else
        die "rt.sh: could not copy report.  Missing file?  I/O error?"
    fi
fi

echo
echo "The test directory is"
echo "    $RUNDIR"
echo "Execution area: $RUNDIR/tmp"
echo "Temporary log files: $RUNDIR/tmp/log"
echo "Subversion log files: $( pwd )/$LOGDIR"
echo "Exit status of rtrun: $status"
echo
echo "TEST RESULT: $SUCCESS"
echo

if [[ "$SUCCESS" == PASS ]] ; then
    echo "Success: rejoice!  All tests passed."
else
    echo "I deeply apologize, but at least one test failed."
    echo "You should fix the problem, and then use rtrewind"
    echo "and rtrun to re-run failed tests."
fi

exit $status
