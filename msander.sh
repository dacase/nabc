# Source this script to add the variables necessary for msander.

# determine file path of this script (credit http://unix.stackexchange.com/questions/96203/find-location-of-sourced-shell-script)
if [ -n "$BASH_SOURCE" ]; then
    this_script="$BASH_SOURCE"
elif [ -n "$DASH_SOURCE" ]; then
    this_script="$DASH_SOURCE"
elif [ -n "$ZSH_VERSION" ]; then
    setopt function_argzero
    this_script="$0"
elif eval '[[ -n ${.sh.file} ]]' 2>/dev/null; then
    eval 'this_script=${.sh.file}'
else
    echo 1>&2 "Unsupported shell. Please use bash, dash, ksh93 or zsh."
    exit 2
fi

# assume that the python in the PATH is the correct version:
pythonv=`python --version 2>&1 | cut -c 8-10`

export MSANDERHOME=$(cd "$(dirname "$this_script")"; pwd)
export PATH="$MSANDERHOME/bin:$PATH"
if [ -z "$PYTHONPATH" ]; then
    export PYTHONPATH="$MSANDERHOME/lib/python${pythonv}/site-packages"
else
    export PYTHONPATH="$MSANDERHOME/lib/python${pythonv}/site-packages:$PYTHONPATH"
fi
