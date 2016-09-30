#! /usr/bin/env ksh
function atparse {
    # Use fancy __ names to avoid clashes in variable expansion {var} blocks.
    typeset __text
    typeset __before
    typeset __after
    typeset __during
    for __text in "$@" ; do
        if [[ $__text =~ ^([a-zA-Z][a-zA-Z0-9_]*)=(.*)$ ]] ; then
            eval "typeset ${.sh.match[1]}"
            eval "${.sh.match[1]}="'"${.sh.match[2]}"'
        else
            echo "ERROR: Ignoring invalid argument $__text\n" 1>&2
        fi
    done
    while IFS= read -r __text ; do
        while [[ $__text =~ ^([^@]*)(@\[[a-zA-Z_][a-zA-Z_0-9]*\]|@\[\'[^\']*\'\]|@\[@\]|@)(.*) ]] ; do
            __before="${.sh.match[1]}"
            __during="${.sh.match[2]}"
            __after="${.sh.match[3]}"
            print -n -- "$__before"
            if [[ "$__during" =~ ^@\[\'(.*)\'\]$ ]] ; then
                print -n -- "${.sh.match[1]}"
            elif [[ "$__during" == '@[@]' ]] ; then
                print -n @
            elif [[ "$__during" =~ ^@\[([a-zA-Z_][a-zA-Z_0-9]*)\] ]] ; then
                eval 'print -n -- "$'"${.sh.match[1]}"'"'
            else
                print -n "$__during"
            fi
            if [[ "$__after" == "$__text" ]] ; then
                break
            fi
            __text="$__after"
        done
        print -- "$__text"
    done
}

function test_atparse {
    # Note that these cannot be typeset since they will be invisible
    # to atparse:
    testvar='[testvar]'
    var1='[var1]'
    var2='[var2]'
    cat<<\EOF | atparse var3='**'
Nothing special here. = @['Nothing special here.']
[testvar] = @[testvar]
[var1] [var2] = @[var1] @[var2]
** = @[var3]
@ = @[@] = @['@']
-n
 eval "export PE$c=\${PE$c:-0}" = @[' eval "export PE$c=\${PE$c:-0}"']
EOF
    echo "After block, \$var3 = \"$var3\" should be empty"
}
