# cols and dhead
## For some reason, vim doesn't understand the apostrophes used in this part and highlights it badly
function cols(){
    head -1 $1 | tr '\t' '\n' | nl -ba
}

function dhead(){
    if [ $# -ne 1 ]
    then
        echo "Usage: dhead <filename>"
        return 1
    fi
    diff -y <(cols ${1}) <(sed '2!d' ${1} | cols)
}

# random alnum string
function randomstr(){
    # The length of the string we want
    strlen=32
    # The number of lines in /dev/urandom to use since the data 
    # can be huge in a non-interactive script and will take forever to process
    search_space=1
    if [ "$#" -gt "2" ]
    then
        echo "Usage: randomstr [<LENGTH=32>] [<SEARCH_SPACE=1>]" >&2
        return 1
    elif [ "$#" -eq "1" ]
    then
        if [ ${1} == "-h" ] || [ ${1} == "--help" ]
        then
            echo "Usage: randomstr [<LENGTH=32>] [<SEARCH_SPACE=1>]" >&2
            return 1
        else
            strlen=${1}
        fi
    elif [ "$#" -eq "2" ]
    then
        strlen=${1}
        search_space=${2}
    fi
    match_found=0
    # Iterate 5 times till we get a match
    for iter in {1..5}
    do
        retval=$( head -n ${search_space} /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w ${strlen} | head -n 1 )
        # TODO: There's some edge case here!!
        if [ "$(expr length ${retval})" -eq "${strlen}" ]
        then
            match_found=1
            break
        fi
    done
    if [ "${match_found}" -eq "1" ]
    then
        echo ${retval}
        return 0
    else
        echo "Could not find a string after 10 attempts. Try increasing 'SEARCH_SPACE' but remember you will have to explicitly specify 'LENGTH'" >&2
        return 1
    fi
}
