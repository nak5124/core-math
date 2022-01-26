#!/bin/sh

MAKE=make

KIND="$1"
shift
FUN="$1"
shift

FILE="$(echo src/*/*/"$FUN".c)"
DIR="$(dirname "$FILE")"

if ! [ -d "$DIR" ]; then
    echo "Could not find $FUN"
    exit 1
fi

case "$KIND" in
    --exhaustive)
        "$MAKE" --quiet -C "$DIR" clean
        "$MAKE" --quiet -C "$DIR" check
        "$DIR/check" "$@"
        ;;
    *)
        echo "Unrecognized command"
        exit 1
esac
