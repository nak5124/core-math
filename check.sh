#!/bin/bash

MAKE=make

KIND="$1"
shift

FUN="${!#}"
ARGS=("${@:1:$#-1}")

FILE="$(echo src/*/*/"$FUN".c)"
ORIG_DIR="$(dirname "$FILE")"

if ! [ -d "$ORIG_DIR" ]; then
    echo "Could not find $FUN"
    exit 1
fi

TMP_DIR="$(mktemp -d --tmpdir core-math.XXXXXX)"

trap 'rm -rf "$TMP_DIR"' EXIT

DIR="$TMP_DIR/$(basename "$ORIG_DIR")"

cp -a "$ORIG_DIR" "$ORIG_DIR/../support" "$TMP_DIR"

case "$KIND" in
    --exhaustive)
        "$MAKE" --quiet -C "$DIR" clean
        "$MAKE" --quiet -C "$DIR" check_exhaustive
        if [ "${#ARGS[@]}" -eq 0 ]; then
            MODES=("--rndn" "--rndz" "--rndu" "--rndd")
        else
            MODES=("${ARGS[@]}")
        fi
        for MODE in "${MODES[@]}"; do
            echo "Running exhaustive check in $MODE mode..."
            "$DIR/check_exhaustive" "$MODE"
        done
        ;;
    --worst)
        "$MAKE" --quiet -C "$DIR" clean
        "$MAKE" --quiet -C "$DIR" check_worst
        if [ "${#ARGS[@]}" -eq 0 ]; then
            MODES=("--rndn" "--rndz" "--rndu" "--rndd")
        else
            MODES=("${ARGS[@]}")
        fi
        for MODE in "${MODES[@]}"; do
            echo "Running worst cases check in $MODE mode..."
            "$DIR/check_worst" "$MODE" < "${FILE%.c}.wc"
        done
        ;;
    --special)
        "$MAKE" --quiet -C "$DIR" clean
        "$MAKE" --quiet -C "$DIR" check_special
        if [ "${#ARGS[@]}" -eq 0 ]; then
            MODES=("--rndn" "--rndz" "--rndu" "--rndd")
        else
            MODES=("${ARGS[@]}")
        fi
        for MODE in "${MODES[@]}"; do
            echo "Running special checks in $MODE mode..."
            "$DIR/check_special" "$MODE"
        done
        ;;
    *)
        echo "Unrecognized command"
        exit 1
esac
