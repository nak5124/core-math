#!/bin/bash

MAKE=make

FUN="${!#}"
ARGS=("${@:1:$#-1}")

MODES=()
for i in "${!ARGS[@]}"; do
    case "${ARGS[i]}" in
        --rnd*)
            MODES+=("${ARGS[i]}")
            unset 'ARGS[i]'
            ;;
    esac
done
if [[ "${#MODES[@]}" -eq 0 ]]; then
    MODES=("--rndn" "--rndz" "--rndu" "--rndd")
fi

FILE="$(echo src/*/*/"$FUN".c)"
ORIG_DIR="$(dirname "$FILE")"

if ! [ -d "$ORIG_DIR" ]; then
    echo "Could not find $FUN"
    exit 1
fi

TMP_DIR="$(mktemp -d -t core-math.XXXXXX)"

trap 'rm -rf "$TMP_DIR"' EXIT

DIR="$TMP_DIR/toto/$(basename "$ORIG_DIR")"

mkdir "$TMP_DIR/toto"
cp -a "$ORIG_DIR" "$ORIG_DIR/../support" "$TMP_DIR/toto"
cp -a "$ORIG_DIR/../../generic" "$TMP_DIR"

if [ -n "${ARGS[0]}" ]; then
    KIND="${ARGS[0]}"
    unset 'ARGS[0]'
else
    SIZE=${FILE#src/binary}
    SIZE=${SIZE%%/*}
    case "$SIZE" in
        32)
            KIND=--exhaustive
            ;;
        *)
            KIND=--worst
    esac
fi

case "$KIND" in
    --exhaustive)
        "$MAKE" --quiet -C "$DIR" clean
        "$MAKE" --quiet -C "$DIR" check_exhaustive
        for MODE in "${MODES[@]}"; do
            echo "Running exhaustive check in $MODE mode..."
            "$DIR/check_exhaustive" "$MODE"
        done
        ;;
    --worst)
        "$MAKE" --quiet -C "$DIR" clean
        "$MAKE" --quiet -C "$DIR" check_worst
        for MODE in "${MODES[@]}"; do
            echo "Running worst cases check in $MODE mode..."
            "$DIR/check_worst" "$MODE" < "${FILE%.c}.wc"
        done
        ;;
    --special)
        "$MAKE" --quiet -C "$DIR" clean
        "$MAKE" --quiet -C "$DIR" check_special
        for MODE in "${MODES[@]}"; do
            echo "Running special checks in $MODE mode..."
            "$DIR/check_special" "$MODE"
        done
        ;;
    --exact)
        "$MAKE" --quiet -C "$DIR" clean
        "$MAKE" --quiet -C "$DIR" check_exact
        for MODE in "${MODES[@]}"; do
            echo "Running exact checks in $MODE mode..."
            "$DIR/check_exact" "$MODE"
        done
        ;;
    *)
        echo "Unrecognized command"
        exit 1
esac
