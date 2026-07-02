#!/usr/bin/env bash

set -uo pipefail

workdir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$workdir"

input_file="input.nml"
executable="./skvp_AtomDiatom"
bspline_values=(80, 100, 120, 140  )

if [[ ! -f "$input_file" ]]; then
    echo "Error: $workdir/$input_file does not exist." >&2
    exit 1
fi

input_backup="$(mktemp "${TMPDIR:-/tmp}/input.nml.BMKP_42.XXXXXX")"
cp "$input_file" "$input_backup"

restore_input() {
    cp "$input_backup" "$input_file"
    rm -f "$input_backup"
}
trap restore_input EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

echo "Building SKVP executable..."
if ! make code; then
    echo "Error: make code failed. No calculations were started." >&2
    exit 1
fi

if [[ ! -x "$executable" ]]; then
    echo "Error: $workdir/skvp_AtomDiatom was not created." >&2
    exit 1
fi

set_basis() {
    local bspline_count="$1"
    local temporary_input
    temporary_input="$(mktemp "${TMPDIR:-/tmp}/input.nml.BMKP_42.edit.XXXXXX")"

    awk -v bspline_count="$bspline_count" '
        /^[[:space:]]*pbasst\(1\)[[:space:]]*=/ {
            sub(/\047bspl_fbr\047[[:space:]]+[0-9]+/,
                "\047bspl_fbr\047   " bspline_count)
        }
        /^[[:space:]]*pbasst\(2\)[[:space:]]*=/ {
            sub(/\047asleg_fbr\047[[:space:]]+[0-9]+[[:space:]]+[0-9]+/,
                "\047asleg_fbr\047    4  2")
        }
        /^[[:space:]]*pbasst\(3\)[[:space:]]*=/ {
            sub(/\047asleg_fbr\047[[:space:]]+[0-9]+[[:space:]]+[0-9]+/,
                "\047asleg_fbr\047    4  2")
        }
        { print }
    ' "$input_file" > "$temporary_input"

    mv "$temporary_input" "$input_file"
}

for bspline_count in "${bspline_values[@]}"; do
    log_file="BMKP_${bspline_count}_42.log"

    echo
    echo "============================================================"
    echo "Starting BMKP: B-spline=${bspline_count}, basis=(4,2,4,2)"
    echo "Log: $workdir/$log_file"
    echo "============================================================"

    set_basis "$bspline_count"

    # These files depend on the radial basis/grid. Rebuild them for every run.
    rm -f matrices.bin previous_calc.dat A_cache.dat

    if "$executable" 2>&1 | tee "$log_file"; then
        printf '\nBMKP %s 42\n' "$bspline_count" >> proba.dat
        echo "Completed BMKP ${bspline_count} 42"
        echo "Marker appended to $workdir/proba.dat"
    else
        run_status=${PIPESTATUS[0]}
        echo "Error: BMKP ${bspline_count} 42 failed with status ${run_status}." \
            | tee -a "$log_file" >&2
        echo "No completion marker was appended. Stopping the sequence." >&2
        exit "$run_status"
    fi
done

echo
echo "All BMKP B-spline tests completed successfully."
