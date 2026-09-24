#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
analysis_script="$repo_root/scripts/convergence_rate_analysis.py"
input_dir="$repo_root/solns/mms_adaptive_temperature_error"

mapfile -d '' org_files < <(
    find "$input_dir" -type f \
        \( -name 'error-global-q1.org' \
        -o -name 'error-global-q2.org' \
        -o -name 'error-global-q3.org' \) \
        -print0 | sort -z
)

if ((${#org_files[@]} == 0)); then
    echo "No matching error-global-q1/q2/q3.org files found under $input_dir" >&2
    exit 1
fi

for org_file in "${org_files[@]}"; do
    python3 "$analysis_script" "$org_file"
done
