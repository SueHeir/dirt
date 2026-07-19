#!/usr/bin/env bash
# Verify the availability/provenance boundary for a retired SPH repose study.
# This is intentionally not a calibration test.
set -euo pipefail

repo_root=$(git rev-parse --show-toplevel)
source_ref=f7fe1a4
retired_case=examples/SPH_glass_sphere_calibration/03_angle_of_repose
# Git-object witness, not a scientific threshold.
expected_case_deleted_lines=3511
soil_sph_ref=b4997642678bf8072baa4d98be60429a4dfc59a9
soil_sph_dir=${SOIL_SPH_DIR:-}
online=false

usage() {
    cat <<'EOF'
Usage: ci/verify-retired-sph-repose.sh [--soil-sph DIR] [--online]

Checks DIRT's historical removal of the SPH repose case and, when an
independent dev_soil_sph checkout is supplied, its pinned snapshot for a
replacement study. --online checks only Crossref bibliographic identity.
Passing is not a calibration, material comparison, or acceptance-gate pass.
EOF
}

while (($#)); do
    case $1 in
        --soil-sph) soil_sph_dir=${2:?--soil-sph needs a directory}; shift 2 ;;
        --online) online=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

cd "$repo_root"
git cat-file -e "${source_ref}^{commit}"
if ! git merge-base --is-ancestor "$source_ref" HEAD; then
    echo "FAIL: current revision does not retain retirement commit $source_ref" >&2
    exit 1
fi
if git ls-tree -r --name-only HEAD -- "$retired_case" | grep -q .; then
    echo "FAIL: retired SPH case is present: $retired_case" >&2
    exit 1
fi
removed=$(git diff --numstat "${source_ref}^1" "$source_ref" -- "$retired_case" |
    awk '{add += $1; del += $2} END {print del}')
if [[ ${removed:-0} -ne $expected_case_deleted_lines ]]; then
    echo "FAIL: expected $expected_case_deleted_lines deleted lines at $source_ref, got ${removed:-0}" >&2
    exit 1
fi
printf 'DIRT retirement boundary verified (%s lines removed from %s at %s).\n' \
    "$removed" "$retired_case" "$source_ref"

if [[ -n $soil_sph_dir ]]; then
    git -C "$soil_sph_dir" cat-file -e "${soil_sph_ref}^{commit}"
    if git -C "$soil_sph_dir" ls-tree -r --name-only "$soil_sph_ref" -- examples |
        grep -Eiq 'repose|calibrat'; then
        echo 'FAIL: pinned dev_soil_sph snapshot unexpectedly exposes a replacement study' >&2
        exit 1
    fi
    if ! git -C "$soil_sph_dir" grep -qi 'mu_r.*0.0' "$soil_sph_ref" -- docs/dem-campaign.md; then
        echo 'FAIL: pinned dev_soil_sph campaign statement changed or is unavailable' >&2
        exit 1
    fi
    printf 'No replacement SPH repose executable found at dev_soil_sph %s.\n' "$soil_sph_ref"
fi

if $online; then
    command -v curl >/dev/null
    command -v jq >/dev/null
    record=$(curl -LfsS 'https://api.crossref.org/works/10.1016/S0378-4371(99)00183-1' |
        jq -r '[.status, .message.DOI, .message.title[0], .message.type,
                ([.message.author[] | "\(.given) \(.family)"] | join(";"))] | @tsv')
    expected=$'ok\t10.1016/s0378-4371(99)00183-1\tRolling friction in the dynamic simulation of sandpile formation\tjournal-article\tY.C. Zhou;B.D. Wright;R.Y. Yang;B.H. Xu;A.B. Yu'
    if [[ $record != "$expected" ]]; then
        echo "FAIL: unexpected Crossref identity: $record" >&2
        exit 1
    fi
    echo 'Crossref identity verified; it remains bibliographic metadata only.'
fi

echo 'PASS: retirement facts verified; no calibration claim has been established.'
