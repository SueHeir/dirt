#!/usr/bin/env bash
# Verify the *negative* evidence boundary for the retired SPH repose study.
#
# This is intentionally not a calibration test: DIRT no longer contains the
# SPH solver or its angle-of-repose executable.  It makes the two independently
# checkable premises of the retirement note reproducible instead.
set -euo pipefail

repo_root=$(git rev-parse --show-toplevel)
source_ref=f7fe1a4
soil_sph_ref=b4997642678bf8072baa4d98be60429a4dfc59a9
soil_sph_dir=${SOIL_SPH_DIR:-}
online=false

usage() {
    cat <<'EOF'
Usage: ci/verify-retired-sph-repose.sh [--soil-sph DIR] [--online]

Checks that the SPH calibration surface was removed from DIRT and, when an
independent dev_soil_sph checkout is supplied, that the pinned snapshot has no
replacement angle-of-repose executable.  --online also verifies the cited
Crossref record's bibliographic identity.  Passing this script is not a
calibration, a material/protocol comparison, or an acceptance-gate pass.
EOF
}

while (($#)); do
    case $1 in
        --soil-sph)
            soil_sph_dir=${2:?--soil-sph needs a directory}
            shift 2
            ;;
        --online) online=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

cd "$repo_root"
git cat-file -e "${source_ref}^{commit}"
if git ls-tree -r --name-only HEAD -- examples crates | grep -q '^examples/SPH_glass_sphere_calibration/'; then
    echo 'FAIL: retired SPH calibration surface is present in this revision' >&2
    exit 1
fi

removed=$(git diff --numstat "${source_ref}^1" "$source_ref" -- \
    examples/SPH_glass_sphere_calibration Cargo.toml | awk '{add += $1; del += $2} END {print del}')
if [[ ${removed:-0} -lt 1 ]]; then
    echo 'FAIL: removal commit has no deleted SPH calibration content' >&2
    exit 1
fi
printf 'DIRT retirement boundary verified (%s SPH lines removed at %s).\n' "$removed" "$source_ref"

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
            ([.message.author[] | "\(.given) \(.family)"] | join("; "))] | @tsv')
    expected=$'ok\t10.1016/s0378-4371(99)00183-1\tRolling friction in the dynamic simulation of sandpile formation\tjournal-article\tY.C. Zhou; B.D. Wright; R.Y. Yang; B.H. Xu; A.B. Yu'
    if [[ $record != "$expected" ]]; then
        echo "FAIL: unexpected Crossref identity: $record" >&2
        exit 1
    fi
    echo 'Crossref identity verified; it remains bibliographic metadata only.'
fi

echo 'PASS: retirement facts verified; no calibration claim has been established.'
