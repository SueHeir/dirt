#!/usr/bin/env bash
# Check DIRT against explicit GRASS and SOIL source trees.  The temporary Cargo
# patch is deliberately not committed as a workspace configuration: normal
# DIRT builds remain pinned to their reviewed dependencies.
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
GRASS_DIR=${GRASS_DIR:-"$HOME/projects/grass"}
SOIL_DIR=${SOIL_DIR:-"$HOME/projects/soil"}

usage() {
  cat <<'EOF'
Usage: ci/ecosystem-head-check.sh [--grass PATH] [--soil PATH]

Checks this DIRT checkout with every GRASS and SOIL dependency patched to the
given source trees.  Defaults are ~/projects/grass and ~/projects/soil.
EOF
}

while (($#)); do
  case "$1" in
    --grass) GRASS_DIR=${2:?--grass requires a path}; shift 2 ;;
    --soil) SOIL_DIR=${2:?--soil requires a path}; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

for repo in "$ROOT" "$GRASS_DIR" "$SOIL_DIR"; do
  # A linked Git worktree has a .git *file*, while a normal checkout has a
  # directory, so test for existence rather than directory type.
  if [[ ! -f "$repo/Cargo.toml" || ! -e "$repo/.git" ]]; then
    echo "ECOSYSTEM_HEAD error: expected a Git Cargo workspace at $repo" >&2
    exit 2
  fi
done

describe_repo() {
  local name=$1 path=$2 require_clean=${3:-true} dirty
  dirty=$(git -C "$path" status --porcelain)
  if [[ -n "$dirty" && "$require_clean" == true ]]; then
    echo "ECOSYSTEM_HEAD error: $name checkout is dirty: $path" >&2
    echo "$dirty" >&2
    exit 2
  fi
  printf 'ECOSYSTEM_HEAD %s commit=%s remote=%s path=%s dirty=%s\n' \
    "$name" "$(git -C "$path" rev-parse HEAD)" \
    "$(git -C "$path" remote get-url origin 2>/dev/null || echo '(no origin)')" "$path" \
    "$([[ -n "$dirty" ]] && echo yes || echo no)"
}

describe_repo dirt "$ROOT" false
describe_repo grass "$GRASS_DIR"
describe_repo soil "$SOIL_DIR"

for crate in grass_app grass_derive grass_io grass_mpi grass_scheduler; do
  test -f "$GRASS_DIR/crates/$crate/Cargo.toml" || {
    echo "ECOSYSTEM_HEAD error: GRASS HEAD lacks crates/$crate" >&2; exit 2;
  }
done
for crate in soil_core soil_deform soil_derive soil_fixes soil_print soil_verlet; do
  test -f "$SOIL_DIR/crates/$crate/Cargo.toml" || {
    echo "ECOSYSTEM_HEAD error: SOIL HEAD lacks crates/$crate" >&2; exit 2;
  }
done

TMP=$(mktemp -d "${TMPDIR:-/tmp}/dirt-ecosystem-head.XXXXXX")
trap 'rm -rf "$TMP"' EXIT
CONFIG="$TMP/config.toml"
METADATA="$TMP/metadata.json"

{
  echo '[patch."https://github.com/SueHeir/grass.git"]'
  for crate in grass_app grass_derive grass_io grass_mpi grass_scheduler; do
    printf '%s = { path = "%s" }\n' "$crate" "$GRASS_DIR/crates/$crate"
  done
  echo
  echo '[patch."https://github.com/SueHeir/soil.git"]'
  for crate in soil_core soil_deform soil_derive soil_fixes soil_print soil_verlet; do
    printf '%s = { path = "%s" }\n' "$crate" "$SOIL_DIR/crates/$crate"
  done
} > "$CONFIG"

cd "$ROOT"
echo 'ECOSYSTEM_HEAD running: cargo metadata --format-version 1'
cargo --config "$CONFIG" metadata --format-version 1 > "$METADATA"

python3 - "$METADATA" "$GRASS_DIR" "$SOIL_DIR" <<'PY'
import json
import pathlib
import sys

metadata, grass, soil = map(pathlib.Path, sys.argv[1:])
packages = {package["name"]: pathlib.Path(package["manifest_path"]).resolve()
            for package in json.loads(metadata.read_text())["packages"]}
expected = {
    **{name: grass / "crates" / name / "Cargo.toml"
       for name in ("grass_app", "grass_derive", "grass_io", "grass_mpi", "grass_scheduler")},
    **{name: soil / "crates" / name / "Cargo.toml"
       for name in ("soil_core", "soil_deform", "soil_derive", "soil_fixes", "soil_print", "soil_verlet")},
}
for name, manifest in expected.items():
    actual = packages.get(name)
    if actual != manifest.resolve():
        raise SystemExit(f"ECOSYSTEM_HEAD error: {name} resolved to {actual}, expected {manifest.resolve()}")
    print(f"ECOSYSTEM_HEAD resolved {name}={actual}")
PY

echo 'ECOSYSTEM_HEAD running: cargo check --workspace --no-default-features --features precision-double'
cargo --config "$CONFIG" check --workspace --no-default-features --features precision-double
echo 'ECOSYSTEM_HEAD PASS: metadata and non-MPI precision-double check used the printed GRASS, SOIL, and DIRT commits.'
