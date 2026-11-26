#!/usr/bin/env bash
set -euo pipefail

if [[ -n "${DIRT_BUILD_ENV:-}" ]]; then
  # Optional site-local compiler/MPI environment. A normal configured shell
  # needs no extra setup.
  source "$DIRT_BUILD_ENV"
fi
cargo build --release --example fake_felted_tension
