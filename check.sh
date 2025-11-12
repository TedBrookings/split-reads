#!/usr/bin/env bash
set -euo pipefail
cargo check
cargo clippy
cargo fmt
