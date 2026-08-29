#!/usr/bin/env bash
set -euo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }
log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }

resolve_path() {
  local value="$1"
  [[ -z "$value" ]] && return 0
  if [[ "$value" = /* ]]; then printf '%s\n' "$value"; else printf '%s/%s\n' "$GMRS_ROOT" "$value"; fi
}

manifest_for() {
  case "$1" in
    gwas) resolve_path "$GWAS_MANIFEST" ;;
    eqtl) resolve_path "$EQTL_MANIFEST" ;;
    pqtl) resolve_path "$PQTL_MANIFEST" ;;
    *) die "unknown layer '$1'" ;;
  esac
}

layers_for() {
  if [[ "$1" == all ]]; then printf '%s\n' gwas eqtl pqtl; else printf '%s\n' "$1"; fi
}

weight_dir() { local sub; case "$1" in gwas) sub="${GWAS_WEIGHT_SUBDIR:-gwas}";; eqtl) sub="${EQTL_WEIGHT_SUBDIR:-eqtl}";; pqtl) sub="${PQTL_WEIGHT_SUBDIR:-pqtl}";; esac; printf '%s/%s\n' "$(resolve_path "${WEIGHT_ROOT:-results/weights}")" "$sub"; }

