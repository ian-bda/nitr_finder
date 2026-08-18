#!/bin/bash
set -euo pipefail
PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "$PIPE/lib.sh"
nitr_init

echo "=== tools on PATH ==="
for cmd in makeblastdb tblastn blastp hmmscan bedtools python3; do
  if command -v "$cmd" >/dev/null 2>&1; then
    printf "  OK  %-12s %s\n" "$cmd" "$(command -v "$cmd")"
  else
    printf "  MISS %-12s\n" "$cmd"
  fi
done

echo "=== optional tools ==="
for cmd in "$SAMTOOLS" "$MINIPROT" "${STAR:-}" "${STRINGTIE:-}" "${INTERPROSCAN:-}"; do
  if [[ -z "$cmd" ]]; then
    continue
  fi
  if [[ -x "$cmd" ]] || command -v "$cmd" >/dev/null 2>&1; then
    printf "  OK  %s\n" "$cmd"
  else
    printf "  MISS %s\n" "$cmd"
  fi
done

echo "=== Python packages ==="
python3 - <<'PY'
import importlib
ok = True
for m in ("Bio", "pyfaidx"):
    try:
        importlib.import_module(m)
        print(f"  OK  {m}")
    except ImportError:
        print(f"  MISS {m}  (pip install -r requirements.txt)")
        ok = False
raise SystemExit(0 if ok else 1)
PY

echo "=== required files ==="
nitr_need "$GENOME" "$I_QUERY"
if [[ ! -s "$SMART_HMM" ]]; then
  echo "  MISS SMART HMM: $SMART_HMM" >&2
  exit 1
fi
echo "  OK  genome, I-query, SMART HMM"
echo "Done: $(date)"
