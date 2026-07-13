#!/usr/bin/env bash
cd /nfsshare/users/P126156127/workspace/gastric_cancer
dl() {
  url="${1%%$'\t'*}"; out="${1##*$'\t'}"
  # skip if exists and non-empty
  [ -s "$out" ] && return 0
  curl -s -f -C - --max-time 300 --retry 3 --retry-delay 5 -o "$out" "$url" || { echo "FAIL $url" >> logs_real/dl_fail.txt; rm -f "$out"; }
}
export -f dl
mkdir -p data/microbiome/raw logs_real
: > logs_real/dl_fail.txt
# tab-safe read via mapfile
mapfile -t LINES < data/microbiome/reprocess/dl_list.txt
printf '%s\n' "${LINES[@]}" | xargs -d '\n' -P 10 -I{} bash -c 'dl "$@"' _ {}
echo "DONE. present: $(ls data/microbiome/raw/*.bz2 2>/dev/null | wc -l)/1888  failures: $(wc -l < logs_real/dl_fail.txt)"
