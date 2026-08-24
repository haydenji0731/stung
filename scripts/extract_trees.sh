#!/usr/bin/env bash

set -euo pipefail

tmp=$(mktemp); trap 'rm -f "$tmp"' EXIT
cat "${1:-/dev/stdin}" > "$tmp"

extract() {
    awk -v hdr="$1" '
        index($0, hdr) == 1 { grab = 1; next }
        grab { buf = buf $0; if (index($0, ";")) exit }
        END { print buf }
    ' "$tmp" \
    | tr -d "[:space:]" \
    | sed "s/;.*/;/" > "$2"

    if [ ! -s "$2" ]; then
        echo "WARNING: no tree found for header \"$1\"" >&2
    else
        n=$(grep -oE '[A-Za-z0-9_]+[:#]' "$2" | wc -l)
        echo "wrote $2  (${n} tips)" >&2
    fi
}

extract "dS tree:"                  dS.nwk
extract "dN tree:"                  dN.nwk
extract "w ratios as node labels:"  omega.nwk