#!/bin/bash

TREE_FILE=$1
SOURCE_PHY=$2
OUT_FILE=$3


taxa_list=$(tail -n +2 "$TREE_FILE" | tr -d '();' | tr ',' '\n' | grep -v '^$')

head -n 1 "$SOURCE_PHY" > "$OUT_FILE"

for taxon in $taxa_list; do
    grep -w "^$taxon" "$SOURCE_PHY" >> "$OUT_FILE"
done